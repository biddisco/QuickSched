/*******************************************************************************
 * This file is part of QuickSched.
 * Coypright (c) 2014 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
* *****************************************************************************/

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include <fenv.h>

/* Local includes. */
#include "quicksched.h"

/* Some local constants. */
#define cell_pool_grow 1000
#define cell_maxparts 16
#define task_limit 5000
#define const_G 6.6738e-8
#define dist_min 0.5 /* Used fpr legacy walk only */
#define iact_pair_direct iact_pair_direct_sorted

#define ICHECK -1

/** Data structure for the particles. */
struct part {
  double x[3];
  // union {
  float a[3];
  float a_legacy[3];
  float a_exact[3];
  //};
  float mass;
  int id;
} __attribute__((aligned(32)));

/** Data structure for the sorted particle positions. */
struct index {
  int ind;
  float d;
};

/** Data structure for the BH tree cell. */
struct cell {
  double loc[3];
  double h;

  int count;
  unsigned short int split, sorted;
  struct part *parts;

  struct cell *firstchild; /* Next node if opening */
  struct cell *sibling;    /* Next node */

  /* We keep both CoMs and masses to make sure the comp_com calculation is
   * correct (use an anonymous union to keep variable names compact).  */
  union {

    /* Information for the legacy walk */
    struct {
      double com[3];
      float mass;
    } legacy;

    /* Information for the QuickShed walk */
    struct {
      double com[3];
      float mass;
    } new;
  };

  int res, com_tid;
  struct index *indices;

} __attribute__((aligned(128)));

/** Task types. */
enum task_type {
  task_type_self = 0,
  task_type_pair,
  task_type_pair_pc,
  task_type_pair_direct,
  task_type_com,
  task_type_count
};

/** Per-type timers. */
ticks task_timers[task_type_count];

/** Global variable for the pool of allocated cells. */
struct cell *cell_pool = NULL;

/** Constants for the sorting axes. */
const float axis_shift[13 * 3] = {
  5.773502691896258e-01, 5.773502691896258e-01, 5.773502691896258e-01,
  7.071067811865475e-01, 7.071067811865475e-01, 0.0, 5.773502691896258e-01,
  5.773502691896258e-01, -5.773502691896258e-01, 7.071067811865475e-01, 0.0,
  7.071067811865475e-01, 1.0, 0.0, 0.0, 7.071067811865475e-01, 0.0,
  -7.071067811865475e-01, 5.773502691896258e-01, -5.773502691896258e-01,
  5.773502691896258e-01, 7.071067811865475e-01, -7.071067811865475e-01, 0.0,
  5.773502691896258e-01, -5.773502691896258e-01, -5.773502691896258e-01, 0.0,
  7.071067811865475e-01, 7.071067811865475e-01, 0.0, 1.0, 0.0, 0.0,
  7.071067811865475e-01, -7.071067811865475e-01, 0.0, 0.0, 1.0,
};
const char axis_flip[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* Map shift vector to sortlist. */
const int axis_sid[27] = {
  /* ( -1 , -1 , -1 ) */ 0,
  /* ( -1 , -1 ,  0 ) */ 1,
  /* ( -1 , -1 ,  1 ) */ 2,
  /* ( -1 ,  0 , -1 ) */ 3,
  /* ( -1 ,  0 ,  0 ) */ 4,
  /* ( -1 ,  0 ,  1 ) */ 5,
  /* ( -1 ,  1 , -1 ) */ 6,
  /* ( -1 ,  1 ,  0 ) */ 7,
  /* ( -1 ,  1 ,  1 ) */ 8,
  /* (  0 , -1 , -1 ) */ 9,
  /* (  0 , -1 ,  0 ) */ 10,
  /* (  0 , -1 ,  1 ) */ 11,
  /* (  0 ,  0 , -1 ) */ 12,
  /* (  0 ,  0 ,  0 ) */ 0,
  /* (  0 ,  0 ,  1 ) */ 12,
  /* (  0 ,  1 , -1 ) */ 11,
  /* (  0 ,  1 ,  0 ) */ 10,
  /* (  0 ,  1 ,  1 ) */ 9,
  /* (  1 , -1 , -1 ) */ 8,
  /* (  1 , -1 ,  0 ) */ 7,
  /* (  1 , -1 ,  1 ) */ 6,
  /* (  1 ,  0 , -1 ) */ 5,
  /* (  1 ,  0 ,  0 ) */ 4,
  /* (  1 ,  0 ,  1 ) */ 3,
  /* (  1 ,  1 , -1 ) */ 2,
  /* (  1 ,  1 ,  0 ) */ 1,
  /* (  1 ,  1 ,  1 ) */ 0
};

/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The indices
 * @param N The number of entries.
 */
void indices_sort(struct index *sort, int N) {

  struct {
    short int lo, hi;
  } qstack[10];
  int qpos, i, j, lo, hi, imin;
  struct index temp;
  float pivot;

  /* Sort parts in cell_i in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    if (hi - lo < 15) {
      for (i = lo; i < hi; i++) {
        imin = i;
        for (j = i + 1; j <= hi; j++)
          if (sort[j].d < sort[imin].d) imin = j;
        if (imin != i) {
          temp = sort[imin];
          sort[imin] = sort[i];
          sort[i] = temp;
        }
      }
    } else {
      pivot = sort[(lo + hi) / 2].d;
      i = lo;
      j = hi;
      while (i <= j) {
        while (sort[i].d < pivot) i++;
        while (sort[j].d > pivot) j--;
        if (i <= j) {
          if (i < j) {
            temp = sort[i];
            sort[i] = sort[j];
            sort[j] = temp;
          }
          i += 1;
          j -= 1;
        }
      }
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }
}

/**
 * @brief Sort the particles in the given cell along the given axis.
 *
 * @param c The #cell.
 * @param axis The normalized axis along which to sort.
 * @param aid The axis ID at which to store the indices.
 */
void cell_sort(struct cell *c, float *axis, int aid) {

  /* Has the indices array even been allocated? */
  if (c->indices == NULL) {
    if ((c->indices = (struct index *)malloc((c->count + 1) * 13 *
                                             sizeof(struct index))) ==
        NULL)
      error("Failed to allocate cell sorting indices.");
  } else if (c->sorted & (1 << aid)) {
    return;
  }

  /* If this cell has been split, merge from the progeny. */
  if (c->split) {

    /* Heap of pointers to the progeny's indices. */
    struct {
      struct index *index;
      int offset;
      } temp, pindex[8];
    int k = 0;

    /* First, make sure all the progeny have been sorted, and get the pointers
       to their first entries. */
    for (struct cell *cp = c->firstchild; cp != c->sibling; cp = cp->sibling) {
      if (!(cp->sorted & (1 << aid))) cell_sort(cp, axis, aid);
      pindex[k].index = &cp->indices[(cp->count + 1) * aid];
      pindex[k].offset = cp->parts - c->parts;
      k++;
    }
    
    /* Fill any remaining gaps with the sentinel. */
    c->indices[(c->count + 1) * aid + c->count].d = FLT_MAX;
    for (; k < 8; k++)
      pindex[k].index = &c->indices[(c->count + 1) * aid + c->count];

    /* Heapify the pindices. */
    for (k = 3; k >= 0; k--) {
      int j = k;
      while (j < 4) {
        int jj = 2*j + 1;
        if (jj + 1 < 8 && pindex[jj + 1].index->d < pindex[jj].index->d) jj = jj + 1;
        if (pindex[jj].index->d < pindex[j].index->d) {
          temp = pindex[jj];
          pindex[jj] = pindex[j];
          pindex[j] = temp;
          j = jj;
        } else
          break;
      }
    }

    /* Copy the indices into the local index in increasing order. */
    struct index *dest = &c->indices[(c->count + 1) * aid];
    for (int k = 0; k < c->count; k++) {

      /* Copy the top of the heap to the destination. */
      *dest = *pindex[0].index;
      dest->ind += pindex[0].offset;

      /* Increase the pointers. */
      dest++;
      pindex[0].index++;

      /* Fix the heap. */
      int j = 0;
      while (j < 4) {
        int jj = 2*j + 1;
        if (jj + 1 < 8 && pindex[jj + 1].index->d < pindex[jj].index->d) jj = jj + 1;
        if (pindex[jj].index->d < pindex[j].index->d) {
          temp = pindex[jj];
          pindex[jj] = pindex[j];
          pindex[j] = temp;
          j = jj;
        } else
          break;
      }
    }

    /* Add a sentinel at the end. */
    dest->ind = -1;
    dest->d = FLT_MAX;

    /* Otherwise, just sort the entries. */
  } else {

    /* Fill the indices. */
    struct index *dest = &c->indices[(c->count + 1) * aid];
    for (int k = 0; k < c->count; k++) {
      dest[k].ind = k;
      dest[k].d = c->parts[k].x[0] * axis[0] + c->parts[k].x[1] * axis[1] +
                  c->parts[k].x[2] * axis[2];
    }

    /* Sort the indices. */
    indices_sort(&c->indices[(c->count + 1) * aid], c->count);

    /* Set the sentinel on the last entry. */
    dest[c->count].ind = -1;
    dest[c->count].d = FLT_MAX;
  }

  /* Mark this cell as sorted in the given direction. */
  c->sorted |= (1 << aid);
  
  /* Verify the sort. */
  /* int offset = (c->count + 1) * aid;
  for (int k = 0; k < c->count; k++)
    if (c->indices[offset + k].d > c->indices[offset + k + 1].d)
      error( "Sorting failed." ); */
} /*  */

/**
 * @brief Get all the data needed for computing a sorted pair.
 *
 * @param ci Pointer to a pointer to the first cell.
 * @param cj Pointer to a pointer to the second cell.
 * @param ind_i Sorted indices of the cell @c ci.
 * @param ind_j Sorted indices of the cell @c cj.
 * @param corr Axis distance correction factor, i.e. the scaling applied
 *        to the distance along the axis.
 */
void get_axis(struct cell **ci, struct cell **cj, struct index **ind_i,
              struct index **ind_j, float *corr) {

  float dx[3], axis[3];
  int aid = 0;

  /* Get the cell pair separation and the axis index. */
  for (int k = 0; k < 3; k++) {
    dx[k] = (*cj)->loc[k] - (*ci)->loc[k];
    aid = 3 * aid + ((dx[k] < 0) ? 0 : (dx[k] > 0) ? 2 : 1);
  }

  /* Flip the cells? */
  if (axis_flip[aid]) {
    struct cell *temp = *ci;
    *ci = *cj;
    *cj = temp;
  }
  aid = axis_sid[aid];

  /* Copy the shift vector to the output parameter. */
  axis[0] = axis_shift[aid * 3 + 0];
  axis[1] = axis_shift[aid * 3 + 1];
  axis[2] = axis_shift[aid * 3 + 2];

  /* Make sure the cells are sorted. */
  if (!((*ci)->sorted & (1 << aid))) cell_sort(*ci, axis, aid);
  if (!((*cj)->sorted & (1 << aid))) cell_sort(*cj, axis, aid);

  /* Set the indices. */
  *ind_i = &(*ci)->indices[aid * ((*ci)->count + 1)];
  *ind_j = &(*cj)->indices[aid * ((*cj)->count + 1)];

  /* Compute the axis scaling correction. */
  *corr = (axis[0] * dx[0] + axis[1] * dx[1] + axis[2] * dx[2]) /
          sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

  /* message("dx=[%.3e,%.3e,%.3e], axis=[%.3e,%.3e,%.3e], corr=%.3e.",
    dx[0], dx[1], dx[2], axis[0], axis[1], axis[2], *corr); */

  /* Make sure the sorts are ok. */
  /* for (int k = 1; k < (*ci)->count; k++)
    if ((*ind_i)[k].d < (*ind_i)[k-1].d)
      error("Sorting failed.");
  for (int k = 1; k < (*cj)->count; k++)
    if ((*ind_j)[k].d < (*ind_j)[k-1].d)
      error("Sorting failed."); */
}

/**
 * @brief Get a #cell from the pool.
 */
struct cell *cell_get() {

  struct cell *res;
  int k;

  /* Allocate a new batch? */
  if (cell_pool == NULL) {

    /* Allocate the cell array. */
    if ((cell_pool =
             (struct cell *)calloc(cell_pool_grow, sizeof(struct cell))) ==
        NULL)
      error("Failed to allocate fresh batch of cells.");

    /* Link them up via their progeny pointers. */
    for (k = 1; k < cell_pool_grow; k++)
      cell_pool[k - 1].firstchild = &cell_pool[k];
  }

  /* Pick a cell off the pool. */
  res = cell_pool;
  cell_pool = cell_pool->firstchild;

  /* Clean up a few things. */
  res->res = qsched_res_none;
  res->firstchild = 0;

  /* Return the cell. */
  return res;
}

/**
 * @brief Sort the parts into eight bins along the given pivots and
 *        fill the multipoles. Also adds the hierarchical resources
 *        to the sched.
 *
 * @param c The #cell to be split.
 * @param N The total number of parts.
 * @param s The #sched to store the resources.
 */
void cell_split(struct cell *c, struct qsched *s) {

  int i, j, k, kk, count = c->count;
  struct part temp, *parts = c->parts;
  struct cell *cp;
  int left[8], right[8];
  double pivot[3];
  static struct cell *root = NULL;
  struct cell *progenitors[8];

  /* Set the root cell. */
  if (root == NULL) {
    root = c;
    c->sibling = 0;
  }

  /* Add a resource for this cell if it doesn't have one yet. */
  if (c->res == qsched_res_none)
    c->res = qsched_addres(s, qsched_owner_none, qsched_res_none);

  /* Attach a center-of-mass task to the cell. */
  if (count > 0)
    c->com_tid = qsched_addtask(s, task_type_com, task_flag_none, &c,
                                sizeof(struct cell *), 1);

  /* Does this cell need to be split? */
  if (count > cell_maxparts) {

    /* Mark this cell as split. */
    c->split = 1;

    /* Create the progeny. */
    for (k = 0; k < 8; k++) {
      progenitors[k] = cp = cell_get();
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->h = c->h * 0.5;
      cp->res = qsched_addres(s, qsched_owner_none, c->res);
      if (k & 4) cp->loc[0] += cp->h;
      if (k & 2) cp->loc[1] += cp->h;
      if (k & 1) cp->loc[2] += cp->h;
    }

    /* Init the pivots. */
    for (k = 0; k < 3; k++) pivot[k] = c->loc[k] + c->h * 0.5;

    /* Split along the x-axis. */
    i = 0;
    j = count - 1;
    while (i <= j) {
      while (i <= count - 1 && parts[i].x[0] < pivot[0]) i += 1;
      while (j >= 0 && parts[j].x[0] >= pivot[0]) j -= 1;
      if (i < j) {
        temp = parts[i];
        parts[i] = parts[j];
        parts[j] = temp;
      }
    }
    left[1] = i;
    right[1] = count - 1;
    left[0] = 0;
    right[0] = j;

    /* Split along the y axis, twice. */
    for (k = 1; k >= 0; k--) {
      i = left[k];
      j = right[k];
      while (i <= j) {
        while (i <= right[k] && parts[i].x[1] < pivot[1]) i += 1;
        while (j >= left[k] && parts[j].x[1] >= pivot[1]) j -= 1;
        if (i < j) {
          temp = parts[i];
          parts[i] = parts[j];
          parts[j] = temp;
        }
      }
      left[2 * k + 1] = i;
      right[2 * k + 1] = right[k];
      left[2 * k] = left[k];
      right[2 * k] = j;
    }

    /* Split along the z axis, four times. */
    for (k = 3; k >= 0; k--) {
      i = left[k];
      j = right[k];
      while (i <= j) {
        while (i <= right[k] && parts[i].x[2] < pivot[2]) i += 1;
        while (j >= left[k] && parts[j].x[2] >= pivot[2]) j -= 1;
        if (i < j) {
          temp = parts[i];
          parts[i] = parts[j];
          parts[j] = temp;
        }
      }
      left[2 * k + 1] = i;
      right[2 * k + 1] = right[k];
      left[2 * k] = left[k];
      right[2 * k] = j;
    }

    /* Store the counts and offsets. */
    for (k = 0; k < 8; k++) {
      progenitors[k]->count = right[k] - left[k] + 1;
      progenitors[k]->parts = &c->parts[left[k]];
    }

    /* Find the first non-empty progenitor */
    for (k = 0; k < 8; k++)
      if (progenitors[k]->count > 0) {
        c->firstchild = progenitors[k];
        // message( "first child= %d", k);
        break;
      }
    if (c->firstchild == NULL)
      error("Cell has been split but all progenitors have 0 particles");

    /* Prepare the pointers. */
    for (k = 0; k < 8; k++) {

      /* Find the next non-empty sibling */
      for (kk = k + 1; kk < 8; ++kk) {
        if (progenitors[kk]->count > 0) {
          progenitors[k]->sibling = progenitors[kk];
          break;
        }
      }

      /* No non-empty sibling ? Go back a level */
      if (kk == 8) progenitors[k]->sibling = c->sibling;
    }

    /* Recurse. */
    for (k = 0; k < 8; k++)
      if (progenitors[k]->count > 0) cell_split(progenitors[k], s);

    /* Link the COM tasks. */
    for (k = 0; k < 8; k++)
      if (progenitors[k]->count > 0)
        qsched_addunlock(s, progenitors[k]->com_tid, c->com_tid);

  } /* does the cell need to be split? */

  /* Set this cell's resources ownership. */
  qsched_res_own(s, c->res,
                 s->nr_queues * (c->parts - root->parts) / root->count);
}

/* -------------------------------------------------------------------------- */
/* New tree walk */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the center of mass of a given cell.
 *
 * @param c The #cell.
 */
void comp_com(struct cell *c) {

  int k, count = c->count;
  struct cell *cp;
  struct part *p, *parts = c->parts;
  double com[3] = {0.0, 0.0, 0.0}, mass = 0.0;

  if (c->split) {

    /* Loop over the projenitors and collect the multipole information. */
    for (cp = c->firstchild; cp != c->sibling; cp = cp->sibling) {
      float cp_mass = cp->new.mass;
      com[0] += cp->new.com[0] * cp_mass;
      com[1] += cp->new.com[1] * cp_mass;
      com[2] += cp->new.com[2] * cp_mass;
      mass += cp_mass;
    }

    /* Otherwise, collect the multipole from the particles. */
  } else {

    for (k = 0; k < count; k++) {
      p = &parts[k];
      float p_mass = p->mass;
      com[0] += p->x[0] * p_mass;
      com[1] += p->x[1] * p_mass;
      com[2] += p->x[2] * p_mass;
      mass += p_mass;
    }
  }

  /* Store the COM data, if any was collected. */
  if (mass > 0.0) {
    float imass = 1.0f / mass;
    c->new.com[0] = com[0] * imass;
    c->new.com[1] = com[1] * imass;
    c->new.com[2] = com[2] * imass;
    c->new.mass = mass;
  } else {
    c->new.com[0] = 0.0;
    c->new.com[1] = 0.0;
    c->new.com[2] = 0.0;
    c->new.mass = 0.0;
  }
}

/**
 * @brief Compute the interactions between all particles in a cell
 *        and the center of mass of another cell.
 *
 * @param ci The #cell containing the particles.
 * @param cj The #cell containing the center of mass.
 */
static inline void iact_pair_pc(struct cell *ci, struct cell *cj) {
  int j, k, count = ci->count;
  double com[3];
  float mcom, dx[3], r2, ir, w;
  struct part *parts = ci->parts;

  /* Early abort? */
  if (count == 0 || cj->count == 0) error("Empty cell!");

  /* Sanity check. */
  if (cj->new.mass == 0.0) {
    message("%e %e %e %d %p %d %p", cj->new.com[0], cj->new.com[1],
            cj->new.com[2], cj->count, cj, cj->split, cj->sibling);

    for (j = 0; j < cj->count; ++j)
      message("part %d mass=%e id=%d", j, cj->parts[j].mass, cj->parts[j].id);

    error("com does not seem to have been set.");
  }

  /* Init the com's data. */
  for (k = 0; k < 3; k++) com[k] = cj->new.com[k];
  mcom = cj->new.mass;

  /* Loop over every particle in ci. */
  for (j = 0; j < count; j++) {

    /* Compute the pairwise distance. */
    for (r2 = 0.0, k = 0; k < 3; k++) {
      dx[k] = com[k] - parts[j].x[k];
      r2 += dx[k] * dx[k];
    }

    /* Apply the gravitational acceleration. */
    ir = 1.0f / sqrtf(r2);
    w = mcom * const_G * ir * ir * ir;
    for (k = 0; k < 3; k++) parts[j].a[k] += w * dx[k];

#if ICHECK >= 0
    if (parts[j].id == ICHECK)
      printf("[NEW] Can interact with the monopole. x= %f %f %f m= %f h= %f\n",
             com[0], com[1], com[2], mcom, cj->h);
#endif

  } /* loop over every particle. */
}

/**
 * @brief Compute the interactions between all particles in a cell
 *        and all particles in the other cell (N^2 algorithm)
 *
 * @param ci The #cell containing the particles.
 * @param cj The #cell containing the other particles
 */
static inline void iact_pair_direct_unsorted(struct cell *ci, struct cell *cj) {

  int i, j, k;
  int count_i = ci->count, count_j = cj->count;
  struct part *parts_i = ci->parts, *parts_j = cj->parts;
  double xi[3];
  float dx[3], ai[3], mi, mj, r2, w, ir;

  /* Bad stuff will happen if cell sizes are different */
   if ( ci->h != cj->h )
     error("Non matching cell sizes !! h_i=%f h_j=%f\n", ci->h, cj->h);

  /* Loop over all particles in ci... */
  for (i = 0; i < count_i; i++) {

    /* Init the ith particle's data. */
    for (k = 0; k < 3; k++) {
      xi[k] = parts_i[i].x[k];
      ai[k] = 0.0f;
    }
    mi = parts_i[i].mass;

    /* Loop over every following particle. */
    for (j = 0; j < count_j; j++) {

      /* Compute the pairwise distance. */
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = xi[k] - parts_j[j].x[k];
        r2 += dx[k] * dx[k];
      }

      /* Apply the gravitational acceleration. */
      ir = 1.0f / sqrtf(r2);
      w = const_G * ir * ir * ir;
      mj = parts_j[j].mass;
      for (k = 0; k < 3; k++) {
        float wdx = w * dx[k];
        parts_j[j].a[k] += wdx * mi;
        ai[k] -= wdx * mj;
      }

#if ICHECK >= 0
      if (parts_i[i].id == ICHECK)
        printf("[NEW] Interaction with particle id= %d (pair i)\n",
               parts_j[j].id);

      if (parts_j[j].id == ICHECK)
        printf("[NEW] Interaction with particle id= %d (pair j) h_i= %f h_j= "
               "%f ci= %p cj= %p count_i= %d count_j= %d d_i= %d d_j= %d\n",
               parts_i[i].id, ci->h, cj->h, ci, cj, count_i, count_j, ci->res,
               cj->res);
#endif

    } /* loop over every other particle. */

    /* Store the accumulated acceleration on the ith part. */
    for (k = 0; k < 3; k++) parts_i[i].a[k] += ai[k];

  } /* loop over all particles. */
}

/**
 * @brief Compute the interactions between all particles in a cell
 *        and all particles in the other cell using osrted interactions
 *
 * @param ci The #cell containing the particles.
 * @param cj The #cell containing the other particles
 */
static inline void iact_pair_direct_sorted(struct cell *ci, struct cell *cj) {

  int i, j, k;
  int count_i, count_j;
  struct part *parts_i, *parts_j;
  double cih = ci->h, cjh = cj->h;
  double xi[3];
  float dx[3], ai[3], mi, mj, r2, w, ir;

   /* Bad stuff will happen if cell sizes are different */
   if ( ci->h != cj->h )
     error("Non matching cell sizes !! h_i=%f h_j=%f\n", ci->h, cj->h);

  /* Get the sorted indices and stuff. */
  struct index *ind_i, *ind_j;
  float corr;
  double com[3] = {0.0, 0.0, 0.0};
  float com_mass = 0.0;
  get_axis(&ci, &cj, &ind_i, &ind_j, &corr);
  count_i = ci->count;
  parts_i = ci->parts;
  cih = ci->h;
  count_j = cj->count;
  parts_j = cj->parts;
  cjh = cj->h;

  /* Distance along the axis as of which we will use a multipole. */
  float d_max =  cjh / dist_min / corr;

  /* Loop over all particles in ci... */
  for (i = count_i - 1; i >= 0; i--) {

    /* Get the sorted index. */
    int pid = ind_i[i].ind;
    float di = ind_i[i].d;

    /* Init the ith particle's data. */
    for (k = 0; k < 3; k++) {
      xi[k] = parts_i[pid].x[k];
      ai[k] = 0.0;
    }
    mi = parts_i[pid].mass;

    /* Loop over every following particle within d_max along the axis. */
    for (j = 0; j < count_j && (ind_j[j].d - di) < d_max; j++) {

      /* Get the sorted index. */
      int pjd = ind_j[j].ind;

      /* Compute the pairwise distance. */
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = xi[k] - parts_j[pjd].x[k];
        r2 += dx[k] * dx[k];
      }

      /* Apply the gravitational acceleration. */
      ir = 1.0f / sqrtf(r2);
      w = const_G * ir * ir * ir;
      mj = parts_j[pjd].mass;
      for (k = 0; k < 3; k++) {
        float wdx = w * dx[k];
        parts_j[pjd].a[k] += wdx * mi;
        ai[k] -= wdx * mj;
      }

#if ICHECK >= 0
      if (parts_i[pid].id == ICHECK)
        printf("[NEW] Interaction with particle id= %d (pair i)\n",
               parts_j[pjd].id);

      if (parts_j[j].id == ICHECK)
        printf("[NEW] Interaction with particle id= %d (pair j) h_i= %f h_j= "
               "%f ci= %p cj= %p count_i= %d count_j= %d d_i= %d d_j= %d\n",
               parts_i[pid].id, cih, cjh, ci, cj, count_i, count_j, ci->res,
               cj->res);
#endif

    } /* loop over every other particle. */

    /* Add any remaining particles to the COM. */
    for (int jj = j; jj < count_j; jj++) {
      int pjd = ind_j[jj].ind;
      mj = parts_j[pjd].mass;
      com[0] += mj * parts_j[pjd].x[0];
      com[1] += mj * parts_j[pjd].x[1];
      com[2] += mj * parts_j[pjd].x[2];
      com_mass += mj;
    }

    /* Shrink count_j to the latest valid particle. */
    count_j = j;

    /* Interact part_i with the center of mass. */
    if (com_mass > 0.0) {
      float icom_mass = 1.0f / com_mass;
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = xi[k] - com[k] * icom_mass;
        r2 += dx[k] * dx[k];
      }
      ir = 1.0f / sqrtf(r2);
      w = const_G * ir * ir * ir;
      for (k = 0; k < 3; k++) ai[k] -= w * dx[k] * com_mass;
    }

    /* Store the accumulated acceleration on the ith part. */
    for (k = 0; k < 3; k++) parts_i[pid].a[k] += ai[k];

  } /* loop over all particles in ci. */

  /* Loop over the particles in cj, catch the COM interactions. */
  count_j = cj->count;
  int last_i = 0;
  com[0] = 0.0;
  com[1] = 0.0;
  com[2] = 0.0;
  com_mass = 0.0f;
  d_max = cih / dist_min / corr;
  for (j = 0; j < count_j; j++) {

    /* Get the sorted index. */
    int pjd = ind_j[j].ind;
    float dj = ind_j[j].d;

    /* Fill the COM with any new particles. */
    for (i = last_i; i < count_i && (dj - ind_i[i].d) > d_max; i++) {
      int pid = ind_i[i].ind;
      mi = parts_i[pid].mass;
      com[0] += parts_i[pid].x[0] * mi;
      com[1] += parts_i[pid].x[1] * mi;
      com[2] += parts_i[pid].x[2] * mi;
      com_mass += mi;
    }

    /* Set the new last_i to the last particle checked. */
    last_i = i;

    /* Interact part_j with the COM. */
    if (com_mass > 0.0) {
      float icom_mass = 1.0f / com_mass;
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = com[k] * icom_mass - parts_j[pjd].x[k];
        r2 += dx[k] * dx[k];
      }
      ir = 1.0f / sqrtf(r2);
      w = const_G * ir * ir * ir;
      for (k = 0; k < 3; k++) parts_j[pjd].a[k] += w * dx[k] * com_mass;
    }
  }
}

/**
 * @brief Decides whether two cells use the direct summation interaction or the
* multipole interactions
 *
 * @param ci The #cell.
 * @param cj The other #cell.
 */
void iact_pair(struct cell *ci, struct cell *cj) {

  int k;
  int count_i = ci->count, count_j = cj->count;
  double center_i, center_j, dx[3];
  double min_dist, cih = ci->h, cjh = cj->h;
  struct cell *cp, *cps;

  /* Early abort? */
  if (count_i == 0 || count_j == 0) error("Empty cell !");

   /* Bad stuff will happen if cell sizes are different */
   if ( ci->h != cj->h )
     error("Non matching cell sizes !! h_i=%f h_j=%f\n", ci->h, cj->h);

  /* Sanity check */
  if (ci == cj)
    error("The impossible has happened: pair interaction between a cell and "
          "itself.");

  /* Distance between the cell centers */
  for (k = 0; k < 3; k++) {
    center_i = ci->loc[k] + 0.5 * cih;
    center_j = cj->loc[k] + 0.5 * cjh;
    dx[k] = fabs(center_i - center_j);
  }

  min_dist = cih + cjh;

  /* Are the cells NOT neighbours ? */
  if ((dx[0] > min_dist) || (dx[1] > min_dist) || (dx[2] > min_dist)) {

    iact_pair_pc(ci, cj);
    iact_pair_pc(cj, ci);

  } else {/* Cells are direct neighbours */

    /* Are both cells split ? */
    if (ci->split && cj->split) {

      /* Let's split both cells and build all possible pairs */
      for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling)
	for (cps = cj->firstchild; cps != cj->sibling; cps = cps->sibling)
	  iact_pair(cp, cps);

    } else {/* Ok one of the cells is not split */

      iact_pair_direct(ci, cj);
    }
  }
}

/**
 * @brief Compute the interactions between all particles in a cell.
 *
 * @param c The #cell.
 */
void iact_self(struct cell *c) {
  int i, j, k, count = c->count;
  double xi[3] = {0.0, 0.0, 0.0};
  float ai[3] = {0.0, 0.0, 0.0}, mi, mj, dx[3] = {0.0, 0.0, 0.0}, r2, ir, w;
  struct part *parts = c->parts;
  struct cell *cp, *cps;

  /* Early abort? */
  if (count == 0) error("Empty cell !");

  /* If the cell is split, interact each progeny with itself, and with
     each of its siblings. */
  if (c->split) {
    for (cp = c->firstchild; cp != c->sibling; cp = cp->sibling) {
      iact_self(cp);
      for (cps = cp->sibling; cps != c->sibling; cps = cps->sibling)
        iact_pair(cp, cps);
    }

    /* Otherwise, compute the interactions directly. */
  } else {

    /* Loop over all particles... */
    for (i = 0; i < count; i++) {

      /* Init the ith particle's data. */
      for (k = 0; k < 3; k++) {
        xi[k] = parts[i].x[k];
        ai[k] = 0.0;
      }
      mi = parts[i].mass;

      /* Loop over every following particle. */
      for (j = i + 1; j < count; j++) {

        /* Compute the pairwise distance. */
        for (r2 = 0.0, k = 0; k < 3; k++) {
          dx[k] = xi[k] - parts[j].x[k];
          r2 += dx[k] * dx[k];
        }

        /* Apply the gravitational acceleration. */
        ir = 1.0f / sqrtf(r2);
        w = const_G * ir * ir * ir;
        mj = parts[j].mass;
        for (k = 0; k < 3; k++) {
          float wdx = w * dx[k];
          parts[j].a[k] += wdx * mi;
          ai[k] -= wdx * mj;
        }

#if ICHECK >= 0
        if (parts[i].id == ICHECK)
          message("[NEW] Interaction with particle id= %d (self i)",
                  parts[j].id);

        if (parts[j].id == ICHECK)
          message("[NEW] Interaction with particle id= %d (self j)",
                  parts[i].id);
#endif

      } /* loop over every other particle. */

      /* Store the accumulated acceleration on the ith part. */
      for (k = 0; k < 3; k++) parts[i].a[k] += ai[k];

    } /* loop over all particles. */

  } /* otherwise, compute interactions directly. */
}

/**
 * @brief Create the tasks for the cell pair/self.
 *
 * @param s The #sched in which to create the tasks.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void create_tasks(struct qsched *s, struct cell *ci, struct cell *cj) {

  int k;
  qsched_task_t tid;
  double center_i, center_j, dx[3];
  double min_dist, cih, cjh;
  struct cell *data[2], *cp, *cps;

  /* If either cell is empty, stop. */
  if (ci->count == 0 || (cj != NULL && cj->count == 0)) error("Empty cell !");

  /* Single cell? */
  if (cj == NULL) {

    /* Is this cell split? */
    if (ci->split && ci->count > task_limit) {

      /* Loop over each of this cell's progeny. */
      for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling) {

        /* Make self-interaction task. */
        create_tasks(s, cp, NULL);

        /* Make all pair-interaction tasks. */
        for (cps = cp->sibling; cps != ci->sibling; cps = cps->sibling)
          create_tasks(s, cp, cps);
      }

      /* Otherwise, add a self-interaction task. */
    } else {

      /* Set the data. */
      data[0] = ci;
      data[1] = NULL;

      /* Create the task. */
      tid =
          qsched_addtask(s, task_type_self, task_flag_none, data,
                         sizeof(struct cell *) * 2, ci->count * ci->count / 2);

      /* Add the resource (i.e. the cell) to the new task. */
      qsched_addlock(s, tid, ci->res);

      /* If this call might recurse, add a dependency on the cell's COM
         task. */
      if (ci->split) qsched_addunlock(s, ci->com_tid, tid);
    }

    /* Otherwise, it's a pair. */
  } else {

    cih = ci->h;
    cjh = cj->h;

    /* Distance between the cell centers */
    for (k = 0; k < 3; k++) {
      center_i = ci->loc[k] + 0.5 * cih;
      center_j = cj->loc[k] + 0.5 * cjh;
      dx[k] = fabs(center_i - center_j);
    }

    min_dist = cih + cjh;

    /* Are the cells NOT neighbours ? */
    if ((dx[0] > min_dist) || (dx[1] > min_dist) || (dx[2] > min_dist)) {

      data[0] = ci;
      data[1] = cj;
      tid = qsched_addtask(s, task_type_pair_pc, task_flag_none, data,
                           sizeof(struct cell *) * 2, ci->count);
      qsched_addlock(s, tid, ci->res);
      qsched_addunlock(s, cj->com_tid, tid);

      data[0] = cj;
      data[1] = ci;
      tid = qsched_addtask(s, task_type_pair_pc, task_flag_none, data,
                           sizeof(struct cell *) * 2, cj->count);
      qsched_addlock(s, tid, cj->res);
      qsched_addunlock(s, ci->com_tid, tid);

    } else {/* Cells are direct neighbours */

      /* Are both cells split ? */
      if (ci->split && cj->split) {

        if (ci->count > task_limit && cj->count > task_limit) {

	  /* Let's split both cells and build all possible pairs */
	  for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling)
	    for (cps = cj->firstchild; cps != cj->sibling; cps = cps->sibling)
	      create_tasks(s, cp, cps);

        } else {

	  /* Let's split both cells and build all possible pairs */
	  for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling) {
	    for (cps = cj->firstchild; cps != cj->sibling; cps = cps->sibling) {
	
	      /* Set the data */
	      data[0] = cp;
	      data[1] = cps;

              /* Create the task. */
              tid = qsched_addtask(s, task_type_pair, task_flag_none, data,
                                   sizeof(struct cell *) * 2,
                                   cp->count * cps->count);

              /* Add the resources. */
              qsched_addlock(s, tid, cp->res);
              qsched_addlock(s, tid, cps->res);

              /* Depend on the COMs in case this task recurses. */
              if (cp->split || cps->split) {
                qsched_addunlock(s, cp->com_tid, tid);
                qsched_addunlock(s, cps->com_tid, tid);
	      }
	    }
	  }
        }
      } else {/* Ok one of the cells is not split */

        /* Set the data. */
        data[0] = ci;
        data[1] = cj;

        /* Create the task. */
        tid = qsched_addtask(s, task_type_pair_direct, task_flag_none, data,
                             sizeof(struct cell *) * 2, ci->count * cj->count);

        /* Add the resources. */
        qsched_addlock(s, tid, ci->res);
        qsched_addlock(s, tid, cj->res);
      }
    } /* Cells are direct neighbours */
  }   /* Otherwise it's a pair */
}

/* -------------------------------------------------------------------------- */
/* Legacy tree walk */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the center of mass of a given cell recursively.
 *
 * @param c The #cell.
 */
void legacy_comp_com(struct cell *c, int *countCoMs) {

  int k, count = c->count;
  struct cell *cp;
  struct part *p, *parts = c->parts;
  double com[3] = {0.0, 0.0, 0.0}, mass = 0.0;

  ++(*countCoMs);

  /* Is the cell split? */
  if (c->split) {

    /* Loop over the progeny. */
    for (cp = c->firstchild; cp != c->sibling; cp = cp->sibling) {
      /* Recurse */
      legacy_comp_com(cp, countCoMs);

      /* Collect multipole information */
      float cp_mass = cp->legacy.mass;
      com[0] += cp->legacy.com[0] * cp_mass;
      com[1] += cp->legacy.com[1] * cp_mass;
      com[2] += cp->legacy.com[2] * cp_mass;
      mass += cp_mass;
    }

    /* Otherwise, collect the multipole from local data. */
  } else {

    for (k = 0; k < count; k++) {
      p = &parts[k];
      float p_mass = p->mass;
      com[0] += p->x[0] * p_mass;
      com[1] += p->x[1] * p_mass;
      com[2] += p->x[2] * p_mass;
      mass += p_mass;
    }
  }

  /* Finish multipole calculation */
  if (mass > 0.0) {
    float imass = 1.0f / mass;
    c->legacy.com[0] = com[0] * imass;
    c->legacy.com[1] = com[1] * imass;
    c->legacy.com[2] = com[2] * imass;
    c->legacy.mass = mass;
  } else {
    c->legacy.com[0] = 0.0;
    c->legacy.com[1] = 0.0;
    c->legacy.com[2] = 0.0;
    c->legacy.mass = 0.0;
  }
}

/**
 * @brief Interacts a particle with a cell recursively using the original B-H
 * tree walk procedure
 *
 * @param parts The array of particles
 * @param i The particle of interest
 * @param root The root of the tree under which we will search.
 * @param monitor If set to @c parts[i].id, will produce debug output when
 *        ICHECK is set.
 * @param cell The cell the particle interacts with
 */
void legacy_interact(struct part *parts, int i, struct cell *root, int monitor,
                     int *countMultipoles, int *countPairs) {

  int j, k;
  float r2, dx[3], ir, w;
  float a[3] = {0.0, 0.0, 0.0};
  double pix[3] = {parts[i].x[0], parts[i].x[1], parts[i].x[2]};
  int pid = parts[i].id;
  struct cell *cell = root;

  /* Traverse the cells of the tree. */
  while (cell != NULL) {

    /* Are we in a leaf ? */
    if (!cell->split) {

      /* Interact the particle with the particles in the leaf */
      for (j = 0; j < cell->count; ++j) {
        if (cell->parts[j].id == pid) continue;

#if ICHECK >= 0
        if (pid == monitor)
          message("[BH_] Interaction with particle id= %d", cell->parts[j].id);
#endif

        /* Compute the pairwise distance. */
        for (r2 = 0.0, k = 0; k < 3; k++) {
          dx[k] = cell->parts[j].x[k] - pix[k];
          r2 += dx[k] * dx[k];
        }

        /* Apply the gravitational acceleration. */
        ir = 1.0f / sqrtf(r2);
        w = cell->parts[j].mass * const_G * ir * ir * ir;
        for (k = 0; k < 3; k++) a[k] += w * dx[k];

        (*countPairs)++;
      }

      cell = cell->sibling;
    } else {

      /* We are in a node */
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = cell->legacy.com[k] - pix[k];
        r2 += dx[k] * dx[k];
      }

#if ICHECK >= 0
      if (pid == monitor)
        message("This is a node with %d particles h= %f. r= %f theta= %f",
                cell->count, cell->h, sqrt(r2), dist_min);
#endif

      /* Is the cell far enough ? */
      if (dist_min * dist_min * r2 < cell->h * cell->h) {

#if ICHECK >= 0
        if (pid == monitor) printf("Recursing...\n");
#endif
        cell = cell->firstchild;
        continue;
      }

#if ICHECK >= 0
      if (pid == monitor)
        message("[BH_] Can interact with the monopole. x= %f %f %f m= %f h= %f",
                cell->legacy.com[0], cell->legacy.com[1], cell->legacy.com[2],
                cell->legacy.mass, cell->h);
#endif

      /* Apply the gravitational acceleration. */
      ir = 1.0f / sqrtf(r2);
      w = cell->legacy.mass * const_G * ir * ir * ir;
      for (k = 0; k < 3; k++) a[k] += w * dx[k];

      (*countMultipoles)++;

      /* Move to the next node */
      cell = cell->sibling;
    }
  }

  /* Store the locally computed acceleration back to the particle. */
  for (k = 0; k < 3; k++) parts[i].a_legacy[k] += a[k];
}

/**
 * @brief Does a tree walk as in the B-H original work for all particles
 *
 * @param N The number of particles
 * @param parts The array of particles
 * @param root The root cell of the tree
 * @param monitor ID of the particle to monitor and output interactions to
 *        stdout
 */
void legacy_tree_walk(int N, struct part *parts, struct cell *root, int monitor,
                      int *countMultipoles, int *countPairs, int *countCoMs) {

  int i;

  /* Compute multipoles (recursively) */
  legacy_comp_com(root, countCoMs);

  //#pragma omp parallel for
  for (i = 0; i < N; ++i) {
    if (parts[i].id == monitor)
      message("tree walk for particle %d x= %f %f %f", parts[i].id,
              parts[i].x[0], parts[i].x[1], parts[i].x[2]);

    legacy_interact(parts, i, root, monitor, countMultipoles, countPairs);

    if (parts[i].id == monitor)
      message("\n[LEGACY] acceleration for particle %d a= %.3e %.3e %.3e",
              parts[i].id, parts[i].a_legacy[0], parts[i].a_legacy[1],
              parts[i].a_legacy[2]);
  }
}

/* -------------------------------------------------------------------------- */
/* Exact interaction */
/* -------------------------------------------------------------------------- */

/**
 * @brief Solve the particle interactions using the stupid N^2 algorithm
 *
 * @param N The number of particles
 * @param parts The array of particles
 */
void interact_exact(int N, struct part *parts, int monitor) {

  int i, j, k;
  float ir, w, r2, dx[3];

  /* Loop over all particles. */
  for (i = 0; i < N; ++i) {

    /* Some things to store locally. */
    double pix[3] = {parts[i].x[0], parts[i].x[1], parts[i].x[2]};
    float mi = parts[i].mass;

    /* Loop over every other particle. */
    for (j = i + 1; j < N; ++j) {

      /* Compute the pairwise distance. */
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = parts[j].x[k] - pix[k];
        r2 += dx[k] * dx[k];
      }

      /* Apply the gravitational acceleration. */
      ir = 1.0f / sqrtf(r2);
      w = const_G * ir * ir * ir;

      for (k = 0; k < 3; k++) {
        float wdx = w * dx[k];
        parts[j].a_exact[k] -= wdx * mi;
        parts[i].a_exact[k] += wdx * parts[j].mass;
      }
    }
  }

  for (i = 0; i < N; ++i)
    if (parts[i].id == monitor)
      message("[EXACT ] acceleration for particle %d a= %.3e %.3e %.3e\n",
              parts[i].id, parts[i].a_exact[0], parts[i].a_exact[1],
              parts[i].a_exact[2]);
}

/**
 * @brief Set up and run a task-based Barnes-Hutt N-body solver.
 *
 * @param N The number of random particles to use.
 * @param nr_threads Number of threads to use.
 * @param runs Number of force evaluations to use as a benchmark.
 * @param fileName Input file name. If @c NULL or an empty string, random
 *        particle positions will be used.
 */
void test_bh(int N, int nr_threads, int runs, char *fileName) {

  int i, k;
  struct cell *root;
  struct part *parts;
  FILE *file;
  struct qsched s;
  ticks tic, toc_run, tot_setup = 0, tot_run = 0;
  int countMultipoles, countPairs, countCoMs;

  /* Runner function. */
  void runner(int type, void * data) {

    ticks tic = getticks();

    /* Decode the data. */
    struct cell **d = (struct cell **)data;

    /* Decode and execute the task. */
    switch (type) {
      case task_type_self:
        iact_self(d[0]);
        break;
      case task_type_pair:
        iact_pair(d[0], d[1]);
        break;
      case task_type_pair_pc:
        iact_pair_pc(d[0], d[1]);
        break;
      case task_type_pair_direct:
        iact_pair_direct(d[0], d[1]);
        break;
      case task_type_com:
        comp_com(d[0]);
        break;
      default:
        error("Unknown task type.");
    }

    atomic_add(&task_timers[type], getticks() - tic);
  }

  /* Initialize the per-task type timers. */
  for (k = 0; k < task_type_count; k++) task_timers[k] = 0;

  /* Initialize the scheduler. */
  qsched_init(&s, nr_threads, qsched_flag_noreown);

  /* Init and fill the particle array. */
  if ((parts = (struct part *)malloc(sizeof(struct part) * N)) == NULL)
    error("Failed to allocate particle buffer.");

  /* If no input file was specified, generate random particle positions. */
  if (fileName == NULL || fileName[0] == 0) {
    for (k = 0; k < N; k++) {
      parts[k].id = k;
      parts[k].x[0] = ((double)rand()) / RAND_MAX;
      parts[k].x[1] = ((double)rand()) / RAND_MAX;
      parts[k].x[2] = ((double)rand()) / RAND_MAX;
      parts[k].mass = ((double)rand()) / RAND_MAX;
      parts[k].a_legacy[0] = 0.0;
      parts[k].a_legacy[1] = 0.0;
      parts[k].a_legacy[2] = 0.0;
    }

    /* Otherwise, read them from a file. */
  } else {
    file = fopen(fileName, "r");
    if (file) {
      for (k = 0; k < N; k++) {
        if (fscanf(file, "%d", &parts[k].id) != 1)
          error("Failed to read ID of part %i.", k);
        if (fscanf(file, "%lf%lf%lf", &parts[k].x[0], &parts[k].x[1],
                   &parts[k].x[2]) !=
            3)
          error("Failed to read position of part %i.", k);
        if (fscanf(file, "%f", &parts[k].mass) != 1)
          error("Failed to read mass of part %i.", k);
      }
      fclose(file);
    }
  }

  /* Init the cells. */
  root = cell_get();
  root->loc[0] = 0.0;
  root->loc[1] = 0.0;
  root->loc[2] = 0.0;
  root->h = 1.0;
  root->count = N;
  root->parts = parts;
  cell_split(root, &s);

#if ICHECK > 0
  printf("----------------------------------------------------------\n");

  /* Do a N^2 interactions calculation */

  ticks tic_exact = getticks();
  interact_exact(N, parts, ICHECK);
  ticks toc_exact = getticks();

  printf("Exact calculation (1 thread) took %lli (= %e) ticks\n",
         toc_exact - tic_exact, (float)(toc_exact - tic_exact));

  printf("----------------------------------------------------------\n");
#endif

  /* Create the tasks. */
  tic = getticks();
  create_tasks(&s, root, NULL);
  tot_setup += getticks() - tic;

  /* Dump the number of tasks. */
  message("total nr of tasks: %i.", s.count);
  message("total nr of deps: %i.", s.count_deps);
  message("total nr of res: %i.", s.count_res);
  message("total nr of locks: %i.", s.count_locks);
  message("total nr of uses: %i.", s.count_uses);
  int counts[task_type_count];
  for (k = 0; k < task_type_count; k++) counts[k] = 0;
  for (k = 0; k < s.count; k++) counts[s.tasks[k].type] += 1;

  char buffer[200];
  sprintf(buffer, "timings_legacy_%d_%d.dat", cell_maxparts, nr_threads);
  FILE *fileTime = fopen(buffer, "w");

  /* Loop over the number of runs. */
  for (k = 0; k < runs; k++) {

    countMultipoles = 0;
    countPairs = 0;
    countCoMs = 0;

    /* Execute the legacy walk. */
    tic = getticks();
    legacy_tree_walk(N, parts, root, ICHECK, &countMultipoles, &countPairs,
                     &countCoMs);
    toc_run = getticks();

    /* Dump some timings. */
    message("%ith run took %lli (= %e) ticks...", k, toc_run - tic,
            (float)(toc_run - tic));
    tot_run += toc_run - tic;
    fprintf(fileTime, "%lli %e\n", toc_run - tic, (float)(toc_run - tic));
  }

  fclose(fileTime);

#if ICHECK >= 0
  message("[check] accel of part %i is [%.3e,%.3e,%.3e]", ICHECK,
          root->parts[ICHECK].a[0], root->parts[ICHECK].a[1],
          root->parts[ICHECK].a[2]);
#endif
  printf("task counts: [ %8s %8s %8s %8s %8s ]\n", "self", "pair", "m-poles",
         "direct", "CoMs");
  printf("task counts: [ %8i %8i %8i %8i %8i ] (legacy).\n", 0, 0,
         countMultipoles, countPairs, countCoMs);
  printf("task counts: [ ");
  for (k = 0; k < task_type_count; k++) printf("%8i ", counts[k]);
  printf("] (new).\n");

  /* Loop over the number of runs. */
  for (k = 0; k < runs; k++) {

    for (i = 0; i < N; ++i) {
      parts[i].a[0] = 0.0;
      parts[i].a[1] = 0.0;
      parts[i].a[2] = 0.0;
    }

    /* Execute the tasks. */
    tic = getticks();
    qsched_run(&s, nr_threads, runner);
    toc_run = getticks();
    message("%ith run took %lli (= %e) ticks...", k, toc_run - tic,
            (float)(toc_run - tic));
    tot_run += toc_run - tic;
  }

  message("[check] root mass= %f %f", root->legacy.mass, root->new.mass);
  message("[check] root CoMx= %f %f", root->legacy.com[0], root->new.com[0]);
  message("[check] root CoMy= %f %f", root->legacy.com[1], root->new.com[1]);
  message("[check] root CoMz= %f %f", root->legacy.com[2], root->new.com[2]);
#if ICHECK >= 0
  message("[check] accel of part %i is [%.3e,%.3e,%.3e]", ICHECK,
          root->parts[ICHECK].a[0], root->parts[ICHECK].a[1],
          root->parts[ICHECK].a[2]);
#endif

  /* Dump the tasks. */
  /* for ( k = 0 ; k < s.count ; k++ ) */
  /*     printf( " %i %i %lli %lli\n" , s.tasks[k].type , s.tasks[k].qid ,
   * s.tasks[k].tic , s.tasks[k].toc ); */

  /* Dump the costs. */
  message("costs: setup=%lli ticks, run=%lli ticks.", tot_setup,
          tot_run / runs);

  /* Dump the timers. */
  for (k = 0; k < qsched_timer_count; k++)
    message("timer %s is %lli ticks.", qsched_timer_names[k],
            s.timers[k] / runs);

  /* Dump the per-task type timers. */
  printf("task timers: [ ");
  for (k = 0; k < task_type_count; k++) printf("%lli ", task_timers[k] / runs);
  printf("] ticks.\n");

  /* Dump the particles to a file */
  file = fopen("particle_dump.dat", "w");
  fprintf(file, "# x y z a_exact.x   a_exact.y    a_exact.z    a_legacy.x    "
                "a_legacy.y    a_legacy.z    a_new.x     a_new.y    a_new.z\n");
  for (k = 0; k < N; ++k)
    fprintf(file, "%d %e %e %e %e %e %e %e %e %e %e %e %e\n", parts[k].id,
            parts[k].x[0], parts[k].x[1], parts[k].x[2], parts[k].a_exact[0],
            parts[k].a_exact[1], parts[k].a_exact[2], parts[k].a_legacy[0],
            parts[k].a_legacy[1], parts[k].a_legacy[2], parts[k].a[0],
            parts[k].a[1], parts[k].a[2]);
  fclose(file);

  /* Clean up. */
  qsched_free(&s);
}

/**
 * @brief Main function.
 */

int main(int argc, char *argv[]) {

  int c, nr_threads;
  int N = 1000, runs = 1;
  char fileName[100] = {0};

  /* Die on FP-exceptions. */
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

/* Get the number of threads. */
#pragma omp parallel shared(nr_threads)
  {
    if (omp_get_thread_num() == 0) nr_threads = omp_get_num_threads();
  }

  /* Parse the options */
  while ((c = getopt(argc, argv, "n:r:t:f:")) != -1) switch (c) {
      case 'n':
        if (sscanf(optarg, "%d", &N) != 1)
          error("Error parsing number of particles.");
        break;
      case 'r':
        if (sscanf(optarg, "%d", &runs) != 1)
          error("Error parsing number of runs.");
        break;
      case 't':
        if (sscanf(optarg, "%d", &nr_threads) != 1)
          error("Error parsing number of threads.");
        omp_set_num_threads(nr_threads);
        break;
      case 'f':
        if (sscanf(optarg, "%s", &fileName[0]) != 1)
          error("Error parsing file name.");
        break;
      case '?':
        fprintf(stderr,
                "Usage: %s [-t nr_threads] [-n N] [-r runs] [-f file]\n",
                argv[0]);
        fprintf(stderr, "Solves the N-body problem using a Barnes-Hut\n"
                        "tree code with N random particles read from a file in "
                        "[0,1]^3 using\n"
                        "nr_threads threads.\n");
        exit(EXIT_FAILURE);
    }

  /* Tree node information */
  printf("Size of cell: %zu bytes.\n", sizeof(struct cell));

  /* Part information */
  printf("Size of part: %zu bytes.\n", sizeof(struct part));

  /* Dump arguments. */
  if (fileName[0] == 0) {
    message("Computing the N-body problem over %i random particles using %i "
            "threads (%i runs).",
            N, nr_threads, runs);
  } else {
    message("Computing the N-body problem over %i particles read from '%s' "
            "using %i threads (%i runs).",
            N, fileName, nr_threads, runs);
  }

  /* Run the test. */
  test_bh(N, nr_threads, runs, fileName);

  return 0;
}

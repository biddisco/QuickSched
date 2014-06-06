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
#include <omp.h>

/* Local includes. */
#include "quicksched.h"

/* Some local constants. */
#define cell_pool_grow 1000
#define cell_maxparts 16
#define task_limit 5000
#define const_G 6.6738e-8
#define dist_min 0.5  // 0.5

#define ICHECK -1

/** Data structure for the particles. */
struct part {
  double x[3];
  double a[3];
  double a_legacy[3];
  double a_exact[3];
  double mass;
  int id;
};

/** Data structure for the BH tree cell. */
struct cell {
  double loc[3];
  double h;

  int split, count;
  struct part *parts;

  struct cell *firstchild; /* Next node if opening */
  struct cell *sibling;    /* Next node */

  /* We keep both CoMs and masses to make sure the comp_com calculation is
   * correct (use an anonymous union to keep variable names compact).  */
  union {

    /* Information for the legacy walk */
    struct {
      double com[3];
      double mass;
    } legacy;

    /* Information for the QuickShed walk */
    struct {
      double com[3];
      double mass;
    } new;
  };

  int res, com_tid;

} __attribute__((aligned(128)));

/** Task types. */
enum task_type {
  task_type_self = 0,
  task_type_pair,
  task_type_pair_pc,
  task_type_com,
  task_type_count
};

/** Per-type timers. */
ticks task_timers[task_type_count];

/** Global variable for the pool of allocated cells. */
struct cell *cell_pool = NULL;

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
    // c->depth = 0;
    // c->parent = 0;
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
      // cp->parent = c;
      // cp->depth = c->depth + 1;
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->h = c->h / 2;
      cp->h = c->h / 2;
      cp->h = c->h / 2;
      cp->res = qsched_addres(s, qsched_owner_none, c->res);
      if (k & 4) cp->loc[0] += cp->h;
      if (k & 2) cp->loc[1] += cp->h;
      if (k & 1) cp->loc[2] += cp->h;
    }

    /* Init the pivots. */
    for (k = 0; k < 3; k++) pivot[k] = c->loc[k] + c->h / 2;

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

    /* Prepare the pointers. */
    for (k = 0; k < 7; k++) {

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

    /* Last progenitor links to the next sibling */
    progenitors[7]->sibling = c->sibling;
    c->firstchild = progenitors[0];

    /* Recurse. */
    for (k = 0; k < 8; k++) cell_split(progenitors[k], s);

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
      double cp_mass = cp->new.mass;
      com[0] += cp->new.com[0] * cp_mass;
      com[1] += cp->new.com[1] * cp_mass;
      com[2] += cp->new.com[2] * cp_mass;
      mass += cp_mass;
    }

    /* Otherwise, collect the multipole from the particles. */
  } else {

    for (k = 0; k < count; k++) {
      p = &parts[k];
      double p_mass = p->mass;
      com[0] += p->x[0] * p_mass;
      com[1] += p->x[1] * p_mass;
      com[2] += p->x[2] * p_mass;
      mass += p_mass;
    }
  }

  /* Store the COM data, if any was collected. */
  if (mass > 0.0) {
    double imass = 1.0 / mass;
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
void iact_pair_pc(struct cell *ci, struct cell *cj) {
  int j, k, count = ci->count;
  double com[3], mcom, dx[3], r2, ir, w;
  struct part *parts = ci->parts;

  /* Early abort? */
  if (count == 0 || cj->count == 0) return;

  /* message( "ci=[%.3e,%.3e,%.3e], cj=[%.3e,%.3e,%.3e], h=%.3e/%.3e.",
      ci->loc[0], ci->loc[1], ci->loc[2],
      cj->loc[0], cj->loc[1], cj->loc[2],
      ci->h, cj->h ); */

  /* Sanity check. */
  if (cj->new.mass == 0.0) {
    message("%e %e %e %d %p %d %p", cj->new.com[0], cj->new.com[1],
            cj->new.com[2], cj->count, cj, cj->split, cj->sibling);

    for (j = 0; j < cj->count; ++j)
      message("part %d mass=%e id=%d", j, cj->parts[j].mass, cj->parts[j].id);

    error("com does not seem to have been set.");
  }

  /* Corectness check */
  if (cj->new.mass != cj->legacy.mass)
    error("Calculation of the CoM is wrong! m_new=%e m_legacy=%e", cj->new.mass,
          cj->legacy.mass);

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
    ir = 1.0 / sqrt(r2);
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
 * @brief Compute the interactions between all particles in a cell.
 *
 * @param ci The #cell.
 * @param cj The other #cell.
 */
void iact_pair(struct cell *ci, struct cell *cj) {

  int i, j, k;
  int count_i = ci->count, count_j = cj->count;
  double dx[3], xi[3], ai[3], mi, mj, r2, r2_i, r2_j, w, ir;
  double cih = ci->h, cjh = cj->h;
  struct part *parts_i = ci->parts, *parts_j = cj->parts;
  struct cell *cp;

  /* Early abort? */
  if (count_i == 0 || count_j == 0) return;

  /* Sanity check */
  if (ci == cj)
    error("The impossible has happened: pair interaction between a cell and "
          "itself.");  // debug

  /* Distance between the CoMs */
  r2 = 0.0;
  r2_i = 0.0;
  r2_j = 0.0;
  for (k = 0; k < 3; k++) {
    //   dx[k] = fabs( ci->new.com[k] - cj->new.com[k] );
    dx[k] = fabs(ci->loc[k] - cj->loc[k]);

    r2 += dx[k] * dx[k];
    r2_i += (dx[k] - 0.5 * cih) * (dx[k] - 0.5 * cih);
    r2_j += (dx[k] - 0.5 * cjh) * (dx[k] - 0.5 * cjh);
  }

  /* If ci and cj are sufficiently well separated, split this interaction
     into a pair of particle-cell interactions. */
  if ((dist_min * dist_min * r2_j > ci->h * ci->h) &&
      (dist_min * dist_min * r2_i > cj->h * cj->h)) {
    iact_pair_pc(ci, cj);
    iact_pair_pc(cj, ci);

    /* Otherwise, if neither cell is split, compute the particle-particle
       interactions directly. */
  } else if (!ci->split && !cj->split) {

    /* Loop over all particles in ci... */
    for (i = 0; i < count_i; i++) {

      /* Init the ith particle's data. */
      for (k = 0; k < 3; k++) {
        xi[k] = parts_i[i].x[k];
        ai[k] = 0.0;
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
        ir = 1.0 / sqrt(r2);
        w = const_G * ir * ir * ir;
        mj = parts_j[j].mass;
        for (k = 0; k < 3; k++) {
          double wdx = w * dx[k];
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

    /* Otherwise, compute the interaction recursively over the progeny. */
  } else {

    /* We can split one of the two cells. Let's try the biggest one first.*/
    if (ci->h > cj->h) {
      if (ci->split) {
        for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling)
          iact_pair(cp, cj);
      } else if (cj->split) {
        for (cp = cj->firstchild; cp != cj->sibling; cp = cp->sibling)
          iact_pair(ci, cp);
      }

    } else {
      if (cj->split) {
        for (cp = cj->firstchild; cp != cj->sibling; cp = cp->sibling)
          iact_pair(ci, cp);
      } else if (ci->split) {
        for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling)
          iact_pair(cp, cj);
      }
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
  double xi[3], ai[3], mi, mj, dx[3], r2, ir, w;
  struct part *parts = c->parts;
  struct cell *cp, *cps;

  /* Early abort? */
  if (count == 0) return;

  /* message( "cell=[%.3e,%.3e,%.3e], h=%.3e.",
      c->loc[0], c->loc[1], c->loc[2], c->h ); */

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
        ir = 1.0 / sqrt(r2);
        w = const_G * ir * ir * ir;
        mj = parts[j].mass;
        for (k = 0; k < 3; k++) {
          double wdx = w * dx[k];
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
  struct cell *data[2], *cp, *cps;

  /* If either cell is empty, stop. */
  if (ci->count == 0 || (cj != NULL && cj->count == 0)) return;

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

    double dx, r2 = 0.0, r2_i = 0.0, r2_j = 0.0;

    /* Distance between the cells */
    for (k = 0; k < 3; k++) {
      dx = fabs(ci->loc[k] - cj->loc[k]);

      r2 += dx * dx;
      r2_i += (dx - 0.5 * ci->h) * (dx - 0.5 * ci->h);
      r2_j += (dx - 0.5 * cj->h) * (dx - 0.5 * cj->h);
    }

    /* Check whether we can use the multipoles. */
    if ((dist_min * dist_min * r2_j > ci->h * ci->h) &&
        (dist_min * dist_min * r2_i > cj->h * cj->h)) {
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

      /* Otherwise, if neither cells are split, generate a part-part task. */
    } else if (!ci->split && !cj->split) {

      /* Set the data. */
      data[0] = ci;
      data[1] = cj;

      /* Create the task. */
      tid = qsched_addtask(s, task_type_pair, task_flag_none, data,
                           sizeof(struct cell *) * 2, ci->count * cj->count);

      /* Add the resources. */
      qsched_addlock(s, tid, ci->res);
      qsched_addlock(s, tid, cj->res);
      qsched_addunlock(s, ci->com_tid, tid);
      qsched_addunlock(s, cj->com_tid, tid);

      /* Otherwise, compute the interaction recursively over the progeny. */
    } else if (ci->count > task_limit && cj->count > task_limit) {

      /* We can split one of the two cells. Let's try the biggest one */
      if (ci->h > cj->h) {

        if (ci->split) {
          for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling)
            create_tasks(s, cp, cj);
        } else if (cj->split) {
          for (cp = cj->firstchild; cp != cj->sibling; cp = cp->sibling)
            create_tasks(s, ci, cp);
        }

      } else {

        if (cj->split) {
          for (cp = cj->firstchild; cp != cj->sibling; cp = cp->sibling)
            create_tasks(s, ci, cp);
        } else if (ci->split) {
          for (cp = ci->firstchild; cp != ci->sibling; cp = cp->sibling)
            create_tasks(s, cp, cj);
        }
      }

      /* Create a task if too few particles */
    } else {
      /* Set the data. */
      data[0] = ci;
      data[1] = cj;

      /* Create the task. */
      tid = qsched_addtask(s, task_type_pair, task_flag_none, data,
                           sizeof(struct cell *) * 2, ci->count * cj->count);

      /* Add the resources. */
      qsched_addlock(s, tid, ci->res);
      qsched_addlock(s, tid, cj->res);

      /* Depend on the COMs in case this task recurses. */
      if (ci->split || cj->split) {
        qsched_addunlock(s, ci->com_tid, tid);
        qsched_addunlock(s, cj->com_tid, tid);
      }
    }

  } /* otherwise, it's a pair. */
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
      double cp_mass = cp->legacy.mass;
      com[0] += cp->legacy.com[0] * cp_mass;
      com[1] += cp->legacy.com[1] * cp_mass;
      com[2] += cp->legacy.com[2] * cp_mass;
      mass += cp_mass;
    }

    /* Otherwise, collect the multipole from local data. */
  } else {

    for (k = 0; k < count; k++) {
      p = &parts[k];
      double p_mass = p->mass;
      com[0] += p->x[0] * p_mass;
      com[1] += p->x[1] * p_mass;
      com[2] += p->x[2] * p_mass;
      mass += p_mass;
    }
  }

  /* Finish multipole calculation */
  if (mass > 0.0) {
    double imass = 1.0 / mass;
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
  double r2, dx[3], ir, w;
  double a[3] = {0.0, 0.0, 0.0};
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
        ir = 1.0 / sqrt(r2);
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
      ir = 1.0 / sqrt(r2);
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
      message("\n[LEGACY] acceleration for particle %d a= %f %f %f",
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
  double ir, w, r2, dx[3];

  /* Loop over all particles. */
  for (i = 0; i < N; ++i) {

    /* Some things to store locally. */
    double pix[3] = {parts[i].x[0], parts[i].x[1], parts[i].x[2]};
    double mi = parts[i].mass;

    /* Loop over every other particle. */
    for (j = i + 1; j < N; ++j) {

      /* Compute the pairwise distance. */
      for (r2 = 0.0, k = 0; k < 3; k++) {
        dx[k] = parts[j].x[k] - pix[k];
        r2 += dx[k] * dx[k];
      }

      /* Apply the gravitational acceleration. */
      ir = 1.0 / sqrt(r2);
      w = const_G * ir * ir * ir;

      for (k = 0; k < 3; k++) {
        double wdx = w * dx[k];
        parts[j].a_exact[k] -= wdx * mi;
        parts[i].a_exact[k] += wdx * parts[j].mass;
      }
    }

  }

  for (i = 0; i < N; ++i)
    if (parts[i].id == monitor)
      message("[EXACT ] acceleration for particle %d a= %f %f %f\n",
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
  ticks tic, toc_run, tot_setup = 0, tot_run = 0, tic_exact, toc_exact;
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
        if (fscanf(file, "%lf", &parts[k].mass) != 1)
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

  printf("----------------------------------------------------------\n");

  /* Do a N^2 interactions calculation */

  tic_exact = getticks();
  //interact_exact( N , parts , ICHECK );
  toc_exact = getticks();

  printf("Exact calculation (1 thread) took %lli (= %e) ticks\n",
         toc_exact - tic_exact, (float)(toc_exact - tic_exact));

  printf("----------------------------------------------------------\n");

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

  printf("task counts: [ %8s %8s %8s %8s ]\n", "self", "direct", "m-poles",
         "CoMs");
  printf("task counts: [ %8i %8i %8i %8i ] (legacy).\n", 0, countPairs,
         countMultipoles, countCoMs);
  printf("task counts: [ ");
  for (k = 0; k < task_type_count; k++) printf("%8i ", counts[k]);
  printf("] (new).\n");

  /* Loop over the number of runs. */
  for (k = 0; k < runs; k++) {

    for (i = 0; i < N; ++i) {
      parts[i].a[0] = 0.;
      parts[i].a[1] = 0.;
      parts[i].a[2] = 0.;
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
  /*file = fopen("particle_dump.dat", "w");
  fprintf(file, "# x y z a_exact.x   a_exact.y    a_exact.z    a_legacy.x    "
                "a_legacy.y    a_legacy.z    a_new.x     a_new.y    a_new.z\n");
  for (k = 0; k < N; ++k)
    fprintf(file, "%d %e %e %e %e %e %e %e %e %e %e %e %e\n", parts[k].id,
            parts[k].x[0], parts[k].x[1], parts[k].x[2], parts[k].a_exact[0],
            parts[k].a_exact[1], parts[k].a_exact[2], parts[k].a_legacy[0],
            parts[k].a_legacy[1], parts[k].a_legacy[2], parts[k].a[0],
            parts[k].a[1], parts[k].a[2]);
  fclose(file);*/

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
        fprintf(stderr, "Solves the N-body problem using a Barnes-Hutt\n"
                        "tree code with N random particles read from a file in "
                        "[0,1]^3 using\n"
                        "nr_threads threads.\n");
        exit(EXIT_FAILURE);
    }

  /* Tree node information */
  printf("Size of cell: %zu bytes.\n", sizeof(struct cell));

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

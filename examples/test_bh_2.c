/*******************************************************************************
 * This file is part of QuickSched.
 * Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
 ******************************************************************************/


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
#define cell_pool_grow 100
#define cell_maxparts 20
#define task_limit 5000
#define const_G 6.6738e-11
#define dist_min 0.5 // 0.5


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
    double h[3];
    int split, count;
    struct part *parts;
    struct cell *progeny[8];
    struct cell *parent;
    double com[3];
    double mass;
    double com_legacy[3];
    double mass_legacy;
    int res, com_tid, depth;
    };
    
    
/** Task types. */
enum task_type {
    task_type_self = 0,
    task_type_pair,
    task_type_pair_pc,
    task_type_com,
    task_type_count
    };
    
/** Per-type timers. */
ticks task_timers[ task_type_count ];
    
    
/** Global variable for the pool of allocated cells. */
struct cell *cell_pool = NULL;


/**
 * @brief Get a #cell from the pool.
 */
 
struct cell *cell_get ( ) {

    struct cell *res;
    int k;

    /* Allocate a new batch? */
    if ( cell_pool == NULL ) {
    
        /* Allocate the cell array. */
        if ( ( cell_pool = (struct cell *)malloc( sizeof(struct cell) * cell_pool_grow ) ) == NULL )
            error( "Failed to allocate fresh batch of cells." );
            
        /* Clear the cells. */
        bzero( cell_pool , sizeof(struct cell) * cell_pool_grow );
        
        /* Link them up via their progeny pointers. */
        for ( k = 1 ; k < cell_pool_grow ; k++ )
            cell_pool[k-1].progeny[0] = &cell_pool[k];
    
        }
        
    /* Pick a cell off the pool. */
    res = cell_pool;
    cell_pool = cell_pool->progeny[0];
    
    /* Clean up a few things. */
    res->res = qsched_res_none;
    
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
 
void cell_split ( struct cell *c , struct qsched *s ) {

    int i, j, k, count = c->count;
    struct part temp, *parts = c->parts;
    struct cell *cp;
    int left[8], right[8];
    double pivot[3];
    static struct cell *root = NULL;
    
    /* Set the root cell. */
    if ( root == NULL ) {
        root = c;
        c->depth = 0;
        }
    
    /* Add a resource for this cell if it doesn't have one yet. */
    if ( c->res == qsched_res_none ) 
        c->res = qsched_addres( s , qsched_owner_none , qsched_res_none );
    
    /* Attach a center-of-mass task to the cell. */
    c->com_tid = qsched_addtask( s , task_type_com , task_flag_none , &c , sizeof(struct cell *) , 1 );
    
    /* Does this cell need to be split? */
    if ( count > cell_maxparts ) {
    
        /* Mark this cell as split. */
        c->split = 1;
    
        /* Create the progeny. */
        for ( k = 0 ; k < 8 ; k++ ) {
            c->progeny[k] = cp = cell_get();
            cp->parent = c;
            cp->depth = c->depth + 1;
            cp->loc[0] = c->loc[0];
            cp->loc[1] = c->loc[1];
            cp->loc[2] = c->loc[2];
            cp->h[0] = c->h[0]/2;
            cp->h[1] = c->h[1]/2;
            cp->h[2] = c->h[2]/2;
            cp->res = qsched_addres( s , qsched_owner_none , c->res );
            if ( k & 4 )
                cp->loc[0] += cp->h[0];
            if ( k & 2 )
                cp->loc[1] += cp->h[1];
            if ( k & 1 )
                cp->loc[2] += cp->h[2];
            }
    
        /* Init the pivots. */
        for ( k = 0 ; k < 3 ; k++ )
            pivot[k] = c->loc[k] + c->h[k]/2;

        /* Split along the x-axis. */
        i = 0; j = count - 1;
        while ( i <= j ) {
            while ( i <= count-1 && parts[i].x[0] < pivot[0] )
                i += 1;
            while ( j >= 0 && parts[j].x[0] >= pivot[0] )
                j -= 1;
            if ( i < j ) {
                temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                }
            }
        left[1] = i; right[1] = count - 1;
        left[0] = 0; right[0] = j;

        /* Split along the y axis, twice. */
        for ( k = 1 ; k >= 0 ; k-- ) {
            i = left[k]; j = right[k];
            while ( i <= j ) {
                while ( i <= right[k] && parts[i].x[1] < pivot[1] )
                    i += 1;
                while ( j >= left[k] && parts[j].x[1] >= pivot[1] )
                    j -= 1;
                if ( i < j ) {
                    temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                    }
                }
            left[2*k+1] = i; right[2*k+1] = right[k];
            left[2*k] = left[k]; right[2*k] = j;
            }

        /* Split along the z axis, four times. */
        for ( k = 3 ; k >= 0 ; k-- ) {
            i = left[k]; j = right[k];
            while ( i <= j ) {
                while ( i <= right[k] && parts[i].x[2] < pivot[2] )
                    i += 1;
                while ( j >= left[k] && parts[j].x[2] >= pivot[2] )
                    j -= 1;
                if ( i < j ) {
                    temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                    }
                }
            left[2*k+1] = i; right[2*k+1] = right[k];
            left[2*k] = left[k]; right[2*k] = j;
            }

        /* Store the counts and offsets. */
        for ( k = 0 ; k < 8 ; k++ ) {
            c->progeny[k]->count = right[k] - left[k] + 1;
            c->progeny[k]->parts = &c->parts[ left[k] ];
            }
            
        /* Recurse. */
        for ( k = 0 ; k < 8 ; k++ )
            cell_split( c->progeny[k] , s );
            
        /* Link the COM tasks. */
        for ( k = 0 ; k < 8 ; k++ )
            qsched_addunlock( s , c->progeny[k]->com_tid , c->com_tid );
            
        } /* does the cell need to be split? */
        
    /* Set this cell's resources ownership. */
    qsched_res_own( s , c->res , s->nr_queues * ( c->parts - root->parts ) / root->count );
        
    }
    
    
/**
 * @brief Compute the center of mass of a given cell.
 *
 * @param c The #cell.
 */
 
void comp_com ( struct cell *c ) {

    int k, count = c->count;
    struct cell *cp;
    struct part *p, *parts = c->parts;

    c->com[0] = c->com[1] = c->com[2] = c->mass = 0.;

    /* Is the cell split? */
    if ( c->split ) {
    
        /* Collect the centre of gravity and mass from the progeny. */
        for ( k = 0 ; k < 8 ; k++ ) {
            cp = c->progeny[k];
            c->com[0] += cp->com[0]*cp->mass;
            c->com[1] += cp->com[1]*cp->mass;
            c->com[2] += cp->com[2]*cp->mass;
            c->mass += cp->mass;
            }

	if ( c-> mass > 0. )
	  {
	    c->com[0] /= c->mass; 
	    c->com[1] /= c->mass; 
	    c->com[2] /= c->mass;
	  }
	else
	  {
	    c->com[0] = 0.;
	    c->com[1] = 0.;
	    c->com[2] = 0.;
	  }
            
        }
        
    /* Otherwise, collect the multipole from local data. */
    else {
    
        for ( k = 0 ; k < count ; k++ ) {
            p = &parts[k];
            c->com[0] += p->x[0]*p->mass;
            c->com[1] += p->x[1]*p->mass;
            c->com[2] += p->x[2]*p->mass;
            c->mass += p->mass;
            }

	if ( c-> mass > 0. )
	  {
	    c->com[0] /= c->mass; 
	    c->com[1] /= c->mass; 
	    c->com[2] /= c->mass;
	  }
	else
	  {
	    c->com[0] = 0.;
	    c->com[1] = 0.;
	    c->com[2] = 0.;
	  }            
        }
        
    }
    
    
/**
 * @brief Compute the interactions between all particles in a cell
 *        and the center of mass of another cell.
 *
 * @param ci The #cell containing the particles.
 * @param cj The #cell containing the center of mass.
 */
 
void iact_pair_pc ( struct cell *ci , struct cell *cj ) {

    int j, k, count = ci->count;
    double com[3], mcom, dx[3], r2, ir, w;
    struct part *parts = ci->parts;
    
    /* Early abort? */
    if ( count == 0 || cj->count == 0 )
        return;
        
    /* message( "ci=[%.3e,%.3e,%.3e], cj=[%.3e,%.3e,%.3e], h=%.3e/%.3e.",
        ci->loc[0], ci->loc[1], ci->loc[2], 
        cj->loc[0], cj->loc[1], cj->loc[2],
        ci->h[0], cj->h[0] ); */
    
    /* Sanity check. */
    if ( cj->mass == 0.0 ){
      printf("%e %e %e %d %p\n", cj->com[0], cj->com[1], cj->com[2], cj->count, cj);

      for ( j = 0 ; j < cj->count ; ++j )
	printf("part %d mass= %e\n", j, cj->parts[j].mass );

        error( "com does not seem to have been set." );
    }
    /* Init the com's data. */
    for ( k = 0 ; k < 3 ; k++ )
        com[k] = cj->com[k];
    mcom = cj->mass;

    /* Loop over every particle in ci. */
    for ( j = 0 ; j < count ; j++ ) {

        /* Compute the pairwise distance. */
        for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
            dx[k] = com[k] - parts[j].x[k];
            r2 += dx[k]*dx[k];
            }

        /* Apply the gravitational acceleration. */
        ir = 1.0 / sqrt( r2 );
        w = mcom * const_G * ir * ir * ir;
        for ( k = 0 ; k < 3 ; k++ )
            parts[j].a[k] += w * dx[k];

#if ICHECK >= 0
	if ( parts[j].id == ICHECK )
	  printf("[NEW] Can interact with the monopole. x= %f %f %f m= %f h= %f\n", com[0], com[1], com[2], mcom, cj->h[0]);
#endif

        } /* loop over every particle. */
                
    }
    

/**
 * @brief Compute the interactions between all particles in a cell.
 *
 * @param ci The #cell.
 * @param cj The other #cell.
 */
 
void iact_pair ( struct cell *ci , struct cell *cj ) {
  
  int i, j, k;
    int count_i = ci->count, count_j = cj->count;
    double dx[3], xi[3], ai[3], mi, mj, r2, w, ir;
    struct part *parts_i = ci->parts, *parts_j = cj->parts;
  
    /* Early abort? */
    if ( count_i == 0 || count_j == 0 )
        return;

      /* Distance between the CoMs */
      for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
        dx[k] = fabs( ci->com[k] - cj->com[k] );	
        r2 += dx[k]*dx[k];
      }

      double s_max_i = ci->h[0]; 
      double s_max_j = cj->h[0]; 
      
      if ( ( dist_min * dist_min * r2 > s_max_i * s_max_i ) && ( dist_min * dist_min * r2 > s_max_j * s_max_j ) )
	{
	  iact_pair_pc( ci, cj );
	  iact_pair_pc( cj, ci );
	}
      else if ( ci->split == 0 && cj->split == 0 )
	{
	  /* Do direct summation */
	
        /* Loop over all particles... */
        for ( i = 0 ; i < count_i ; i++ ) {

            /* Init the ith particle's data. */
            for ( k = 0 ; k < 3 ; k++ ) {
                xi[k] = parts_i[i].x[k];
                ai[k] = 0.0;
                }
            mi = parts_i[i].mass;

            /* Loop over every following particle. */
            for ( j = 0 ; j < count_j ; j++ ) {

                /* Compute the pairwise distance. */
                for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
                    dx[k] = xi[k] - parts_j[j].x[k];
                    r2 += dx[k]*dx[k];
                    }

                /* Apply the gravitational acceleration. */
                ir = 1.0 / sqrt( r2 );
                w = const_G * ir * ir * ir;
                mj = parts_j[j].mass;
                for ( k = 0 ; k < 3 ; k++ ) {
                    parts_j[j].a[k] += w * dx[k] * mi;
                    ai[k] -= w * dx[k] * mj;
                    }

#if ICHECK >= 0
		if ( parts_i[i].id == ICHECK )
		  printf("[NEW] Interaction with particle id= %d (pair i)\n", parts_j[j].id);

		if ( parts_j[j].id == ICHECK )
		  printf("[NEW] Interaction with particle id= %d (pair j) h_i= %f h_j= %f ci= %p cj= %p count_i= %d count_j= %d d_i= %d d_j= %d\n", parts_i[i].id, ci->h[0], cj->h[0], ci, cj, count_i, count_j,
			 ci->res, cj->res ) ;
#endif
                } /* loop over every other particle. */

            /* Store the accumulated acceleration on the ith part. */
            for ( k = 0 ; k < 3 ; k++ )
                parts_i[i].a[k] += ai[k];

            } /* loop over all particles. */
	  

	}
      else {

      if ( ci == cj )
	printf("The impossible has happened\n");  //debug


      /* We can split one of the two cells. Let's try the biggest one */
      if ( ci->h[0] > cj->h[0] ) {

	if (  ci->split ) 
	  for ( k = 0 ; k < 8 ; k++ )
	    iact_pair ( ci->progeny[k], cj );
	
	else if ( cj->split ) 
	  for ( k = 0 ; k < 8 ; k++ )
	    iact_pair ( ci, cj->progeny[k] );

	else
	  printf("oO\n");

      }
      else {

	if ( cj->split ) 
	  for ( k = 0 ; k < 8 ; k++ )
	    iact_pair ( ci, cj->progeny[k] );
	
	else if (  ci->split ) 
	  for ( k = 0 ; k < 8 ; k++ )
	    iact_pair ( ci->progeny[k], cj );

	else
	  printf("Oo\n");
	
	    }
    }
 }
    

/**
 * @brief Compute the interactions between all particles in a cell.
 *
 * @param c The #cell.
 */
 
void iact_self ( struct cell *c ) {

    int i, j, k, count = c->count;
    double xi[3], ai[3], mi, mj, dx[3], r2, ir, w;
    struct part *parts = c->parts;
    
    /* Early abort? */
    if ( count == 0 )
        return;

    //    printf("entering cell at x= %f %f %f h= %f\n",  c->loc[0], c->loc[1], c->loc[2], c->h[0] );
    
    /* message( "cell=[%.3e,%.3e,%.3e], h=%.3e.",
        c->loc[0], c->loc[1], c->loc[2], c->h[0] ); */
    
    /* Recurse? */
    if ( c->split )
      {
        for ( j = 0 ; j < 8 ; j++ ) {
	  iact_self( c->progeny[j] );
            for ( k = j+1 ; k < 8 ; k++ )
	      iact_pair( c->progeny[j] , c->progeny[k] );
            }
      }
    /* Otherwise, compute interactions directly. */
    else {

        /* Loop over all particles... */
        for ( i = 0 ; i < count ; i++ ) {
        
            /* Init the ith particle's data. */
            for ( k = 0 ; k < 3 ; k++ ) {
                xi[k] = parts[i].x[k];
                ai[k] = 0.0;
                }
            mi = parts[i].mass;
                
            /* Loop over every following particle. */
            for ( j = i+1 ; j < count ; j++ ) {
            
                /* Compute the pairwise distance. */
                for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
                    dx[k] = xi[k] - parts[j].x[k];
                    r2 += dx[k]*dx[k];
                    }
                    
                /* Apply the gravitational acceleration. */
                ir = 1.0 / sqrt( r2 );
                w = const_G * ir * ir * ir;
                mj = parts[j].mass;
                for ( k = 0 ; k < 3 ; k++ ) {
                    parts[j].a[k] += w * dx[k] * mi;
                    ai[k] -= w * dx[k] * mj;
                    }

#if ICHECK >= 0
		if ( parts[i].id == ICHECK )
		  printf("[NEW] Interaction with particle id= %d (self i)\n", parts[j].id);

		if ( parts[j].id == ICHECK )
		  printf("[NEW] Interaction with particle id= %d (self j)\n", parts[i].id);
#endif
            
                } /* loop over every other particle. */
                
            /* Store the accumulated acceleration on the ith part. */
            for ( k = 0 ; k < 3 ; k++ )
                parts[i].a[k] += ai[k];
        
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
 
void create_tasks ( struct qsched *s , struct cell *ci , struct cell *cj ) {

    int j, k;
    qsched_task_t tid;
    struct cell *data[2];
    double dx, r2;
    
    /* If either cell is empty, stop. */
    if ( ci->count == 0 || ( cj != NULL && cj->count == 0 ) )
        return;

    /* Single cell? */
    if ( cj == NULL ) {
    
        /* Is this cell split? */
        if ( ci->split && ci->count > task_limit ) {

            /* Recurse over the progeny. */
            for ( j = 0 ; j < 8 ; j++ ) {
            
                /* Make self-task. */
                create_tasks( s , ci->progeny[j] , NULL );
                
                /* Make all pair tasks. */
                for ( k = j+1 ; k < 8 ; k++ )
                    create_tasks( s , ci->progeny[j] , ci->progeny[k] );
                }
        
            }
            
        /* Otherwise, add a self-interaction task. */
        else {
            
            /* Set the data. */
            data[0] = ci; data[1] = NULL;
            
            /* Create the task. */
            tid = qsched_addtask( s , task_type_self , task_flag_none , data , sizeof(struct cell *) * 2 , ci->count*ci->count/2 );
            
            /* Add the resource. */
            qsched_addlock( s , tid , ci->res );
            
            /* If this call might recurse, add a dependency on the cell's COM. */
            if ( ci->split )
                qsched_addunlock( s , ci->com_tid , tid );
        
            }
    
        }
        
    /* Otherwise, it's a pair. */
    else {


      //      printf("Pair !\n");

      /* Distance between the CoMs */
      for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
        dx = fabs( ci->loc[k] - cj->loc[k] );	
        r2 += dx*dx;
      }
      
      const double s_max_i = ci->h[0]; 
      const double s_max_j = cj->h[0]; 
      
      /* Check whether we can use the multipoles. */
      if ( ( dist_min * dist_min * r2 > s_max_i * s_max_i ) && ( dist_min * dist_min * r2 > s_max_j * s_max_j ) )
	{  
	  data[0] = ci; data[1] = cj;
	  tid = qsched_addtask( s , task_type_pair_pc , task_flag_none , data , sizeof(struct cell *) * 2 , ci->count );
	  qsched_addlock( s , tid , ci->res );
	  qsched_addunlock( s , cj->com_tid , tid );
	  
	  data[0] = cj; data[1] = ci;
	  tid = qsched_addtask( s , task_type_pair_pc , task_flag_none , data , sizeof(struct cell *) * 2 , cj->count );
	  qsched_addlock( s , tid , cj->res );
	  qsched_addunlock( s , ci->com_tid , tid );
	}
      
            
      /* Otherwise, generate a part-part task. */
      else if ( ci->split == 0 && cj->split == 0 )
	{
        
	  /* Set the data. */
	  data[0] = ci; data[1] = cj;

	  /* Create the task. */
	  tid = qsched_addtask( s , task_type_pair , task_flag_none , data , sizeof(struct cell *) * 2 , ci->count*cj->count );

	  /* Add the resources. */
	  qsched_addlock( s , tid , ci->res );
	  qsched_addlock( s , tid , cj->res );

	  /* qsched_addunlock( s , ci->com_tid , tid ); */
	  /* qsched_addunlock( s , cj->com_tid , tid ); */


	}


      else if ( ci->count > task_limit && cj->count > task_limit )
	{
	  
	  /* if ( ci == cj ) */
	  /*   printf("The impossible has happened\n");  //debug */
	  
	  /* We can split one of the two cells. Let's try the biggest one */
	  if ( ci->h[0] > cj->h[0] ) {
	    
	    if (  ci->split ) 
	      for ( k = 0 ; k < 8 ; k++ )
		create_tasks ( s , ci->progeny[k], cj );
	    
	    else if ( cj->split ) 
	      for ( k = 0 ; k < 8 ; k++ )
		create_tasks ( s , ci, cj->progeny[k] );
	    
	    /* else */
	    /*   printf("oO\n"); */
	    
	  }
	  else {
	    
	    if ( cj->split ) 
	      for ( k = 0 ; k < 8 ; k++ )
		create_tasks ( s , ci, cj->progeny[k] );
	    
	    else if (  ci->split ) 
	      for ( k = 0 ; k < 8 ; k++ )
		create_tasks ( s , ci->progeny[k], cj );
	    
	    /* else */
	    /*   printf("Oo\n"); */
	    
	  }
	}

      /* /\* Create a task if too few particles *\/ */
      else 
      	{
      	  /* Set the data. */
      	  data[0] = ci; data[1] = cj;

      	  /* Create the task. */
      	  tid = qsched_addtask( s , task_type_pair , task_flag_none , data , sizeof(struct cell *) * 2 , ci->count*cj->count );

      	  /* Add the resources. */
      	  qsched_addlock( s , tid , ci->res );
      	  qsched_addlock( s , tid , cj->res );

      	  /* Depend on the COMs in case this task recurses. */
      	  if ( ci->split || cj->split ) {
      	    qsched_addunlock( s , ci->com_tid , tid );
      	    qsched_addunlock( s , cj->com_tid , tid );
      	  }
      	}

	

      } /* otherwise, it's a pair. */
      
    }










/**
 * @brief Compute the center of mass of a given cell recursively.
 *
 * @param c The #cell.
 */
 
void legacy_comp_com ( struct cell *c , int* countCoMs ) {

    int k, count = c->count;
    struct cell *cp;
    struct part *p, *parts = c->parts;

    ++(*countCoMs);

    /* Is the cell split? */
    if ( c->split ) {
    
      
      for ( k = 0 ; k < 8 ; k++ ) {
	legacy_comp_com( c->progeny[k] , countCoMs );
	}

        /* Collect the centre of gravity and mass from the progeny. */
        for ( k = 0 ; k < 8 ; k++ ) {
            cp = c->progeny[k];
            c->com_legacy[0] += cp->com_legacy[0]*cp->mass_legacy;
            c->com_legacy[1] += cp->com_legacy[1]*cp->mass_legacy;
            c->com_legacy[2] += cp->com_legacy[2]*cp->mass_legacy;
            c->mass_legacy += cp->mass_legacy;
            }

	if ( c-> mass_legacy != 0. )
	  {
	    c->com_legacy[0] /= c->mass_legacy;
	    c->com_legacy[1] /= c->mass_legacy;
	    c->com_legacy[2] /= c->mass_legacy;
	  }
	else
	  {
	    c->com_legacy[0] = 0.;
	    c->com_legacy[1] = 0.;
	    c->com_legacy[2] = 0.;
	  }
            
        }
        
    /* Otherwise, collect the multipole from local data. */
    else {
    
        for ( k = 0 ; k < count ; k++ ) {
            p = &parts[k];
            c->com_legacy[0] += p->x[0]*p->mass;
            c->com_legacy[1] += p->x[1]*p->mass;
            c->com_legacy[2] += p->x[2]*p->mass;
            c->mass_legacy += p->mass;
            }

	if ( c-> mass_legacy != 0. )
	  {
	    c->com_legacy[0] /= c->mass_legacy;
	    c->com_legacy[1] /= c->mass_legacy;
	    c->com_legacy[2] /= c->mass_legacy;
	  }
	else
	  {
	    c->com_legacy[0] = 0.;
	    c->com_legacy[1] = 0.;
	    c->com_legacy[2] = 0.;
	  }
        }
        
    }




/**
 * @brief Interacts a particle with a cell recursively using the original B-H tree walk procedure
 *
 * @param parts The array of particles
 * @param i The particle of interest
 * @param cell The cell the particle interacts with
 */
void legacy_interact( struct part* parts, int i , struct cell* cell , int monitor,  int* countMultipoles, int* countPairs) {

  int j, k;
  double r2, ir, w, dx[3];

  /* Early abort */
  if( cell->count == 0 )
    return;

  if( parts[i].id == monitor )
    printf( "Particle %d entering cell at x= %f %f %f h= %f\n", parts[i].id , cell->loc[0], cell->loc[1], cell->loc[2], cell->h[0] );

  /* Are we in a leaf ? */
  if ( ! cell->split )
    {
      if( parts[i].id == monitor )
	printf( "This is a leaf with %d particles.\n", cell->count);

      /* Interact the particle with the particles in the leaf */
      
      for( j = 0 ; j < cell->count ; ++ j )
	{
	  if( cell->parts[j].id == parts[i].id )
	    continue;

	  if( parts[i].id == monitor )
	    printf( "[BH_] Interaction with particle id= %d\n", cell->parts[j].id );

	  /* Compute the pairwise distance. */
	  for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
	    dx[k] = cell->parts[j].x[k] - parts[i].x[k];
	    r2 += dx[k]*dx[k];
	  }
	  
	  /* Apply the gravitational acceleration. */
	  ir = 1.0 / sqrt( r2 );
	  w = cell->parts[j].mass * const_G * ir * ir * ir;
	  for ( k = 0 ; k < 3 ; k++ )
	    parts[i].a_legacy[k] += w * dx[k];

	  (*countPairs)++;
	}
    }
  else
    {
      /* We are in a node */
      for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
        dx[k] = cell->loc[k] - parts[i].x[k];
        r2 += dx[k]*dx[k];
      }
      
      if( parts[i].id == monitor )
	printf( "This is a node with %d particles h= %f. r= %f theta= %f\n", cell->count, cell->h[0], sqrt( r2 ), dist_min );


      /* Is the cell far enough ? */
      if( dist_min*dist_min*r2 > cell->h[0]*cell->h[0]  )
	{
	  /* Interact with the monopole */
	  if( parts[i].id == monitor )
	    printf( "[BH_] Can interact with the monopole. x= %f %f %f m= %f h= %f\n", cell->com_legacy[0] , cell->com_legacy[1] , cell->com_legacy[2] , cell->mass_legacy , cell->h[0]);

	  for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
	    dx[k] = cell->com_legacy[k] - parts[i].x[k];
	    r2 += dx[k]*dx[k];
	  }

	  /* Apply the gravitational acceleration. */
	  ir = 1.0 / sqrt( r2 );
	  w = cell->mass_legacy * const_G * ir * ir * ir;
	  for ( k = 0 ; k < 3 ; k++ )
            parts[i].a_legacy[k] += w * dx[k];

	  (*countMultipoles)++;
	  
	}
      else
	{
	  if( parts[i].id == monitor )
	    printf( "Recursing...\n" );

	  /* OK, we need to recurse */
	  for( k = 0 ; k < 8 ; k++ )
	    legacy_interact( parts, i, cell->progeny[k] , monitor , countMultipoles , countPairs );
	}
    }
}
  
  
  
/**
 * @brief Does a tree walk as in the B-H original work for all particles
 *
 * @param N The number of particles
 * @param parts The array of particles
 * @param root The root cell of the tree
 * @param monitor ID of the particle to monitor and output interactions to stdout
 */
void legacy_tree_walk( int N , struct part* parts , struct cell* root , int monitor , int* countMultipoles, int* countPairs , int* countCoMs ) {
  
  int i;

  legacy_comp_com( root , countCoMs );
  
  for( i = 0 ; i < N ; ++i )
    {
      if ( parts[i].id == monitor )
	printf( "tree walk for particle %d x= %f %f %f \n", parts[i].id , parts[i].x[0] , parts[i].x[1], parts[i].x[2] );

      legacy_interact( parts , i , root , monitor , countMultipoles , countPairs );

      if ( parts[i].id == monitor )
	printf( "\n[LEGACY] acceleration for particle %d a= %f %f %f \n", parts[i].id , parts[i].a_legacy[0] , parts[i].a_legacy[1], parts[i].a_legacy[2] );

    }

  
}
  

  






/**
 * @brief Solve the particle interactions using the stupid N^2 algorithm
 *
 * @param N The number of particles
 * @param parts The array of particles
 */
void interact_exact( int N , struct part* parts , int monitor) {

  int i, j, k;
  double ir, w, r2, dx[3] = {0., 0., 0.};

  for ( i = 0 ; i < N ; ++i ) {
      
    for ( j = i+1 ; j < N ; ++j ) {
	
      /* Compute the pairwise distance. */
      for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
	dx[k] = parts[i].x[k] - parts[j].x[k];
	r2 += dx[k]*dx[k];
      }
      
      
      /* Apply the gravitational acceleration. */
      ir = 1.0 / sqrt( r2 );
      w = const_G * ir * ir * ir;

      for ( k = 0 ; k < 3 ; k++ ) {
	parts[j].a_exact[k] += w * dx[k] * parts[i].mass;
	parts[i].a_exact[k] -= w * dx[k] * parts[j].mass;
      }
      
      
    }
  }

  for( i = 0 ; i < N ; ++i )
    if( parts[i].id == monitor )
      	printf( "[EXACT ] acceleration for particle %d a= %f %f %f \n\n", parts[i].id , parts[i].a_exact[0] , parts[i].a_exact[1], parts[i].a_exact[2] );

}



 
/**
 * @brief Set up and run a task-based Barnes-Hutt N-body solver.
 *
 * @param N The number of random particles to use.
 * @param nr_threads Number of threads to use.
 */
 
void test_bh ( int N , int nr_threads , int runs ) {

  int i , k;
    struct cell *root;
    struct part *parts;
    struct qsched s;
    ticks tic, toc_run, tot_setup = 0, tot_run = 0, tic_legacy;
    int countMultipoles, countPairs, countCoMs;

    
    /* Runner function. */
    void runner ( int type , void *data ) {
    
        ticks tic = getticks();
    
        /* Decode the data. */
        struct cell **d = (struct cell **)data;
    
        /* Decode and execute the task. */
        switch ( type ) {
            case task_type_self:
                iact_self( d[0] );
                break;
            case task_type_pair:
                iact_pair( d[0] , d[1] );
                break;
            case task_type_pair_pc:
                iact_pair_pc( d[0] , d[1] );
                break;
            case task_type_com:
                comp_com( d[0] );
                break;
            default:
                error( "Unknown task type." );
            }
            
        atomic_add( &task_timers[type] , getticks() - tic );
        
        }
        
    /* Initialize the per-task type timers. */
    for ( k = 0 ; k < task_type_count ; k++ )
        task_timers[k] = 0;
    
    /* Initialize the scheduler. */
    qsched_init( &s , nr_threads , qsched_flag_noreown );
    
    /* Init and fill the particle array. */
    if ( ( parts = (struct part *)malloc( sizeof(struct part) * N ) ) == NULL )
        error( "Failed to allocate particle buffer." );

    // int hN = ceil( cbrt( N ) );
    // double h = 1.0 / hN;
    for ( k = 0 ; k < N ; k++ ) {
        parts[k].id = k;
        parts[k].x[0] = ((double)rand())/RAND_MAX; // h*(k % hN);
        parts[k].x[1] = ((double)rand())/RAND_MAX; // h*((k / hN) % hN);
        parts[k].x[2] = ((double)rand())/RAND_MAX; // h*(k / hN / hN);
        parts[k].mass = ((double)rand())/RAND_MAX; // 1.0;
        parts[k].a[0] = 0.0;
        parts[k].a[1] = 0.0;
        parts[k].a[2] = 0.0;
        }
        
    /* Init the cells. */
    root = cell_get();
    root->loc[0] = 0.0; root->loc[1] = 0.0; root->loc[2] = 0.0;
    root->h[0] = 1.0; root->h[1] = 1.0; root->h[2] = 1.0;
    root->count = N;
    root->parts = parts;
    cell_split( root , &s );

    printf("----------------------------------------------------------\n");

    /* Do a traditional (serial) B-H tree walk */
    tic_legacy = getticks();
    countMultipoles = 0;
    countPairs = 0;
    countCoMs = 0;
    legacy_tree_walk( N , parts , root , ICHECK , &countMultipoles , &countPairs , &countCoMs );
    printf( "Legacy B-H walk (1 thread) took %lli (= %e) ticks\n", getticks() - tic_legacy , (float)(getticks() - tic_legacy) );
    


    /* Do a N^2 interactions calculation */
    interact_exact( N , parts , ICHECK );

    printf("----------------------------------------------------------\n");

    /* Create the tasks. */
    tic = getticks();
    create_tasks( &s , root , NULL );
    tot_setup += getticks() - tic;

    /* Dump the number of tasks. */
    message( "total nr of tasks: %i." , s.count );
    message( "total nr of deps: %i." , s.count_deps );
    message( "total nr of res: %i." , s.count_res );
    message( "total nr of locks: %i." , s.count_locks );
    message( "total nr of uses: %i." , s.count_uses );
    int counts[ task_type_count ];
    for ( k = 0 ; k < task_type_count ; k++ )
        counts[k] = 0;
    for ( k = 0 ; k < s.count ; k++ )
        counts[ s.tasks[k].type ] += 1;
    printf( "task counts: [ %8s %8s %8s %8s ]\n" , "self", "direct" , "m-poles" , "CoMs" );
    printf( "task counts: [ " );
    for ( k = 0 ; k < task_type_count ; k++ )
        printf( "%8i " , counts[k] );
    printf( "].\n" );

    printf( "task counts: [ %8i %8i %8i %8i ].\n" , 0 , countPairs , countMultipoles , countCoMs );
        
    /* Loop over the number of runs. */
    for ( k = 0 ; k < runs ; k++ ) {

      for ( i = 0 ; i < N ; ++i ) {
      	parts[i].a[0] = 0.;
      	parts[i].a[1] = 0.;
      	parts[i].a[2] = 0.;
      }

        /* Execute the tasks. */
        tic = getticks();
        qsched_run( &s , nr_threads , runner );
	toc_run = getticks();
	message( "%ith run took %lli (= %e) ticks..." , k , toc_run - tic , (float)(toc_run - tic) );
        tot_run += toc_run - tic;
        
        }
        
    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ )
        printf( " %i %i %lli %lli\n" , s.tasks[k].type , s.tasks[k].qid , s.tasks[k].tic , s.tasks[k].toc ); */
        
    /* Dump the costs. */
    message( "costs: setup=%lli ticks, run=%lli ticks." ,
        tot_setup , tot_run/runs );
        
    /* Dump the timers. */
    for ( k = 0 ; k < qsched_timer_count ; k++ )
        message( "timer %s is %lli ticks." , qsched_timer_names[k] , s.timers[k]/runs );
    
    /* Dump the per-task type timers. */
    printf( "task timers: [ " );
    for ( k = 0 ; k < task_type_count ; k++ )
        printf( "%lli " , task_timers[k]/runs );
    printf( "] ticks.\n" );


    /* Dump the particles to a file */
    FILE* file = fopen( "particle_dump.dat" , "w" );
    fprintf(file, "# a_exact.x   a_exact.y    a_exact.z    a_legacy.x    a_legacy.y    a_legacy.z    a_new.x     a_new.y    a_new.z\n");
    for ( k = 0 ; k < N ; ++k )
      fprintf ( file , "%d %e %e %e %e %e %e %e %e %e\n" ,
    		parts[k].id,
    		parts[k].a_exact[0] , parts[k].a_exact[1] , parts[k].a_exact[2] ,
    		parts[k].a_legacy[0] , parts[k].a_legacy[1] , parts[k].a_legacy[2] ,
    		parts[k].a[0] , parts[k].a[1] , parts[k].a[2] );
    fclose( file );
    
    /* Clean up. */
    qsched_free( &s );
    
    }
    

/**
 * @brief Main function.
 */
 
int main ( int argc , char *argv[] ) {

    int c, nr_threads;
    int N = 1000, runs = 1;
    
    /* Get the number of threads. */
    #pragma omp parallel shared(nr_threads)
    {
        if ( omp_get_thread_num() == 0 )
            nr_threads = omp_get_num_threads();
    }
    
    /* Parse the options */
    while ( ( c = getopt( argc , argv  , "n:r:t:" ) ) != -1 )
        switch( c ) {
	        case 'n':
	            if ( sscanf( optarg , "%d" , &N ) != 1 )
	                error( "Error parsing number of particles." );
	            break;
            case 'r':
                if ( sscanf( optarg , "%d" , &runs ) != 1 )
                    error( "Error parsing number of runs." );
                break;
	        case 't':
	            if ( sscanf( optarg , "%d" , &nr_threads ) != 1 )
	                error( "Error parsing number of threads." );
	            omp_set_num_threads( nr_threads );
	            break;
	        case '?':
                fprintf( stderr , "Usage: %s [-t nr_threads] [-n N] [-r runs]\n" , argv[0] );
                fprintf( stderr , "Solves the N-body problem using a Barnes-Hutt\n"
                                  "tree code with N random particles in [0,1]^3 using\n"
                                  "nr_threads threads.\n" );
	            exit( EXIT_FAILURE );
	        }
            
    /* Dump arguments. */
    message( "Computing the N-body problem over %i particles using %i threads (%i runs)." ,
        N , nr_threads , runs );
        
    /* Run the test. */
    test_bh( N , nr_threads, runs );
    
    return 0;
    
    }
    
    

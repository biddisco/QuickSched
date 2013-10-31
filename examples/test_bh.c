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
#define cell_maxparts 100
#define task_limit 4000
#define const_G 6.6738e-11
#define dist_min 0.5


/** Data structure for the particles. */
struct part {
    double x[3];
    double a[3];
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
    int res, com_tid;
    };
    
    
/** Task types. */
enum task_type {
    task_type_self,
    task_type_pair,
    task_type_pair_pc,
    task_type_com
    };
    
    
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
        for ( k = 1 ; k < 100 ; k++ )
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
    if ( root == NULL )
        root = c;
    
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
            while ( i <= count-1 && parts[i].x[0] <= pivot[0] )
                i += 1;
            while ( j >= 0 && parts[j].x[0] > pivot[0] )
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
                while ( i <= right[k] && parts[i].x[1] <= pivot[1] )
                    i += 1;
                while ( j >= left[k] && parts[j].x[1] > pivot[1] )
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
                while ( i <= right[k] && parts[i].x[2] <= pivot[2] )
                    i += 1;
                while ( j >= left[k] && parts[j].x[2] > pivot[2] )
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
        c->com[0] /= c->mass; c->com[1] /= c->mass; c->com[2] /= c->mass;
            
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
        c->com[0] /= c->mass; c->com[1] /= c->mass; c->com[2] /= c->mass;
            
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
    if ( count == 0 )
        return;
    
    /* Init the com's data. */
    for ( k = 0 ; k < 3 ; k++ )
        com[k] = cj->com[k];
    mcom = cj->mass;

    /* Loop over every following particle. */
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

        } /* loop over every other particle. */
                
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
    double xi[3], ai[3], mi, mj, dx[3], r2, ir, w;
    struct part *parts_i = ci->parts, *parts_j = cj->parts;
    
    /* Early abort? */
    if ( count_i == 0 || count_j == 0 )
        return;
    
    /* Recurse? */
    if ( ci->split && cj->split )
        for ( j = 0 ; j < 8 ; j++ )
            for ( k = j+1 ; k < 8 ; k++ )
                iact_pair( ci->progeny[j] , cj->progeny[k] );
            
    /* Get the minimum distance between both cells. */
    for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
        dx[k] = fabs( ci->loc[k] - cj->loc[k] );
        if ( dx[k] > 0 )
            dx[k] -= ci->h[k];
        r2 += dx[k]*dx[k];
        }

    /* Sufficiently well-separated? */
    if ( r2/ci->h[0] > dist_min*dist_min ) {

        /* Compute the center of mass interactions. */
        iact_pair_pc( ci , cj );
        iact_pair_pc( cj , ci );

        }

    /* Recurse? */
    else if ( ci->split && cj->split )
        for ( j = 0 ; j < 8 ; j++ )
            for ( k = j+1 ; k < 8 ; k++ )
                iact_pair( ci->progeny[j] , cj->progeny[k] );
            
    /* Otherwise, do direct interactions. */
    else {

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

                } /* loop over every other particle. */

            /* Store the accumulated acceleration on the ith part. */
            for ( k = 0 ; k < 3 ; k++ )
                parts_i[i].a[k] += ai[k];

            } /* loop over all particles. */

        } /* otherwise, compute interactions directly. */

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
    
    /* Recurse? */
    if ( c->split )
        for ( j = 0 ; j < 8 ; j++ ) {
            iact_self( c->progeny[j] );
            for ( k = j+1 ; k < 8 ; k++ )
                iact_pair( c->progeny[j] , c->progeny[k] );
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
        
            }
    
        }
        
    /* Otherwise, it's a pair. */
    else {
    
        /* Get the minimum distance between both cells. */
        for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
            dx = fabs( ci->loc[k] - cj->loc[k] );
            if ( dx > 0 )
                dx -= ci->h[k];
            r2 += dx*dx;
            }
            
        /* Are the cells sufficiently well separated? */
        if ( r2/ci->h[0] > dist_min*dist_min ) {

            /* Interact ci's parts with cj as a cell. */
            data[0] = ci; data[1] = cj;
            tid = qsched_addtask( s , task_type_pair_pc , task_flag_none , data , sizeof(struct cell *) * 2 , ci->count );
            qsched_addlock( s , tid , ci->res );
            qsched_addunlock( s , cj->com_tid , tid );

            /* Interact cj's parts with ci as a cell. */
            data[0] = cj; data[1] = ci;
            tid = qsched_addtask( s , task_type_pair_pc , task_flag_none , data , sizeof(struct cell *) * 2 , ci->count );
            qsched_addlock( s , tid , cj->res );
            qsched_addunlock( s , ci->com_tid , tid );

            }
            
        /* Does this task need to be broken-down further? */
        else if ( ci->split && cj->split &&
                  ci->count > task_limit && cj->count > task_limit ) {
             
            /* Loop over all pairs between ci and cj's progeny. */
            for ( j = 0 ; j < 8 ; j++ )
                for ( k = 0 ; k < 8 ; k++ )
                    create_tasks( s , ci->progeny[j] , cj->progeny[k] );
        
            }
            
        /* Otherwise, generate a part-part task. */
        else {
        
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
 * @brief Set up and run a task-based Barnes-Hutt N-body solver.
 *
 * @param N The number of random particles to use.
 * @param nr_threads Number of threads to use.
 */
 
void test_bh ( int N , int nr_threads , int runs ) {

    int k;
    struct cell *root;
    struct part *parts;
    struct qsched s;
    ticks tic, toc_run, tot_setup = 0, tot_run = 0;
    
    /* Runner function. */
    void runner ( int type , void *data ) {
    
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
        }
    
    /* Initialize the scheduler. */
    qsched_init( &s , nr_threads , qsched_flag_noreown );
    
    /* Init and fill the particle array. */
    if ( ( parts = (struct part *)malloc( sizeof(struct part) * N ) ) == NULL )
        error( "Failed to allocate particle buffer." );
    for ( k = 0 ; k < N ; k++ ) {
        parts[k].id = k;
        parts[k].x[0] = ((double)rand())/RAND_MAX;
        parts[k].x[1] = ((double)rand())/RAND_MAX;
        parts[k].x[2] = ((double)rand())/RAND_MAX;
        parts[k].mass = ((double)rand())/RAND_MAX;
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

    /* Create the tasks. */
    tic = getticks();
    create_tasks( &s , root , NULL );
    tot_setup += getticks() - tic;

    /* Loop over the number of runs. */
    for ( k = 0 ; k < runs ; k++ ) {
    
        /* Execute the tasks. */
        tic = getticks();
        qsched_run( &s , nr_threads , runner );
	    toc_run = getticks();
	    message( "%ith run took %lli ticks..." , k , toc_run - tic );
        tot_run += toc_run - tic;
        
        }
        
    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ )
        printf( " %i %i %lli %lli\n" , s.tasks[k].type , s.tasks[k].qid , s.tasks[k].tic , s.tasks[k].toc ); */
        
    /* Dump the costs. */
    message( "costs: setup=%lli ticks, run=%lli ticks." ,
        tot_setup , tot_run/runs );
    
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
    
    

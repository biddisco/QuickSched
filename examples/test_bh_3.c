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
#define cell_maxparts 1
#define task_limit 5000
#define const_G 6.6738e-8
#define dist_min 0.5 // 0.5


#define ICHECK -1

/** Data structure for the particles. */
struct part {
    double x[3];
    double a_legacy[3];
    double a_exact[3];
    double mass;
    int id;
    };
   
    
/** Data structure for the BH tree cell. */
struct cell {
    double loc[3];
    double h[1];
    int split, count;
    struct part *parts;
    //struct cell *progeny[8];
    //struct cell *parent;
    struct cell* firstchild; /* Next node if opening */
    struct cell* sibling;  /* Next node */
    double com_legacy[3];
    double mass_legacy;
    int res, com_tid; //, depth; 
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
            cell_pool[k-1].firstchild = &cell_pool[k];
    
        }
        
    /* Pick a cell off the pool. */
    res = cell_pool;
    cell_pool = cell_pool->firstchild;
    
    /* /\* Clean up a few things. *\/ */
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
 
void cell_split ( struct cell *c , struct qsched *s ) {

    int i, j, k, kk, count = c->count;
    struct part temp, *parts = c->parts;
    struct cell *cp;
    int left[8], right[8];
    double pivot[3];
    static struct cell *root = NULL;
    struct cell *progenitors[8];
    
    /* Set the root cell. */
    if ( root == NULL ) {
        root = c;
        //c->depth = 0;
	//c->parent = 0;
	c->sibling = 0;
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
            progenitors[k] = cp = cell_get();
            //cp->parent = c;
            //cp->depth = c->depth + 1;
            cp->loc[0] = c->loc[0];
            cp->loc[1] = c->loc[1];
            cp->loc[2] = c->loc[2];
            cp->h[0] = c->h[0]/2;
            cp->h[0] = c->h[0]/2;
            cp->h[0] = c->h[0]/2;
            cp->res = qsched_addres( s , qsched_owner_none , c->res );
            if ( k & 4 )
                cp->loc[0] += cp->h[0];
            if ( k & 2 )
                cp->loc[1] += cp->h[0];
            if ( k & 1 )
                cp->loc[2] += cp->h[0];
            }
    
        /* Init the pivots. */
        for ( k = 0 ; k < 3 ; k++ )
            pivot[k] = c->loc[k] + c->h[0]/2;

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
	    progenitors[k]->count = right[k] - left[k] + 1;
            progenitors[k]->parts = &c->parts[ left[k] ];
            }


        /* Prepare the pointers. */
        for ( k = 0 ; k < 7 ; k++ ) {

	  /* Find the next non-empty sibling */
	  for ( kk = k+1 ; kk < 8; ++kk ) {
	    if ( progenitors[kk]->count > 0 ) {
	      progenitors[k]->sibling = progenitors[kk];
	      break;
	      }
	  }

	  /* No non-empty sibling ? Go back a level */
	  if( kk == 8 )
	    progenitors[k]->sibling = c->sibling;

	}

	/* Last progenitor links to the next sibling */
	progenitors[7]->sibling = c->sibling;
	c->firstchild = progenitors[0];

        /* Recurse. */
        for ( k = 0 ; k < 8 ; k++ )
            cell_split( progenitors[k] , s );
            
        /* Link the COM tasks. */
        for ( k = 0 ; k < 8 ; k++ )
            qsched_addunlock( s , progenitors[k]->com_tid , c->com_tid );
            
        } /* does the cell need to be split? */
        
    /* Set this cell's resources ownership. */
    qsched_res_own( s , c->res , s->nr_queues * ( c->parts - root->parts ) / root->count );
        
    }
    
    









/**
 * @brief Compute the center of mass of a given cell recursively.
 *
 * @param c The #cell.
 */
 
void legacy_comp_com ( struct cell *c , int* countCoMs ) {

    int k, count = c->count;
    struct cell *cp, *nextsib;
    struct part *p, *parts = c->parts;

    ++(*countCoMs);

    /* Initialize the monopole */
    c->com_legacy[0] = 0.;
    c->com_legacy[1] = 0.;
    c->com_legacy[2] = 0.;
    c->mass_legacy = 0.;

    /* Is the cell split? */
    if ( c->split ) {

      /* Loop over the projenitors */
      cp = c->firstchild;

      while ( cp != c->sibling )
	{  
	    /* Recurse */
	    legacy_comp_com( cp , countCoMs );

	    /* Collect multipole information */
            c->com_legacy[0] += cp->com_legacy[0]*cp->mass_legacy;
            c->com_legacy[1] += cp->com_legacy[1]*cp->mass_legacy;
            c->com_legacy[2] += cp->com_legacy[2]*cp->mass_legacy;
            c->mass_legacy += cp->mass_legacy;

	    /* Move to next child */
	    cp = cp->sibling;
	}

      
        /* Finish multipole calculation */
	if ( c->mass_legacy != 0. )
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

	if ( c->mass_legacy != 0. )
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

void legacy_interact( struct part* parts, int i , struct cell* root , int monitor,  int* countMultipoles, int* countPairs) {
  
  int j,k;
  double r2, dx[3], ir, w;
  struct cell* cell, *currentcell;

  cell = root;

  while( cell != NULL )
    {

      /* Are we in a leaf ? */
      if( !cell->split )
	{

	  /* Interact the particle with the particles in the leaf */    
	  for( j = 0 ; j < cell->count ; ++ j )
	    {
	      if( cell->parts[j].id == parts[i].id )
		continue;
	      
#if ICHECK >= 0
	      if( parts[i].id == monitor )
		printf( "[BH_] Interaction with particle id= %d\n", cell->parts[j].id );
#endif

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


	  cell = cell->sibling;
	}
      else
	{

#if ICHECK >= 0
      if( parts[i].id == monitor )
	printf( "This is a node with %d particles h= %f. r= %f theta= %f\n", cell->count, cell->h[0], sqrt( r2 ), dist_min );
#endif


	  /* We are in a node */	  
	  for ( r2 = 0.0, k = 0 ; k < 3 ; k++ ) {
	    dx[k] = cell->com_legacy[k] - parts[i].x[k];
	    r2 += dx[k]*dx[k];
	  }
      
	  /* Is the cell far enough ? */
	  if( dist_min*dist_min*r2 < cell->h[0]*cell->h[0]  )
	    {

#if ICHECK >= 0
	  if( parts[i].id == monitor )
	    printf( "Recursing...\n" );
#endif
	      cell = cell->firstchild;
	      continue;
	    }


#if ICHECK >= 0
	  if( parts[i].id == monitor )
	    printf( "[BH_] Can interact with the monopole. x= %f %f %f m= %f h= %f\n", cell->com_legacy[0] , cell->com_legacy[1] , cell->com_legacy[2] , cell->mass_legacy , cell->h[0]);
#endif


	  /* Apply the gravitational acceleration. */
	  ir = 1.0 / sqrt( r2 );
	  w = cell->mass_legacy * const_G * ir * ir * ir;
	  for ( k = 0 ; k < 3 ; k++ )
            parts[i].a_legacy[k] += w * dx[k];
	  
	  (*countMultipoles)++;


	  /* Move to the next node */
	  cell = cell->sibling;


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
  struct cell* last = 0;

  /* Compute multipoles (recursively) */
  legacy_comp_com( root , countCoMs );
  
  //#pragma omp parallel for
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
  double ir, w, r2, dx[3];

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
 
void test_bh ( int N , int nr_threads , int runs , char* fileName ) {

  int i , k;
    struct cell *root;
    struct part *parts;
    struct qsched s;
    ticks tic, toc_run, tot_setup = 0, tot_run = 0, tic_exact, toc_exact;
    int countMultipoles, countPairs, countCoMs;

    /* Initialize the scheduler. */
    qsched_init( &s , nr_threads , qsched_flag_noreown );
   
    /* Init and fill the particle array. */
    if ( ( parts = (struct part *)malloc( sizeof(struct part) * N ) ) == NULL )
        error( "Failed to allocate particle buffer." );

    if ( fileName[0] == 0 ) {
      for ( k = 0 ; k < N ; k++ ) {
        parts[k].id = k;
        parts[k].x[0] = ((double)rand())/RAND_MAX; 
        parts[k].x[1] = ((double)rand())/RAND_MAX; 
        parts[k].x[2] = ((double)rand())/RAND_MAX; 
        parts[k].mass = ((double)rand())/RAND_MAX; 
        parts[k].a_legacy[0] = 0.0;
        parts[k].a_legacy[1] = 0.0;
        parts[k].a_legacy[2] = 0.0;
      }
    }
    else {

      FILE* file = fopen( fileName , "r" );
      if ( file ) {
	for ( k = 0 ; k < N ; k++ ) {
	  fscanf( file , "%d" , &parts[k].id );
	  fscanf( file , "%lf" , &parts[k].x[0] ); 
	  fscanf( file , "%lf" , &parts[k].x[1] ); 
	  fscanf( file , "%lf" , &parts[k].x[2] ); 
	  fscanf( file , "%lf" , &parts[k].mass ); 
	  parts[k].a_legacy[0] = 0.0;
	  parts[k].a_legacy[1] = 0.0;
	  parts[k].a_legacy[2] = 0.0;
	}
      }

    }
      


    /* Init the cells. */
    root = cell_get();
    root->loc[0] = 0.0; root->loc[1] = 0.0; root->loc[2] = 0.0;
    root->h[0] = 1.0; root->h[0] = 1.0; root->h[0] = 1.0;
    root->count = N;
    root->parts = parts;
    cell_split( root , &s );

    printf("----------------------------------------------------------\n");

    /* Do a N^2 interactions calculation */

    tic_exact = getticks();
    //interact_exact( N , parts , ICHECK );
    toc_exact = getticks();

    printf( "Exact calculation (1 thread) took %lli (= %e) ticks\n", toc_exact - tic_exact , (float)(toc_exact - tic_exact) );
    
    printf("----------------------------------------------------------\n");

    /* /\* Create the tasks. *\/ */
    /* tic = getticks(); */
    /* create_tasks( &s , root , NULL ); */
    /* tot_setup += getticks() - tic; */

    /* /\* Dump the number of tasks. *\/ */
    /* message( "total nr of tasks: %i." , s.count ); */
    /* message( "total nr of deps: %i." , s.count_deps ); */
    /* message( "total nr of res: %i." , s.count_res ); */
    /* message( "total nr of locks: %i." , s.count_locks ); */
    /* message( "total nr of uses: %i." , s.count_uses ); */
    /* int counts[ task_type_count ]; */
    /* for ( k = 0 ; k < task_type_count ; k++ ) */
    /*     counts[k] = 0; */
    /* for ( k = 0 ; k < s.count ; k++ ) */
    /*     counts[ s.tasks[k].type ] += 1; */
    /* printf( "task counts: [ %8s %8s %8s %8s ]\n" , "self", "direct" , "m-poles" , "CoMs" ); */
    /* printf( "task counts: [ " ); */
    /* for ( k = 0 ; k < task_type_count ; k++ ) */
    /*     printf( "%8i " , counts[k] ); */
    /* printf( "].\n" ); */


    char buffer[200];
    sprintf( buffer, "timings_legacy_%d_%d.dat", cell_maxparts, nr_threads );
    FILE* fileTime = fopen(buffer, "w");
        
    /* Loop over the number of runs. */
    for ( k = 0 ; k < runs ; k++ ) {

	for ( i = 0 ; i < N ; i++ ) {
	  parts[i].a_legacy[0] = 0.0;
	  parts[i].a_legacy[1] = 0.0;
	  parts[i].a_legacy[2] = 0.0;
	}


	countMultipoles = 0;
	countPairs = 0;
	countCoMs = 0;


        /* Execute the legacy walk. */
        tic = getticks();
	legacy_tree_walk( N , parts , root , ICHECK , &countMultipoles , &countPairs , &countCoMs );
	toc_run = getticks();

	message( "%ith run took %lli (= %e) ticks..." , k , toc_run - tic , (float)(toc_run - tic) );
        tot_run += toc_run - tic;
        
	fprintf(fileTime, "%lli %e\n", toc_run - tic, (float)(toc_run - tic) );

        }
        
    fclose(fileTime);


    printf( "task counts: [ %8i %8i %8i %8i ].\n" , 0 , countPairs , countMultipoles , countCoMs );

    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ ) */
    /*     printf( " %i %i %lli %lli\n" , s.tasks[k].type , s.tasks[k].qid , s.tasks[k].tic , s.tasks[k].toc ); */
        
    /* /\* Dump the costs. *\/ */
    /* message( "costs: setup=%lli ticks, run=%lli ticks." , */
    /*     tot_setup , tot_run/runs ); */
        
    /* /\* Dump the timers. *\/ */
    /* for ( k = 0 ; k < qsched_timer_count ; k++ ) */
    /*     message( "timer %s is %lli ticks." , qsched_timer_names[k] , s.timers[k]/runs ); */
    
    /* /\* Dump the per-task type timers. *\/ */
    /* printf( "task timers: [ " ); */
    /* for ( k = 0 ; k < task_type_count ; k++ ) */
    /*     printf( "%lli " , task_timers[k]/runs ); */
    /* printf( "] ticks.\n" ); */


    /* Dump the particles to a file */
    FILE* file = fopen( "particle_dump.dat" , "w" );
    fprintf(file, "# a_exact.x   a_exact.y    a_exact.z    a_legacy.x    a_legacy.y    a_legacy.z    a_new.x     a_new.y    a_new.z\n");
    for ( k = 0 ; k < N ; ++k )
      fprintf ( file , "%d %e %e %e %e %e %e %e %e %e %e %e %e\n" ,
    		parts[k].id,
		parts[k].x[0], parts[k].x[1], parts[k].x[2],
		parts[k].a_exact[0] , parts[k].a_exact[1] , parts[k].a_exact[2] ,
    		parts[k].a_legacy[0] , parts[k].a_legacy[1] , parts[k].a_legacy[2] ,
    		0., 0., 0. );//parts[k].a[0] , parts[k].a[1] , parts[k].a[2] );
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
    char fileName[100] = {0};
    
    /* Get the number of threads. */
    #pragma omp parallel shared(nr_threads)
    {
        if ( omp_get_thread_num() == 0 )
            nr_threads = omp_get_num_threads();
    }
    
    /* Parse the options */
    while ( ( c = getopt( argc , argv  , "n:r:t:f:" ) ) != -1 )
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
      case 'f':
	if ( sscanf( optarg , "%s" , &fileName[0] ) != 1 )
	  error( "Error parsing file name." );
	break;
      case '?':
	fprintf( stderr , "Usage: %s [-t nr_threads] [-n N] [-r runs] [-f file]\n" , argv[0] );
	fprintf( stderr , "Solves the N-body problem using a Barnes-Hutt\n"
		 "tree code with N random particles read from a file in [0,1]^3 using\n"
		 "nr_threads threads.\n" );
	exit( EXIT_FAILURE );
	        }

    /* Tree node information */
    printf( "Size of cell: %zu bytes.\n" , sizeof( struct cell ) );

            
    /* Dump arguments. */
    if (fileName[0] == 0)
      {
	message( "Computing the N-body problem over %i random particles using %i threads (%i runs)." ,
		 N , nr_threads , runs );
      }
    else
      {
	message( "Computing the N-body problem over %i particles read from '%s' using %i threads (%i runs)." ,
		 N , fileName , nr_threads , runs );
      }

        
    /* Run the test. */
    test_bh( N , nr_threads, runs, fileName );
    
    return 0;
    
    }
    
    

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


/* Standard includes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* Local includes. */
#include "quicksched.h"

/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s():%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


/**
 * @brief Matrix multiplication kernel.
 */
 
void matmul ( int m , int n , int k , double *a , int lda , double *b , int ldb , double *c , int ldc ) {

    int ii, jj, kk;
    double acc;
    
    // printf( "matmul: m=%i, n=%i, k=%i, lda=%i, ldb=%i, ldc=%i.\n" ,
    //     m , n , k , lda , ldb , ldc ); fflush(stdout);
    
    for ( ii = 0 ; ii < m ; ii++ )
        for ( jj = 0 ; jj < n ; jj++ ) {
            for ( acc = 0.0, kk = 0 ; kk < k ; kk++ )
                acc += a[ ii + lda*kk ] * b[ kk + ldb*jj ];
            c[ ii + ldc*jj ] += acc;
            }

    }
    
    
/**
 * @brief First test: Just tasks, no dependencies or conflicts.
 *
 * Computes a tiled matrix multiplication of the form
 * C_ij = A_i: * B_:j, with k taskx C_ij += A_ik*B_kj.
 */
 
void test2 ( int m , int n , int k , int nr_threads ) { 

    int i, j, kk, qid, data[3], *d, tid, rid;
    struct sched s;
    struct task *t;
    double *a, *b, *c, *res, err = 0.0, irm = 1.0/RAND_MAX;
    ticks tic_task, toc_task, tic_ref, toc_ref;
    
    /* Init the sched. */
    bzero( &s , sizeof(struct sched) );
    sched_init( &s , nr_threads , m * n );

    /* Allocate the matrices. */
    if ( ( a = (double *)malloc( sizeof(double) * m * k * 32 * 32 ) ) == NULL ||
         ( b = (double *)malloc( sizeof(double) * k * n * 32 * 32 ) ) == NULL ||
         ( c = (double *)malloc( sizeof(double) * m * n * 32 * 32 ) ) == NULL ||
         ( res = (double *)malloc( sizeof(double) * m * n * 32 * 32 ) ) == NULL )
        error( "Failed to allocate matrices." );
        
    /* Fill the matrices. */
    for ( i = 0 ; i < m * k * 32 * 32 ; i++ )
        a[i] = rand() * irm;
    for ( i = 0 ; i < k * n * 32 * 32 ; i++ )
        b[i] = rand() * irm;
    bzero( c , sizeof(double) * m * n * 32 * 32 );
    bzero( res , sizeof(double) * m * n * 32 * 32 );
    
    /* Build a task for each tile of the matrix c. */
    for ( i = 0 ; i < m ; i++ )
        for ( j = 0 ; j < n ; j++ ) {
            rid = sched_addres( &s );
            data[0] = i; data[1] = j;
            for ( kk = 0 ; kk < k ; kk++ ) {
                data[2] = kk;
                tid = sched_newtask( &s , 1 , 0 , 0 , data , 3*sizeof(int) );
                sched_addlock( &s , tid , rid );
                }
            }
            
    /* Prepare the sched for execution. */
    sched_prepare( &s );
            
    /* Parallel loop. */
    tic_task = getticks();
    #pragma omp parallel private(t,qid,d)
    {
    
        /* Get the ID of this runner. */
        if ( ( qid = omp_get_thread_num() ) < nr_threads ) {
    
            /* Main loop. */
            while ( 1 ) {

                /* Get a task, break if unsucessful. */
                if ( ( t = sched_gettask( &s , qid ) ) == NULL )
                    break;

                /* Decode and execute the task. */
                switch ( t->type ) {
                    case 1:
                        d = sched_getdata( &s , t );
                        // printf( "test2[%02i]: working on block [ %i , %i ] with k=%i, lock[0]=%i.\n" , qid , d[0] , d[1] , d[2] , t->locks[0] ); fflush(stdout);
                        matmul( 32 , 32 , 32 , &a[ d[2]*32*m*32 + d[0]*32 ] , m*32 , &b[ k*32*d[1]*32 + d[2]*32 ] , k*32 , &c[ d[0]*32 + m*32*d[1]*32 ] , m*32 );
                        break;
                    default:
                        error( "Unknown task type." );
                    }

                /* Clean up afterwards. */
                sched_done( &s , t );

                } /* main loop. */
                
            } /* valid queue? */
    
        }
    toc_task = getticks();
    
    /* Verify the result. */
    tic_ref = getticks();
    matmul( m*32 , n*32 , k*32 , a , m*32 , b , k*32 , res , m*32 );
    toc_ref = getticks();
    for ( i = 0 ; i < m * n * 32 * 32 ; i++ )
        err += ( res[i] - c[i] ) * ( res[i] - c[i] );
    printf( "test2: Frob. norm of error is %.3e.\n" , sqrt( err ) );
    printf( "test2: tasks took %lli ticks.\n" , toc_task - tic_task );
    printf( "test2: ref.  took %lli ticks.\n" , toc_ref - tic_ref );
    
    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ ) {
        d = (int *)&s.data[ s.tasks[k].data ];
        printf( " %i %i %i %i %lli %lli\n" , k , s.tasks[k].qid , d[0] , d[1] , s.tasks[k].tic , s.tasks[k].toc );
        } */
    
    /* Clean up. */
    sched_free( &s );
    free( a );
    free( b );
    free( c );
    free( res );
    
    }
    
    
/**
 * @brief First test: Just tasks, no dependencies or conflicts.
 *
 * Computes a tiled matrix multiplication of the form
 * C_ij = A_i: * B_:j, with a single task per C_ij.
 */
 
void test1 ( int m , int n , int k , int nr_threads ) { 

    int i, j, qid, data[2], *d, tid, rid;
    struct sched s;
    struct task *t;
    double *a, *b, *c, *res, err = 0.0, irm = 1.0/RAND_MAX;
    ticks tic_task, toc_task, tic_ref, toc_ref;
    
    /* Init the sched. */
    bzero( &s , sizeof(struct sched) );
    sched_init( &s , nr_threads , m * n );

    /* Allocate the matrices. */
    if ( ( a = (double *)malloc( sizeof(double) * m * k * 32 * 32 ) ) == NULL ||
         ( b = (double *)malloc( sizeof(double) * k * n * 32 * 32 ) ) == NULL ||
         ( c = (double *)malloc( sizeof(double) * m * n * 32 * 32 ) ) == NULL ||
         ( res = (double *)malloc( sizeof(double) * m * n * 32 * 32 ) ) == NULL )
        error( "Failed to allocate matrices." );
        
    /* Fill the matrices. */
    for ( i = 0 ; i < m * k * 32 * 32 ; i++ )
        a[i] = rand() * irm;
    for ( i = 0 ; i < k * n * 32 * 32 ; i++ )
        b[i] = rand() * irm;
    bzero( c , sizeof(double) * m * n * 32 * 32 );
    bzero( res , sizeof(double) * m * n * 32 * 32 );
    
    /* Build a task for each tile of the matrix c. */
    for ( i = 0 ; i < m ; i++ )
        for ( j = 0 ; j < n ; j++ ) {
            data[0] = i; data[1] = j;
            rid = sched_addres( &s );
            tid = sched_newtask( &s , 1 , 0 , 0 , data , 2*sizeof(int) );
            sched_addlock( &s , tid , rid );
            }
            
    /* Prepare the sched for execution. */
    sched_prepare( &s );
            
    /* Parallel loop. */
    tic_task = getticks();
    #pragma omp parallel private(t,qid,d)
    {
    
        /* Get the ID of this runner. */
        if ( ( qid = omp_get_thread_num() ) < nr_threads ) {
    
            /* Main loop. */
            while ( 1 ) {

                /* Get a task, break if unsucessful. */
                if ( ( t = sched_gettask( &s , qid ) ) == NULL )
                    break;

                /* Decode and execute the task. */
                switch ( t->type ) {
                    case 1:
                        d = sched_getdata( &s , t );
                        // printf( "test1[%02i]: working on block [ %i , %i ].\n" , qid , d[0] , d[1] ); fflush(stdout);
                        matmul( 32 , 32 , k*32 , &a[ d[0]*32 ] , m*32 , &b[ k*32*d[1]*32 ] , k*32 , &c[ d[0]*32 + m*32*d[1]*32 ] , m*32 );
                        break;
                    default:
                        error( "Unknown task type." );
                    }

                /* Clean up afterwards. */
                sched_done( &s , t );

                } /* main loop. */
                
            } /* valid thread. */
    
        }
    toc_task = getticks();
    
    /* Verify the result. */
    tic_ref = getticks();
    matmul( m*32 , n*32 , k*32 , a , m*32 , b , k*32 , res , m*32 );
    toc_ref = getticks();
    for ( i = 0 ; i < m * n * 32 * 32 ; i++ )
        err += ( res[i] - c[i] ) * ( res[i] - c[i] );
    printf( "test1: Frob. norm of error is %.3e.\n" , sqrt( err ) );
    printf( "test1: tasks took %lli ticks.\n" , toc_task - tic_task );
    printf( "test1: ref.  took %lli ticks.\n" , toc_ref - tic_ref );
    
    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ ) {
        d = (int *)&s.data[ s.tasks[k].data ];
        printf( " %i %i %i %i %lli %lli\n" , k , s.tasks[k].qid , d[0] , d[1] , s.tasks[k].tic , s.tasks[k].toc );
        } */
    
    /* Clean up. */
    sched_free( &s );
    free( a );
    free( b );
    free( c );
    free( res );
    
    }
    
    
/**
 * @brief Main function.
 */
 
int main ( int argc , char *argv[] ) {

    int nr_threads;
    int M = 4, N = 4, K = 4;
    
    /* Get the number of threads. */
    #pragma omp parallel shared(nr_threads)
    {
        #pragma omp single
        nr_threads = omp_get_num_threads();
    }
    
    /* Call the first test. */
    test1( M , N , K , nr_threads );

    /* Call the second test. */
    test2( M , N , K , nr_threads );

    }
    
    

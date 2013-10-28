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


/**
 * @brief Matrix multiplication kernel.
 */
 
void matmul ( int m , int n , int k , double *a , int lda , double *b , int ldb , double *c , int ldc ) {

    int ii, jj, kk;
    double acc;
    
    // message( "matmul: m=%i, n=%i, k=%i, lda=%i, ldb=%i, ldc=%i." ,
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
 * C_ij = A_i: * B_:j, with tasks C_ij += A_ik*B_kj.
 */
 
void test2 ( int m , int n , int k , int nr_threads ) { 

    int i, j, kk, data[3];
    qsched_task_t tid;
    qsched_res_t rid;
    struct qsched s;
    double *a, *b, *c, *res, err = 0.0, irm = 1.0/RAND_MAX;
    ticks tic_task, toc_task, tic_ref, toc_ref;
    
    
    /* Runner function to pass to the scheduler. */
    void runner ( int type , void *data ) {
    
        /* Decode the task data. */
        int *d = (int *)data;
        
        /* Decode and execute the task. */
        switch ( type ) {
            case 1:
                // message( "thread %i working on block [ %i , %i ] with k=%i, lock[0]=%i." , qid , d[0] , d[1] , d[2] , t->locks[0] ); fflush(stdout);
                matmul( 32 , 32 , 32 , &a[ d[2]*32*m*32 + d[0]*32 ] , m*32 , &b[ k*32*d[1]*32 + d[2]*32 ] , k*32 , &c[ d[0]*32 + m*32*d[1]*32 ] , m*32 );
                break;
            default:
                error( "Unknown task type." );
            }

        }
        
    
    
    /* Tell the user something about the test. */
    message( "computing a tiled matrix multiplication of the form "
             "C_ij = A_i: * B_:j, with tasks for each k where C_ij += A_ik*B_kj." );
    
    /* Init the sched. */
    bzero( &s , sizeof(struct qsched) );
    qsched_init( &s , nr_threads, qsched_flag_none );

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
            rid = qsched_addres( &s , -1 );
            data[0] = i; data[1] = j;
            for ( kk = 0 ; kk < k ; kk++ ) {
                data[2] = kk;
                tid = qsched_addtask( &s , 1 , task_flag_none , data , 3*sizeof(int) , 1 );
                qsched_addlock( &s , tid , rid );
                }
            }
            

    /* Run the scheduler. */
    tic_task = getticks();
    qsched_run( &s , nr_threads , runner );
    toc_task = getticks();
            

    /* Verify the result. */
    tic_ref = getticks();
    matmul( m*32 , n*32 , k*32 , a , m*32 , b , k*32 , res , m*32 );
    toc_ref = getticks();
    for ( i = 0 ; i < m * n * 32 * 32 ; i++ )
        err += ( res[i] - c[i] ) * ( res[i] - c[i] );
    message( "Frob. norm of error is %.3e." , sqrt( err ) );
    message( "tasks took %lli ticks." , toc_task - tic_task );
    message( "ref.  took %lli ticks." , toc_ref - tic_ref );
    
    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ ) {
        d = (int *)&s.data[ s.tasks[k].data ];
        printf( " %i %i %i %i %lli %lli\n" , k , s.tasks[k].qid , d[0] , d[1] , s.tasks[k].tic , s.tasks[k].toc );
        } */
    
    /* Clean up. */
    qsched_free( &s );
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

    int i, j, data[2];
    qsched_task_t tid;
    qsched_res_t rid;
    struct qsched s;
    double *a, *b, *c, *res, err = 0.0, irm = 1.0/RAND_MAX;
    ticks tic_task, toc_task, tic_ref, toc_ref;
    
    
    /* Runner function to pass to the scheduler. */
    void runner ( int type , void *data ) {
    
        /* Decode the task data. */
        int *d = (int *)data;
        
        /* Decode and execute the task. */
        switch ( type ) {
            case 1:
                matmul( 32 , 32 , k*32 , &a[ d[0]*32 ] , m*32 , &b[ k*32*d[1]*32 ] , k*32 , &c[ d[0]*32 + m*32*d[1]*32 ] , m*32 );
                break;
            default:
                error( "Unknown task type." );
            }

        }
        
    
    /* Tell the user something about the test. */
    message( "computing a tiled matrix multiplication of the form "
             "C_ij = A_i: * B_:j, with a single task per C_ij." );
    
    /* Init the sched. */
    bzero( &s , sizeof(struct qsched) );
    qsched_init( &s , nr_threads , qsched_flag_none );

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
            rid = qsched_addres( &s , -1 );
            tid = qsched_addtask( &s , 1 , task_flag_none , data , 2*sizeof(int) , 1 );
            qsched_addlock( &s , tid , rid );
            }
            
            
    /* Run the scheduler. */
    tic_task = getticks();
    qsched_run( &s , nr_threads , runner );
    toc_task = getticks();
            

    /* Verify the result. */
    tic_ref = getticks();
    matmul( m*32 , n*32 , k*32 , a , m*32 , b , k*32 , res , m*32 );
    toc_ref = getticks();
    for ( i = 0 ; i < m * n * 32 * 32 ; i++ )
        err += ( res[i] - c[i] ) * ( res[i] - c[i] );
    message( "Frob. norm of error is %.3e." , sqrt( err ) );
    message( "tasks took %lli ticks." , toc_task - tic_task );
    message( "ref.  took %lli ticks." , toc_ref - tic_ref );
    
    /* Dump the tasks. */
    /* for ( k = 0 ; k < s.count ; k++ ) {
        d = (int *)&s.data[ s.tasks[k].data ];
        printf( " %i %i %i %i %lli %lli\n" , k , s.tasks[k].qid , d[0] , d[1] , s.tasks[k].tic , s.tasks[k].toc );
        } */
    
    /* Clean up. */
    qsched_free( &s );
    free( a );
    free( b );
    free( c );
    free( res );
    
    }
    
    
/**
 * @brief Main function.
 */
 
int main ( int argc , char *argv[] ) {

    int c, nr_threads;
    int M = 4, N = 4, K = 4;
    
    /* Get the number of threads. */
    #pragma omp parallel shared(nr_threads)
    {
        if ( omp_get_thread_num() == 0 )
            nr_threads = omp_get_num_threads();
    }
    
    /* Parse the options */
    while ( ( c = getopt( argc , argv  , "m:n:k:t:" ) ) != -1 )
        switch( c ) {
	        case 'k':
	            if ( sscanf( optarg , "%d" , &K ) != 1 )
	                error( "Error parsing dimension M." );
	            break;
	        case 'm':
	            if ( sscanf( optarg , "%d" , &M ) != 1 )
	                error( "Error parsing dimension M." );
	            break;
	        case 'n':
	            if ( sscanf( optarg , "%d" , &N ) != 1 )
	                error( "Error parsing dimension M." );
	            break;
	        case 't':
	            if ( sscanf( optarg , "%d" , &nr_threads ) != 1 )
	                error( "Error parsing number of threads." );
	            omp_set_num_threads( nr_threads );
	            break;
	        case '?':
                fprintf( stderr , "Usage: %s [-t nr_threads] [-m M] [-n N] [-k K]\n" , argv[0] );
                fprintf( stderr , "Computes tests with nr_threads threads for the multiplication\n"
                                  "of a matrix of size MxK and of size KxN tiles of size 32x32.\n" );
	            exit( EXIT_FAILURE );
	        }
            
    /* Dump arguments. */
    message( "multiplying two matrices of size %ix%i and %ix%i using %i threads." ,
        32*M , 32*K , 32*K , 32*N , nr_threads );
    
    /* Call the first test. */
    test1( M , N , K , nr_threads );

    /* Call the second test. */
    test2( M , N , K , nr_threads );

    }
    
    

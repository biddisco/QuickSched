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

/* LAPACKE header. */
#include <lapacke.h>
#include <cblas.h>


/* Local includes. */
#include "quicksched.h"


/*
 * Sam's routines for the tiled QR decomposition.
 */

/*
  \brief Computes 2-norm of a vector \f$x\f$
  
  Computes the 2-norm by computing the following: \f[\textrm{2-norm}=\sqrt_0^lx(i)^2\f]
 */
double do2norm(double* x, int l)
{
	double sum = 0, norm;
	int i;

	for(i = 0; i < l; i++)
		sum += x[i] * x[i];

	norm = sqrt(sum);

	return norm;
}

/**
 * \brief Computes a Householder reflector from a pair of vectors from coupled blocks
 *
 * Calculates the Householder vector of the vector formed by a column in a pair of coupled blocks.
 * There is a single non-zero element, in the first row, of the top vector. This is passed as topDiag
 *
 * \param topDiag The only non-zero element of the incoming vector in the top block
 * \param ma The number of elements in the top vector
 * \param xb Pointer to the lower vector
 * \param l The number of elements in the whole vector
 * \param vk A pointer to a pre-allocated array to store the householder vector of size l
 *
 * \returns void
 */
void calcvkDouble	(double topDiag,
			int ma,
			double* xb,
			int l,
			double* vk)
{
	int sign, i;
	double norm, div;
	//same non-standard normalisation as for single blocks above, but organised without a temporary beta veriable

	sign = topDiag >= 0.0 ? 1 : -1;
	vk[0] = topDiag;
	//use vk[0] as beta
	for(i = 1; i < ma; i++)
		vk[i] = 0;

	for(; i < l; i++)
		vk[i] = xb[i - ma];

	norm = do2norm(vk, l);
	vk[0] += norm * sign;

	if(norm != 0.0)
	{
		div = 1/vk[0];
	
		for(i = 1; i < l; i++)
			vk[i] *= div;
	}
}


void updateDoubleQ_WY	(double* blockA,
			double* blockB,
			double* blockTau,
			int k, int ma, int mb, int n,
			int ldm,
			double* hhVector)//bottom, essential part.
{
	int i, j;

	double tau = 1.0, beta;

	/* Compute tau = 2/v'v */
	for(i = 0; i < mb; i ++)
		tau += hhVector[i] * hhVector[i];

	tau = 2/tau;

	for(j = k; j < n; j ++)
	{
		/* Compute v'*b_j */
		beta = blockA[(j*ldm) + k];

		/* Then for lower half */
		for(i = 0; i < mb; i ++)
			beta += blockB[(j*ldm) + i] * hhVector[i];

		beta *= tau;

		/* Compute b_j = b_j - beta*v_k */
		blockA[(j*ldm) + k] -= beta;
		
		for(i = 0; i < mb; i ++)
			blockB[(j*ldm) + i] -= beta * hhVector[i];
	}

	/* Insert vector below diagonal. */
	for(i = 0; i < mb; i ++)
		blockB[(k*ldm) + i] = hhVector[i];

	blockTau[k] = tau;
}

#pragma omp task in( blockA[0] ) inout( blockB[0] ) in( blockTau[0] )
void DTSQRF	(double* blockA,
		double* blockB,
		double* blockTau,
		int ma,
		int mb,
		int n,
		int ldm )
{
	int k;
	double* xVectA, *xVectB;
    double hhVector[ 2*ma ];
	
	xVectA = blockA;
	xVectB = blockB;

	for(k = 0; k < n; k++)
	{
		//vk = sign(x[1])||x||_2e1 + x
		//vk = vk/vk[0]
		calcvkDouble(xVectA[0], ma - k, xVectB, (ma + mb) - k, hhVector);//returns essential

		//matA(k:ma,k:na) = matA(k:ma,k:na) - (2/(vk.T*vk))*vk*(vk.T*matA(k:ma,k:na)
		//update both blocks, preserving the vectors already stored below the diagonal in the top block and treating them as if they were zeros.
		updateDoubleQ_WY	(blockA, blockB,
					blockTau,
					k, ma, mb, n,
					ldm,
					hhVector + ma - k);

		xVectA += ldm + 1;
		xVectB += ldm;
	}
}

#pragma omp task in( blockV[0] ) in( blockB[0] ) inout( blockA[0] ) in( blockTau[0] )
void DSSRFT	(double* blockV,
		double* blockA, double* blockB,
		double* blockTau,
		int b, int n, int ldm)
{
	int i, j, k;

	double tau, beta;

	/* Compute b_j = b_j - tau*v*v'*b_j for each column j of blocks A & B,
	   and for each householder vector v of blockV */

	/* For each column of B */
	for(j = 0; j < n; j ++)
	{
		/* For each householder vector. */
		for(k = 0; k < n; k ++)
		{
			/* tau = 2/v'v, computed earlier, stored in T(k,k). */
			tau = blockTau[k];

			/* Compute beta = v_k'b_j. */
			/* v_k is >0 (=1) only at position k in top half. */
			beta = blockA[(j*ldm) + k];

			/* For lower portion of v_k, aligning with the lower block */
			for(i = 0; i < b; i ++)
				beta += blockB[(j*ldm) + i] * blockV[(k*ldm) + i];

			beta *= tau;
			
			/* Compute b_j = b_j - beta * v */
			/* v_k = 1 at (k) in top half again */
			blockA[(j*ldm) + k] -= beta;

			/* Apply to bottom block. */
			for(i = 0; i < b; i ++)
				blockB[(j*ldm) + i] -= beta * blockV[(k*ldm) + i];
		}
	}
}


/**
 * @breif Wrapper to get the dependencies right.
 */
 
#pragma omp task inout( a[0] ) inout( tau[0] )
void DGEQRF ( int matrix_order, lapack_int m, lapack_int n,
                                double* a, lapack_int lda, double* tau ) {
    LAPACKE_dgeqrf( matrix_order, m, n, a, lda, tau );
    }


/**
 * @breif Wrapper to get the dependencies right.
 */
 
#pragma omp task inout( v[0] ) in( tau[0] ) inout( t[0] )
void DLARFT ( int matrix_order, char direct, char storev,
                                lapack_int n, lapack_int k, const double* v,
                                lapack_int ldv, const double* tau, double* t,
                                lapack_int ldt ) {
    LAPACKE_dlarft_work( matrix_order, direct, storev, n, k, v, ldv, tau, t, ldt );
    }
			
			
/**
 * @brief Computed a tiled QR factorization using QuickSched.
 *
 * @param m Number of tile rows.
 * @param n Number of tile columns.
 * @param nr_threads Number of threads to use.
 */
 
void test_qr ( int m , int n , int K , int nr_threads , int runs ) {

    int k, j, i, r;
    double *A, *A_orig, *tau;
    ticks tic, toc_run, tot_setup, tot_run = 0;


    /* Allocate and fill the original matrix. */
    if ( ( A = (double *)malloc( sizeof(double) * m * n * K * K ) ) == NULL ||
         ( tau = (double *)malloc( sizeof(double) * m * n * K ) ) == NULL ||
         ( A_orig = (double *)malloc( sizeof(double) * m * n * K * K ) ) == NULL )
        error( "Failed to allocate matrices." );
    for ( k = 0 ; k < m * n * K * K ; k++ )
        A_orig[k] = 2*((double)rand()) / RAND_MAX - 1.0;
    memcpy( A , A_orig , sizeof(double) * m * n * K * K );
    bzero( tau , sizeof(double) * m * n * K );
    
    /* Dump A_orig. */
    /* message( "A_orig = [" );
    for ( k = 0 ; k < m*K ; k++ ) {
        for ( j = 0 ; j < n*K ; j++ )
            printf( "%.3f " , A_orig[ j*m*K + k ] );
        printf( "\n" );
        }
    printf( "];\n" ); */
    
    /* Loop over the number of runs. */
    for ( r = 0 ; r < runs ; r++ ) {
    
        /* Start the clock. */
        tic = getticks();
        
        /* Launch the tasks. */
        for ( k = 0 ; k < m && k < n ; k++ ) {

            /* Add kth corner task. */
            DGEQRF( LAPACK_COL_MAJOR , K, K ,
                            &A[ k*m*K*K + k*K ] , m*K , &tau[ k*m*K + k*K ] );

            /* Add column tasks on kth row. */
            for ( j = k+1 ; j < n ; j++ ) {
                DLARFT( LAPACK_COL_MAJOR , 'F' , 'C' ,
                                K , K , &A[ k*m*K*K + k*K ] ,
                                m*K , &tau[ k*m*K + k*K ] , &A[ j*m*K*K + k*K ] ,
                                m*K );
                }

            /* For each following row... */
            for ( i = k+1 ; i < m ; i++ ) {

                /* Add the row taks for the kth column. */
                DTSQRF( &A[ k*m*K*K + k*K ] , &A[ k*m*K*K + i*K ] , &tau[ k*m*K + i*K ] , K , K , K , K*m );

                /* Add the inner tasks. */
                for ( j = k+1 ; j < n ; j++ ) {
                    DSSRFT(	&A[ k*m*K + i*K ] , &A[ j*m*K*K + k*K ] , &A[ j*m*K*K + i*K ] , &tau[ k*m*K + i*K ] , K , K , K*m );
                    }

                }

            } /* build the tasks. */
    
        /* Collect timers. */
        toc_run = getticks(); 
	    message( "%ith run took %lli ticks..." , r , toc_run - tic );
        tot_run += toc_run - tic;
        
        }
    
        
    /* Dump the costs. */
    message( "costs: setup=%lli ticks, run=%lli ticks." ,
        tot_setup , tot_run/runs );
    
    }


/**
 * @brief Main function.
 */
 
int main ( int argc , char *argv[] ) {

    int c, nr_threads;
    int M = 4, N = 4, runs = 1, K = 32;
    
    /* Get the number of threads. */
    #pragma omp parallel shared(nr_threads)
    {
        if ( omp_get_thread_num() == 0 )
            nr_threads = omp_get_num_threads();
    }
    
    /* Parse the options */
    while ( ( c = getopt( argc , argv  , "m:n:k:r:t:" ) ) != -1 )
        switch( c ) {
	        case 'm':
	            if ( sscanf( optarg , "%d" , &M ) != 1 )
	                error( "Error parsing dimension M." );
	            break;
	        case 'n':
	            if ( sscanf( optarg , "%d" , &N ) != 1 )
	                error( "Error parsing dimension M." );
	            break;
	        case 'k':
	            if ( sscanf( optarg , "%d" , &K ) != 1 )
	                error( "Error parsing tile size." );
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
                fprintf( stderr , "Usage: %s [-t nr_threads] [-m M] [-n N] [-k K]\n" , argv[0] );
                fprintf( stderr , "Computes the tiled QR decomposition of an MxN tiled\n"
                                  "matrix using nr_threads threads.\n" );
	            exit( EXIT_FAILURE );
	        }
            
    /* Dump arguments. */
    message( "Computing the tiled QR decomposition of a %ix%i matrix using %i threads (%i runs)." ,
        32*M , 32*N , nr_threads , runs );
        
    test_qr( M , N , K , nr_threads , runs );
    
    }
    
    

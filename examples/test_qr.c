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
#include <unistd.h>
#include <math.h>
#include <omp.h>

/* Local includes. */
#include "quicksched.h"


/* Prototypes for BLAS function. */
void dtrmm_ ( char * , char * , char * , char * , int * , int * , double * , double * , int * , double * , int * );


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


/**
 * \brief Computes a Householder reflector \f$v\f$ of a vector \f$x\f$ for a single block
 *
 * Computes: \f[v = \textrm{sign}(x_1)||x||_2e_1+x\f]
 * Then does a non-standard normalisation \f$v\f$: \f[v = \frac{v}{v_1}\f]
 * 
 * \param x Pointer to an array containing a column vector to compute the Householder reflector of
 * \param l The number of elements in \f$x\f$
 * \param vk A pointer to an allocated array to store the resulting vector of size l - 1
 * due to the implied 1 as the first element
 *
 * \returns void
 */
void calcvkSingle	(double* x,
			int l,
			double* vk)
{
	int sign, i;
	double norm, div, beta;

	sign = x[0] >= 0.0 ? 1 : -1;
	beta = x[0];
	//copy the values
	for(i = 1; i < l; i++)
		vk[i-1] = x[i];

	//take the euclidian norm of the original vector
	norm = do2norm(x, l);
	//calculate the new normalisation
	beta += norm * sign;

	if(norm != 0.0)
	{
		//normalise 
		div = 1/beta;
	
		for(i = 0; i < l-1; i++)
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

void updatekthSingleWY	(double* blockV,
			double* tauBlock,
			double beta,
			int k,
			int m, int n, int ldm,
			double* w)
{
	/* Insert beta on the diagonal of Tau */
	tauBlock[k] = beta;
}

void updateSingleQ_WY	(double* block,
			double* tauBlock,
			int k,
			int m, int n, int ldm,//dims of block
			double* workVector)
{
	/* Compute A = A - 2/v'v*vv'A */
	int i, j;
	double beta = 1.0f, prod;
	
	for(i = k + 1; i < m; i ++)
	{
		beta += workVector[i - k - 1] * workVector[i - k - 1];
	}
	/* Finish computation of 2/v'v */
	beta = (-2)/beta;
	
	for(j = k; j < 32; j ++)
	{
		/* Compute prod = v'A_j */
		prod = block[(j*ldm) + k];//(k,k) to (k,n)

		for(i = k + 1; i < m; i ++)
			prod += block[(j*ldm) + i] * workVector[i - k - 1];

		/* Compute A_j = A_j - beta*v*prod */
		block[(j*ldm) + k] += beta * prod;

		for(i = k + 1; i < m; i ++)
			block[(j*ldm) + i] += beta * prod * workVector[i - k - 1];
	}

	/* Insert nonessential vector below diagonal. */
	for(i = k + 1; i < m; i ++)
		block[(k*ldm) + i] = workVector[i - k - 1];
	
	updatekthSingleWY	(block,
				tauBlock,
				-beta,
				k, m, n, ldm,
				workVector);
}

void DTSQRF	(double* blockA,
		double* blockB,
		double* blockTau,
		int ma,
		int mb,
		int n,
		int ldm,
		double* hhVector)
{
	int k;
	double* xVectA, *xVectB;
	
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
			
			
void DGEQRF	(double* block,
		double* tauBlock,
		int m, int n, int ldm,
		double* workVector)
{
	int k;
	double* xVect;
	
	xVect = block;

	for(k = 0; k < n; k ++)
	{
		/* Get kth householder vector into position starting at workVector */
		calcvkSingle(xVect, m-k, workVector);

		/* Apply householder vector (with an implied 1 in first element to block,
		   generating WY matrices in the process.
		   Stores vector below the diagonal. */
		updateSingleQ_WY	(block, tauBlock,
					k, m, n, ldm,
					workVector);

		/* Shift one along & one down */
		xVect += ldm + 1;
	}
}

void DLARFT	(double* block,
		double* blockV,
		double* tauBlock,
		int m, int n, int ldm)
{
	/* 	Perform the transformation block = block - blockV*(tauBlock*(blockV^T*block)) 
	 	Equivalent to B = B - V(T(V^TB))
		Noting that T is upper triangular, and V is unit lower triangular. */
	int i, j, k;

	double tau, beta;

	/* For each column of the block. */
	for(j = 0; j < n; j ++)
	{
		/* Apply successive reflectors with b_j - tau_k*v_k*v_k'b_j */
		for(k = 0; k < n; k ++)
		{
			/* tau_k is at blockV(k) */
			tau = tauBlock[k];
	
			/* Compute v_k'*b_j, with v_k,k = 1 implied */
			beta = block[(j*ldm) + k];//*1.0

			/* Rest of vector. */
			for(i = k+1; i < m; i ++)
				beta += blockV[(k*ldm) + i] * block[(j*ldm) + i];

			beta *= tau;

			/* Compute b_j = b_j - beta*v_k, again with an implied 1 at v_kk */
			block[(j*ldm) + k] -= beta;/* *1.0 */
			
			/* Compute for rest of b_j */
			for(i = k+1; i < m; i ++)
				block[(j*ldm) + i] -= beta * blockV[(k*ldm) + i];
		}
	}
}


/**
 * @brief Computed a tiled QR factorization using QuickSched.
 *
 * @param m Number of tile rows.
 * @param n Number of tile columns.
 * @param nr_threads Number of threads to use.
 */
 
void test_qr ( int m , int n , int nr_threads ) {

    int k, j, i;
    double *A, *A_orig, *tau;
    struct qsched s;
    qsched_task_t *tid, tid_new;
    qsched_res_t *rid;
    int data[3];
    
    enum task_types { task_DGEQRF , task_DLARFT , task_DTSQRF , task_DSSRFT };
    
    /* Allocate and fill the original matrix. */
    if ( ( A = (double *)malloc( sizeof(double) * m * n * 32 * 32 ) ) == NULL ||
         ( tau = (double *)malloc( sizeof(double) * m * n * 32 ) ) == NULL ||
         ( A_orig = (double *)malloc( sizeof(double) * m * n * 32 * 32 ) ) == NULL )
        error( "Failed to allocate matrices." );
    for ( k = 0 ; k < m * n * 32 * 32 ; k++ )
        A_orig[k] = 2*((double)rand()) / RAND_MAX - 1.0;
    memcpy( A , A_orig , sizeof(double) * m * n * 32 * 32 );
    bzero( tau , sizeof(double) * m * n * 32 );
    
    /* Dump A_orig. */
    /* message( "A_orig = [" );
    for ( k = 0 ; k < m*32 ; k++ ) {
        for ( j = 0 ; j < n*32 ; j++ )
            printf( "%.3f " , A_orig[ j*m*32 + k ] );
        printf( "\n" );
        }
    printf( "];\n" ); */
    
    /* Initialize the scheduler. */
    qsched_init( &s , nr_threads , m*n );
    
    /* Allocate and init the task ID and resource ID matrix. */
    if ( ( tid = (qsched_task_t *)malloc( sizeof(qsched_task_t) * m * n ) ) == NULL ||
         ( rid = (qsched_res_t *)malloc( sizeof(qsched_res_t) * m * n ) ) == NULL )
        error( "Failed to allocate tid/rid matrix." );
    for ( k = 0 ; k < m * n ; k++ ) {
        tid[k] = -1;
        rid[k] = qsched_addres( &s , -1 );
        }
    
    /* Build the tasks. */
    for ( k = 0 ; k < m && k < n ; k++ ) {
    
        /* Add kth corner task. */
        data[0] = k; data[1] = k; data[2] = k;
        tid_new = qsched_newtask( &s , task_DGEQRF , task_flag_none , data , sizeof(int)*3 , 2 );
        qsched_addlock( &s , tid_new , rid[ k*m + k ] );
        if ( tid[ k*m + k ] != -1 )
            qsched_addunlock( &s , tid[ k*m + k ] , tid_new );
        tid[ k*m + k ] = tid_new;
            
        /* Add column tasks on kth row. */
        for ( j = k+1 ; j < n ; j++ ) {
            data[0] = k; data[1] = j; data[2] = k;
            tid_new = qsched_newtask( &s , task_DLARFT , task_flag_none , data , sizeof(int)*3 , 3 );
            qsched_addlock( &s , tid_new , rid[ j*m + k ] );
            qsched_adduse( &s , tid_new , rid[ k*m + k ] );
            qsched_addunlock( &s , tid[ k*m + k ] , tid_new );
            if ( tid[ j*m + k ] != -1 )
                qsched_addunlock( &s , tid[ j*m + k ] , tid_new );
            tid[ j*m + k ] = tid_new;
            }
            
        /* For each following row... */
        for ( i = k+1 ; i < m ; i++ ) {
        
            /* Add the row taks for the kth column. */
            data[0] = i; data[1] = k; data[2] = k;
            tid_new = qsched_newtask( &s , task_DTSQRF , task_flag_none , data , sizeof(int)*3 , 3 );
            qsched_addlock( &s , tid_new , rid[ k*m + i ] );
            qsched_adduse( &s , tid_new , rid[ k*m + k ] );
            qsched_addunlock( &s , tid[ k*m + (i-1) ] , tid_new );
            if ( tid[ k*m + i ] != -1 )
                qsched_addunlock( &s , tid[ k*m + i ] , tid_new );
            tid[ k*m + i ] = tid_new;
            
            /* Add the inner tasks. */
            for ( j = k+1 ; j < n ; j++ ) {
                data[0] = i; data[1] = j; data[2] = k;
                tid_new = qsched_newtask( &s , task_DSSRFT , task_flag_none , data , sizeof(int)*3 , 5 );
                qsched_addlock( &s , tid_new , rid[ j*m + i ] );
                qsched_adduse( &s , tid_new , rid[ k*m + i ] );
                qsched_adduse( &s , tid_new , rid[ j*m + k ] );
                qsched_addunlock( &s , tid[ k*m + i ] , tid_new );
                qsched_addunlock( &s , tid[ j*m + k ] , tid_new );
                if ( tid[ j*m + i ] != -1 )
                    qsched_addunlock( &s , tid[ j*m + i ] , tid_new );
                tid[ j*m + i ] = tid_new;
                }
        
            }
    
        } /* build the tasks. */
        
    /* Prepare the scheduler. */
    qsched_prepare( &s );

    /* Parallel loop. */
    #pragma omp parallel
    {
    
        int *d, qid;
        double buff[ 2*32*32 ];
        struct task *t;
    
        /* Get the ID of this runner. */
        if ( ( qid = omp_get_thread_num() ) < nr_threads ) {
    
            /* Main loop. */
            while ( 1 ) {

                /* Get a task, break if unsucessful. */
                if ( ( t = qsched_gettask( &s , qid ) ) == NULL )
                    break;
                    
                /* Get the task's data. */
                d = qsched_getdata( &s , t );
                i = d[0]; j = d[1]; k = d[2];

                /* Decode and execute the task. */
                switch ( t->type ) {
                    case task_DGEQRF:
                        DGEQRF( &A[ j*m*32*32 + i*32 ] , &tau[ j*m*32 + i*32 ] , 32 , 32 , 32*m , buff );
                        break;
                    case task_DLARFT:
                        DLARFT( &A[ j*m*32*32 + i*32 ] , &A[ i*m*32*32 + i*32 ] , &tau[ i*m*32 + i*32 ] , 32 , 32 , 32*m );
                        break;
                    case task_DTSQRF:
                        DTSQRF( &A[ j*m*32*32 + j*32 ] , &A[ j*m*32*32 + i*32 ] , &tau[ j*m*32 + i*32 ] , 32 , 32 , 32 , 32*m , buff );
                        break;
                    case task_DSSRFT:
                        DSSRFT(	&A[ k*m*32 + i*32 ] , &A[ j*m*32*32 + k*32 ] , &A[ j*m*32*32 + i*32 ] , &tau[ k*m*32 + i*32 ] , 32 , 32 , 32*m );
                        break;
                    default:
                        error( "Unknown task type." );
                    }

                /* Clean up afterwards. */
                qsched_done( &s , t );

                } /* main loop. */
                
            } /* valid thread. */
    
        } /* parallel loop. */
        
    /* Dump A. */
    /* message( "A = [" );
    for ( k = 0 ; k < m*32 ; k++ ) {
        for ( j = 0 ; j < n*32 ; j++ )
            printf( "%.3f " , A[ j*m*32 + k ] );
        printf( "\n" );
        }
    printf( "];\n" ); */
    
    /* Dump tau. */
    /* message( "tau = [" );
    for ( k = 0 ; k < m*32 ; k++ ) {
        for ( j = 0 ; j < n ; j++ )
            printf( "%.3f " , tau[ j*m*32 + k ] );
        printf( "\n" );
        }
    printf( "];\n" ); */
    
    /* Dump the tasks. */
    for ( k = 0 ; k < s.count ; k++ ) {
        int *d = (int *)&s.data[ s.tasks[k].data ];
        printf( " %i %i %i %i %lli %lli\n" , s.tasks[k].type , s.tasks[k].qid , d[0] , d[1] , s.tasks[k].tic , s.tasks[k].toc );
        }
        
    }


/**
 * @brief Main function.
 */
 
int main ( int argc , char *argv[] ) {

    int c, nr_threads;
    int M = 4, N = 4;
    
    /* Get the number of threads. */
    #pragma omp parallel shared(nr_threads)
    {
        if ( omp_get_thread_num() == 0 )
            nr_threads = omp_get_num_threads();
    }
    
    /* Parse the options */
    while ( ( c = getopt( argc , argv  , "m:n:k:t:" ) ) != -1 )
        switch( c ) {
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
                fprintf( stderr , "Usage: %s [-t nr_threads] [-m M] [-n N]\n" , argv[0] );
                fprintf( stderr , "Computes the tiled QR decomposition of an MxN tiled\n"
                                  "matrix using nr_threads threads.\n" );
	            exit( EXIT_FAILURE );
	        }
            
    /* Dump arguments. */
    message( "Computing the tiled QR decomposition of a %ix%i matrix using %i threads." ,
        32*M , 32*N , nr_threads );
        
    test_qr( M , N , nr_threads );
    
    }
    
    

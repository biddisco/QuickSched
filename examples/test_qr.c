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
* *****************************************************************************/

/* Config parameters. */
#include "config.h"

/* Standard includes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

#include <cblas.h>

/* Local includes. */
#include "quicksched.h"

/**
 * Takes a column major matrix, NOT tile major. size is length of a side of the
 * matrix. Only works for square matrices.
 * This function is simply for validation and is implemented naively as we know
 * of no implementation to retrieve Q from the tiled QR.
 */
double* computeQ(double* HR, int size, int tilesize, double* tau, int tauNum) {
  double* Q = malloc(sizeof(double) * size * size);
  double* Qtemp = malloc(sizeof(double) * size * size);
  double* w = malloc(sizeof(double) * size);
  double* ww = malloc(sizeof(double) * size * size);
  double* temp = malloc(sizeof(double) * size * size);
  int i, k, l, j, n;
  bzero(Q, sizeof(double) * size * size);
  bzero(Qtemp, sizeof(double) * size * size);
  bzero(ww, sizeof(double) * size * size);
  for (i = 0; i < size; i++) {
    Q[i * size + i] = 1.0;
  }
  int numcoltile = size / tilesize;
  int numrowtile = size / tilesize;
  for (k = 0; k < numrowtile; k++) {
    for (l = 0; l < tilesize; l++) {
      bzero(Qtemp, sizeof(double) * size * size);
      for (i = 0; i < size; i++) {
        Qtemp[i * size + i] = 1.0;
      }

      for (i = k; i < numcoltile; i++) {
        bzero(w, sizeof(double) * size);

        for (j = 0; j < tilesize; j++) {
          w[i * tilesize + j] =
              HR[(k * tilesize + l) * size + i * tilesize + j];
        }
        w[k * tilesize + l] = 1.0;
        if (k * tilesize + l > i * tilesize) {
          for (j = 0; j < k * tilesize + l; j++) w[j] = 0.0;
        }

        /* Compute (I - tau*w*w')' */
        for (j = 0; j < size; j++) {
          for (n = 0; n < size; n++) {
            if (j != n)
              ww[n * size + j] =
                  -tau[(k * tilesize + l) * tauNum + i] * w[j] * w[n];
            else
              ww[n * size + j] =
                  1.0 - tau[(k * tilesize + l) * tauNum + i] * w[j] * w[n];
          }
        }

        /* Qtemp = Qtemp * (I-tau*w*w')' */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size,
                    1.0, Qtemp, size, ww, size, 0.0, temp, size);
        double* b = Qtemp;
        Qtemp = temp;
        temp = b;
      }
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size,
                  1.0, Q, size, Qtemp, size, 0.0, temp, size);
      double* b = Q;
      Q = temp;
      temp = b;
    }
  }

  free(Qtemp);
  free(w);
  free(ww);
  free(temp);
  return Q;
}

double* getR(double* HR, int size) {
  double* R = malloc(sizeof(double) * size * size);
  int i, j;
  bzero(R, sizeof(double) * size * size);
  for (i = 0; i < size; i++) {
    for (j = 0; j <= i; j++) {
      R[i * size + j] = HR[i * size + j];
    }
  }
  return R;
}

void printMatrix(double* Matrix, int m, int n, int tilesize) {
  int i, j;

  for (i = 0; i < m * tilesize; i++) {
    for (j = 0; j < n * tilesize; j++) {
      printf(" %.3f ", Matrix[j * m * tilesize + i]);
    }
    printf("\n");
  }
}

double* columnToTile(double* columnMatrix, int size, int m, int n,
                     int tilesize) {
  double* TileMatrix;
  TileMatrix = malloc(sizeof(double) * size);
  if (TileMatrix == NULL) error("failed to allocate TileMatrix");
  int i, j, k, l;

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      double* tileStart =
          &columnMatrix[i * m * tilesize * tilesize + j * tilesize];
      double* tilePos =
          &TileMatrix[i * m * tilesize * tilesize + j * tilesize * tilesize];
      for (k = 0; k < tilesize; k++) {
        tileStart = &columnMatrix[i * m * tilesize * tilesize +
                                  k * m * tilesize + j * tilesize];
        for (l = 0; l < tilesize; l++) {
          tilePos[k * tilesize + l] = tileStart[l];
        }
      }
    }
  }

  return TileMatrix;
}

double* tileToColumn(double* tileMatrix, int size, int m, int n, int tilesize) {
  double* ColumnMatrix;
  ColumnMatrix = (double*)malloc(sizeof(double) * size);
  if (ColumnMatrix == NULL) error("failed to allocate ColumnMatrix");
  int i, j, k, l;
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      /* Tile on ith column is at i*m*32*32.*/
      /* Tile on jth is at j*32*32 */
      double* tile =
          &tileMatrix[i * m * tilesize * tilesize + j * tilesize * tilesize];
      /* Column starts at same position as tile. */
      /* Row j*32.*/
      double* tilePos =
          &ColumnMatrix[i * m * tilesize * tilesize + j * tilesize];
      for (k = 0; k < tilesize; k++) {
        for (l = 0; l < tilesize; l++) {
          tilePos[l] = tile[l];
        }
        /* Next 32 elements are the position of the tile in the next column.*/
        tile = &tile[tilesize];
        /* Move to the j*32th position in the next column. */
        tilePos = &tilePos[tilesize * m];
      }
    }
  }
  return ColumnMatrix;
}

/* Routines for the tiled QR decomposition.*/

/**
 *
 * @brief Computes the QR decomposition of a tile.
 *
 * @param cornerTile A pointer to the tile for which the decomposition is
* computed.
 * @param tilesize The number of elements on a row/column of the tile.
 * @param tauMatrix A pointer to the tau Matrix.
 * @param k The value of k for the QR decomposition
 * @param tauNum The number of tau values stored for each row of the matrix
* (this is equal to the number of tiles on each column).
 *
 *
 */
void DGEQRF(double* restrict cornerTile, int tileSize,
            double* restrict tauMatrix, int k, int tauNum) {
  int i, j, n;
  double norm = 0.0, sign, u1, tau, z;
  double w[tileSize];

  /*Find the householder vector for each row. */
  for (i = 0; i < tileSize; i++) {
    norm = 0.0;
    /*Fill w with the vector.*/
    for (j = i; j < tileSize; j++) {
      /* ith row is i*tileSize, only want elements on diagonal or below.*/
      w[j] = cornerTile[i * tileSize + j];
      /*Find the norm as well */
      norm = norm + w[j] * w[j];
    }
    if (w[i] >= 0.0)
      sign = -1;
    else
      sign = 1;

    norm = sqrt(norm);

    u1 = w[i] - sign * norm;

    if (u1 != 0.0) {
      for (j = i + 1; j < tileSize; j++) w[j] = w[j] / u1;
    } else {
      for (j = i + 1; j < tileSize; j++) w[j] = 0.0;
    }

    if (norm != 0.0)
      tau = -sign * u1 / norm;
    else
      tau = 0.0;

    /*Store the below diagonal vector */

    for (j = i + 1; j < tileSize; j++) {
      cornerTile[i * tileSize + j] = w[j];
    }
    cornerTile[i * tileSize + i] = sign * norm;
    w[i] = 1.0;
    /* Apply the householder transformation to the rest of the tile, for
     * everything to the right of the diagonal. */
    for (j = i + 1; j < tileSize; j++) {
      /*Compute w'*A_j*/
      z = cornerTile[j * tileSize + i];
      for (n = i + 1; n < tileSize; n++) {
        z = z + cornerTile[j * tileSize + n] * w[n];
      }
      /* Tile(m,n) = Tile(m,n) - tau*w[n]* w'A_j */
      for (n = i; n < tileSize; n++) {
        cornerTile[j * tileSize + n] =
            cornerTile[j * tileSize + n] - tau * w[n] * z;
      }
    }
    /* Store tau. We're on k*tileSize+ith row. kth column.*/
    tauMatrix[(k * tileSize + i) * tauNum + k] = tau;
  }
}

/**
 *
 * @brief Applies the householder factorisation of the corner to the row tile.
 *
 * @param cornerTile A pointer to the tile for which the householder is stored.
 * @param rowTiles The tile to which the householder is applied.
 * @param tilesize The number of elements on a row/column of the tile.
 * @param tauMatrix A pointer to the tau Matrix.
 * @param jj The value of j for the QR decomposition.
 * @param kk The value of k for the QR decomposition.
 * @param tauNum The number of tau values stored for each row of the matrix
* (this is equal to the number of tiles on each column).
 *
 *
 */
void DLARFT(double* restrict cornerTile, double* restrict rowTile, int tileSize,
            int jj, int kk, double* restrict tauMatrix, int tauNum) {
  int i, j, n;
  double z = 0.0;
  double w[tileSize];

  /* For each row in the corner Tile*/
  for (i = 0; i < tileSize; i++) {
    /*Get w for row i */
    for (j = i; j < tileSize; j++) {
      w[j] = cornerTile[i * tileSize + j];
    }
    w[i] = 1.0;

    /* Apply to the row Tile */
    for (j = 0; j < tileSize; j++) {
      z = 0.0;
      /* Below Diagonal!*/
      /*Compute w'*A_j*/
      for (n = i; n < tileSize; n++) {
        z = z + w[n] * rowTile[j * tileSize + n];
      }
      for (n = i; n < tileSize; n++) {
        rowTile[j * tileSize + n] =
            rowTile[j * tileSize + n] -
            tauMatrix[(kk * tileSize + i) * tauNum + kk] * w[n] * z;
      }
    }
  }
}

/**
 *
 * @brief Applies the householder factorisation of the corner to the row tile.
 *
 * @param cornerTile The corner tile for this value of k.
 * @param columnTile The tile for which householders are computed.
 * @param tilesize The number of elements on a row/column of the tile.
 * @param tauMatrix A pointer to the tau Matrix.
 * @param ii The value of i for the QR decomposition.
 * @param kk The value of k for the QR decomposition.
 * @param tauNum The number of tau values stored for each row of the matrix
* (this is equal to the number of tiles on each column).
 *
 *
 */
void DTSQRF(double* restrict cornerTile, double* restrict columnTile,
            int tilesize, int ii, int kk, double* restrict tauMatrix,
            int tauNum) {
  int i, j, n;
  double norm = 0.0, sign, u1, tau, z;
  double w[2*tilesize];

  /* For each column compute the householder vector. */
  for (i = 0; i < tilesize; i++) {
    norm = 0.0;
    w[i] = cornerTile[i * tilesize + i];
    norm = norm + w[i] * w[i];
    for (j = i + 1; j < tilesize; j++) {
      w[j] = 0.0;
    }
    for (j = 0; j < tilesize; j++) {
      w[tilesize + j] = columnTile[i * tilesize + j];
      norm = norm + w[tilesize + j] * w[tilesize + j];
    }

    norm = sqrt(norm);
    if (w[i] >= 0.0)
      sign = -1;
    else
      sign = 1;

    u1 = w[i] - sign * norm;
    if (u1 != 0.0) {
      for (j = i + 1; j < 2 * tilesize; j++) {
        w[j] = w[j] / u1;
      }
    } else {
      for (j = i + 1; j < 2 * tilesize; j++) w[j] = 0.0;
    }

    if (norm != 0)
      tau = -sign * u1 / norm;
    else
      tau = 0.0;

    /* Apply to each row to the right.*/
    for (j = i; j < tilesize; j++) {
      /* Find w'*A_j, w is 0s except for first value with upper tile.*/
      z = 1.0 * cornerTile[j * tilesize + i];
      for (n = 0; n < tilesize; n++) {
        z = z + w[tilesize + n] * columnTile[j * tilesize + n];
      }
      /* Apply to upper tile.*/
      cornerTile[j * tilesize + i] =
          cornerTile[j * tilesize + i] - tau * 1.0 * z;
      for (n = i + 1; n < tilesize; n++) {
        cornerTile[j * tilesize + n] =
            cornerTile[j * tilesize + n] - tau * w[n] * z;
      }
      /* Apply to lower tile.*/
      for (n = 0; n < tilesize; n++) {
        columnTile[j * tilesize + n] =
            columnTile[j * tilesize + n] - tau * w[tilesize + n] * z;
      }
    }
    /* Store w*/
    for (j = 0; j < tilesize; j++) {
      columnTile[i * tilesize + j] = w[tilesize + j];
    }
    tauMatrix[(kk * tilesize + i) * tauNum + ii] = tau;
  }
}

/**
 *
 * @brief Applies the householder factorisation of the corner to the row tile.
 *
 * @param cornerTile A pointer to the tile to have the householder applied.
 * @param columnTile The tile in which the householders are stored.
 * @param rowTile The upper tile to have the householders applied.
 * @param tilesize The number of elements on a row/column of the tile.
 * @param tauMatrix A pointer to the tau Matrix.
 * @param ii The value of i for the QR decomposition.
 * @param jj The value of j for the QR decomposition.
 * @param kk The value of k for the QR decomposition.
 * @param tauNum The number of tau values stored for each row of the matrix
* (this is equal to the number of tiles on each column).
 *
 *
 */
void DSSRFT(double* restrict cornerTile, double* restrict columnTile,
            double* restrict rowTile, int tilesize, int ii, int jj, int kk,
            double* restrict tauMatrix, int tauNum) {
  int i, j, n;
  double z;
  double w[2*tilesize];

  for (i = 0; i < tilesize; i++) {
    for (j = 0; j < i; j++) w[j] = 0.0;
    w[i] = 1.0;
    for (j = i + 1; j < tilesize; j++) w[j] = 0.0;
    for (j = 0; j < tilesize; j++)
      w[j + tilesize] = columnTile[i * tilesize + j];

    /* Apply householder vector (w) to the tiles.*/
    for (j = 0; j < tilesize; j++) {
      z = 0.0;
      /* Compute w' * A_j */
      for (n = 0; n < tilesize; n++) {
        z += w[n] * rowTile[j * tilesize + n];
        z += w[n + tilesize] * cornerTile[j * tilesize + n];
      }
      for (n = 0; n < tilesize; n++) {
        rowTile[j * tilesize + n] =
            rowTile[j * tilesize + n] -
            tauMatrix[(kk * tilesize + i) * tauNum + ii] * w[n] * z;
        cornerTile[j * tilesize + n] =
            cornerTile[j * tilesize + n] -
            tauMatrix[(kk * tilesize + i) * tauNum + ii] * w[tilesize + n] * z;
      }
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

void test_qr(int m, int n, int K, int nr_threads, int runs, double* matrix) {

  int k, j, i;
  double* A, *A_orig, *tau;
  struct qsched s;
  qsched_task_t* tid, tid_new;
  qsched_res_t* rid;
  int data[3];
  ticks tic, toc_run, tot_setup, tot_run = 0;

  enum task_types {
    task_DGEQRF,
    task_DLARFT,
    task_DTSQRF,
    task_DSSRFT
  };

  /* Runner function to pass to the scheduler. */
  void runner(int type, void * data) {

    /* Decode the task data. */
    int* idata = (int*)data;
    int i = idata[0], j = idata[1], k = idata[2];

    /* Decode and execute the task. */
    switch (type) {
      case task_DGEQRF:
        DGEQRF(&A[(k * m + k) * K * K], K, tau, k, m);
        break;
      case task_DLARFT:
        DLARFT(&A[(k * m + k) * K * K], &A[(j * m + k) * K * K], K, j, k, tau,
               m);
        break;
      case task_DTSQRF:
        DTSQRF(&A[(k * m + k) * K * K], &A[(k * m + i) * K * K], K, i, k, tau,
               m);
        break;
      case task_DSSRFT:
        DSSRFT(&A[(j * m + i) * K * K], &A[(k * m + i) * K * K],
               &A[(j * m + k) * K * K], K, i, j, k, tau, m);
        break;
      default:
        error("Unknown task type.");
    }
  }

  /* Allocate and fill the original matrix. */
  if ((A = (double*)malloc(sizeof(double)* m* n* K* K)) == NULL ||
      (tau = (double*)malloc(sizeof(double)* m* n* K)) == NULL ||
      (A_orig = (double*)malloc(sizeof(double) * m * n * K * K)) == NULL)
    error("Failed to allocate matrices.");
  for (k = 0; k < m * n * K * K; k++) {
    if (matrix == NULL)
      A_orig[k] = 2 * ((double)rand()) / RAND_MAX - 1.0;
    else
      A_orig[k] = matrix[k];
  }
  memcpy(A, A_orig, sizeof(double) * m * n * K * K);
  bzero(tau, sizeof(double) * m * n * K);

  /* Dump A_orig. */
  /* message( "A_orig = [" );
  for ( k = 0 ; k < m*K ; k++ ) {
      for ( j = 0 ; j < n*K ; j++ )
          printf( "%.3f " , A_orig[ j*m*K + k ] );
      printf( "\n" );
      }
  printf( "];\n" ); */

  /* Initialize the scheduler. */
  qsched_init(&s, nr_threads, qsched_flag_none);

  /* Allocate and init the task ID and resource ID matrix. */
  tic = getticks();
  if ((tid = (qsched_task_t*)malloc(sizeof(qsched_task_t)* m* n)) == NULL ||
      (rid = (qsched_res_t*)malloc(sizeof(qsched_res_t) * m * n)) == NULL)
    error("Failed to allocate tid/rid matrix.");
  for (k = 0; k < m * n; k++) {
    tid[k] = qsched_task_none;
    rid[k] = qsched_addres(&s, qsched_owner_none, qsched_res_none);
  }

  /* Build the tasks. */
  for (k = 0; k < m && k < n; k++) {

    /* Add kth corner task. */
    data[0] = k;
    data[1] = k;
    data[2] = k;
    tid_new = qsched_addtask(&s, task_DGEQRF, task_flag_none, data,
                             sizeof(int) * 3, 2);
    qsched_addlock(&s, tid_new, rid[k * m + k]);
    if (tid[k * m + k] != -1) qsched_addunlock(&s, tid[k * m + k], tid_new);
    tid[k * m + k] = tid_new;

    /* Add column tasks on kth row. */
    for (j = k + 1; j < n; j++) {
      data[0] = k;
      data[1] = j;
      data[2] = k;
      tid_new = qsched_addtask(&s, task_DLARFT, task_flag_none, data,
                               sizeof(int) * 3, 3);
      qsched_addlock(&s, tid_new, rid[j * m + k]);
      qsched_adduse(&s, tid_new, rid[k * m + k]);
      qsched_addunlock(&s, tid[k * m + k], tid_new);
      if (tid[j * m + k] != -1) qsched_addunlock(&s, tid[j * m + k], tid_new);
      tid[j * m + k] = tid_new;
    }

    /* For each following row... */
    for (i = k + 1; i < m; i++) {

      /* Add the row taks for the kth column. */
      data[0] = i;
      data[1] = k;
      data[2] = k;
      tid_new = qsched_addtask(&s, task_DTSQRF, task_flag_none, data,
                               sizeof(int) * 3, 3);
      qsched_addlock(&s, tid_new, rid[k * m + i]);
      qsched_adduse(&s, tid_new, rid[k * m + k]);
      qsched_addunlock(&s, tid[k * m + (i - 1)], tid_new);
      if (tid[k * m + i] != -1) qsched_addunlock(&s, tid[k * m + i], tid_new);
      tid[k * m + i] = tid_new;

      /* Add the inner tasks. */
      for (j = k + 1; j < n; j++) {
        data[0] = i;
        data[1] = j;
        data[2] = k;
        tid_new = qsched_addtask(&s, task_DSSRFT, task_flag_none, data,
                                 sizeof(int) * 3, 5);
        qsched_addlock(&s, tid_new, rid[j * m + i]);
        qsched_adduse(&s, tid_new, rid[k * m + i]);
        qsched_addlock(&s, tid_new, rid[j * m + k]);
        qsched_addunlock(&s, tid[k * m + i], tid_new);
        qsched_addunlock(&s, tid[j * m + i - 1], tid_new);
        if (tid[j * m + i] != -1) qsched_addunlock(&s, tid[j * m + i], tid_new);

        tid[j * m + i] = tid_new;
      }
    }

  } /* build the tasks. */
  tot_setup = getticks() - tic;

  /* Dump the number of tasks. */
  message("total nr of tasks: %i.", s.count);
  message("total nr of deps: %i.", s.count_deps);
  message("total nr of res: %i.", s.count_res);
  message("total nr of locks: %i.", s.count_locks);
  message("total nr of uses: %i.", s.count_uses);

  /* Loop over the number of runs. */
  for (k = 0; k < runs; k++) {

    /* Execute the the tasks. */
    tic = getticks();
    qsched_run(&s, nr_threads, runner);
    toc_run = getticks();
    message("%ith run took %lli ticks...", k, toc_run - tic);
    tot_run += toc_run - tic;
  }

  /* Dump A. */
  /* message( "A = [" );
  for ( k = 0 ; k < m*K ; k++ ) {
      for ( j = 0 ; j < n*K ; j++ )
          printf( "%.3f " , A[ j*m*K + k ] );
      printf( "\n" );
      }
  printf( "];\n" ); */

  /* Dump tau. */
  /* message( "tau = [" );
  for ( k = 0 ; k < m*K ; k++ ) {
      for ( j = 0 ; j < n ; j++ )
          printf( "%.3f " , tau[ j*m*K + k ] );
      printf( "\n" );
      }
  printf( "];\n" ); */

  /* Dump the tasks. */
  /* for ( k = 0 ; k < s.count ; k++ ) {
      int *d = (int *)&s.data[ s.tasks[k].data ];
      printf( " %i %i %i %i %lli %lli\n" , s.tasks[k].type , s.tasks[k].qid ,
     d[0] , d[1] , s.tasks[k].tic , s.tasks[k].toc );
      } */

  /* Dump the costs. */
  message("costs: setup=%lli ticks, run=%lli ticks.", tot_setup,
          tot_run / runs);

  /* Dump the timers. */
  for (k = 0; k < qsched_timer_count; k++)
    message("timer %s is %lli ticks.", qsched_timer_names[k],
            s.timers[k] / runs);

  if (matrix != NULL) {
    for (k = 0; k < m * n * K * K; k++) matrix[k] = A[k];
  }

  /* Test if the decomposition was correct.*/
  /*double *tempMatrix = tileToColumn(A, m*n*K*K, m, n, K);
   double *Q = computeQ(tempMatrix, m*K, K, tau, m);
   double *R = getR(tempMatrix, m*K);
   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m*K, m*K, m*K, 1.0, Q,
   m*K, R, m*K, 0.0, tempMatrix, m*K);
       free(Q);
       Q = tileToColumn(A_orig, m*n*K*K, m, n, K);
   for(i = 0; i < m * n * K * K; i++)
   {
       if(Q[i] != 0 && (Q[i] / tempMatrix[i] > 1.005 || Q[i] / tempMatrix[i] <
   0.995))
           printf("Not correct at value %i %.3f %.3e %.3e\n", i, A[i], Q[i],
   tempMatrix[i]);
   }
   free(tempMatrix);
   free(Q);
   free(R);*/

  /* Clean up. */
  free(A);
  free(A_orig);
  free(tau);
  free(tid);
  free(rid);
  qsched_free(&s);
}

/* Generates a random matrix. */
double* generateColumnMatrix(int size) {
  double* matrix = malloc(sizeof(double) * size * size);
  if (matrix == NULL) error("Failed to allocate matrix");

  unsigned long int m_z = 35532;
  int i;
  for (i = 0; i < size * size; i++) {
    m_z = (1664525 * m_z + 1013904223) % 4294967296;
    matrix[i] = m_z % 100;
    if (matrix[i] < 0) matrix[i] += 100;
  }
  return matrix;
}

/**
 * @brief Main function.
 */

int main(int argc, char* argv[]) {

  int c, nr_threads;
  int M = 4, N = 4, runs = 1, K = 32;

/* Get the number of threads. */
#pragma omp parallel shared(nr_threads)
  {
    if (omp_get_thread_num() == 0) nr_threads = omp_get_num_threads();
  }

  /* Parse the options */
  while ((c = getopt(argc, argv, "m:n:k:r:t:")) != -1) switch (c) {
      case 'm':
        if (sscanf(optarg, "%d", &M) != 1) error("Error parsing dimension M.");
        break;
      case 'n':
        if (sscanf(optarg, "%d", &N) != 1) error("Error parsing dimension M.");
        break;
      case 'k':
        if (sscanf(optarg, "%d", &K) != 1) error("Error parsing tile size.");
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
      case '?':
        fprintf(stderr, "Usage: %s [-t nr_threads] [-m M] [-n N] [-k K]\n",
                argv[0]);
        fprintf(stderr, "Computes the tiled QR decomposition of an MxN tiled\n"
                        "matrix using nr_threads threads.\n");
        exit(EXIT_FAILURE);
    }

  /* Dump arguments. */
  message("Computing the tiled QR decomposition of a %ix%i matrix using %i "
          "threads (%i runs).",
          32 * M, 32 * N, nr_threads, runs);

  test_qr(M, N, K, nr_threads, runs, NULL);
}

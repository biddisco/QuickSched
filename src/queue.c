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

/* Local includes. */
#include "error.h"
#include "cycle.h"
#include "atomic.h"
#include "lock.h"
#include "task.h"
#include "qsched.h"
#include "queue.h"


/**
 * @brief Get a task index from the given #queue.
 *
 * @param q The #queue.
 * @param s The #sched in which this queue's tasks lives.
 *
 * @return The task ID or -1 if no available task could be found.
 */
 
int queue_get ( struct queue *q , struct qsched *s ) {

    int k, j, temp, tid, *inds, w, count;
    struct task *tasks = s->tasks;

    /* Should we even try? */
    if ( q->count == 0 )
        return -1;
        
    /* Lock this queue. */
    if ( lock_lock( &q->lock ) != 0 )
        error( "Failed to lock queue." );
        
    /* Get a pointer to the indices. */
    inds = q->inds;
    count = q->count;
        
    /* Loop over the queue entries. */
    for ( k = 0 ; k < count ; k++ ) {
    
        /* Get the task ID. */
        tid = inds[k];
        
        /* If the task can be locked, break. */
        if ( qsched_locktask( s , tid ) )
            break;
    
        }
        
    /* Did we get a task? */
    if ( k < count ) {
    
        /* Swap the last element to the new heap position. */
        q->count = ( count -= 1 );
        inds[k] = inds[ count ];
    
        /* Fix the heap. */
        w = tasks[ inds[k] ].weight;
        while ( 1 ) {
            if ( ( j = 2*k + 1 ) >= count )
                break;
            if ( j+1 < count && tasks[ inds[j] ].weight < tasks[ inds[j+1] ].weight )
                j = j + 1;
            if ( tasks[ inds[j] ].weight > w ) {
                temp = inds[j];
                inds[j] = inds[k];
                inds[k] = temp;
                k = j;
                }
            else
                break;
            }
        
        } /* did we get a task? */
        
    /* Otherwise, clear the task ID. */
    else
        tid = -1;
        
    /* Unlock the queue. */
    lock_unlock_blind( &q->lock );

    /* Return the task ID. */
    return tid;

    }


/**
 * @brief Add a task index to the given #queue.
 * 
 * @param q The #queue.
 * @param s The #sched in which the tasks live.
 * @param tid The task index.
 */
 
void queue_put ( struct queue *q , struct qsched *s , int tid ) {

    int ind, j, temp;
    struct task *tasks = s->tasks;
    int *inds, *inds_new, w;
    
    /* Lock this queue. */
    if ( lock_lock( &q->lock ) != 0 )
        error( "Failed to lock queue." );
        
    /* Get a pointer to the indices. */
    inds = q->inds;
        
    /* Get the index of the new task. */
    ind = q->count;
    
    /* Does the queue need to be extended? */
    if ( ind >= q->size ) {
    
        /* Increase the queue size. */
        q->size *= queue_grow;
        
        /* Allocate the new indices. */
        if ( ( inds_new = (int *)malloc( sizeof(int) * q->size ) ) == NULL )
            error( "Failed to allocate new indices." );
            
        /* Copy the old indices. */
        memcpy( inds_new , inds , sizeof(int) * q->count );
        
        /* Clear the old indices and replace them with the new. */
        free( inds );
        q->inds = ( inds = inds_new );
    
        }
        
    /* Store the task index. */
    q->count += 1;
    inds[ind] = tid;
    
    /* Bubble up the new entry. */
    w = tasks[ inds[ind] ].weight;
    while ( ind > 0 ) {
        j = (ind - 1)/2;
        if ( tasks[ inds[j] ].weight < w ) {
            temp = inds[j];
            inds[j] = inds[ind];
            inds[ind] = temp;
            ind = j;
            }
        else
            break;
        }
        
    /* Unlock the queue. */
    lock_unlock_blind( &q->lock );

    }


/**
 * @brief Clean up a queue and free its memory.
 */
 
void queue_free ( struct queue *q ) {

    /* Free the inds. */
    if ( q->inds != NULL )
        free( (void *)q->inds );

    }


/** 
 * @brief Initialize the given #queue.
 *
 * @param q The #queue.
 * @param size The maximum size of the queue.
 */
 
void queue_init ( struct queue *q , int size ) {
    
    /* Allocate the task list if needed. */
    if ( q->inds == NULL || q->size < size ) {
        if ( q->inds != NULL )
            free( (int *)q->inds );
        q->size = size;
        if ( ( q->inds = (int *)malloc( sizeof(int) * size ) ) == NULL )
            error( "Failed to allocate queue inds." );
        }
    q->size = size;
    
    /* Init the lock. */
    lock_init( &q->lock );
        
    /* Init the count. */
    q->count = 0;
    
    }


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

/* Local includes. */
#include "atomic.h"
#include "lock.h"
#include "task.h"
#include "queue.h"
#include "sched.h"

/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s():%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


/**
 * @brief Get a task index from the given #queue.
 *
 * @param q The #queue.
 *
 * @return The task ID or -1 if the queue is empty.
 *
 * This function will block until there is either an available
 * index or the queue is empty.
 */
 
int queue_get ( struct queue *q ) {

    int ind, tid = -1;

    /* Should we even try? */
    if ( q->count == 0 )
        return -1;
        
    /* Get the next index. */
    ind = atomic_inc( &q->first ) % q->size;
    
    /* Wait for either there to be something at that index, or
       for the queue to be empty. */
    while ( q->count && ( tid = q->inds[ind] ) != -1 );
    
    /* If we got something, clear the field. */
    if ( tid >= 0 )
        q->inds[ind] = -1;
    
    /* Return whatever it is whe got. */
    return tid;

    }


/**
 * @brief Add a task index to the given #queue.
 * 
 * @param q The #queue.
 * @param tid The task index.
 */
 
void queue_put ( struct queue *q , int tid ) {

    int ind;
    
    /* Let everybody know a task will be inserted. */
    atomic_inc( &q->count );
    
    /* Get the next free index. */
    ind = atomic_inc( &q->last ) % q->size;
    
    /* Wait for the data at that position to be free. */
    while ( q->inds[ind] < 0 );
    
    /* Drop the task index in there. */
    q->inds[ind] = tid;

    }


/** 
 * @brief Initialize the given #queue.
 *
 * @param q The #queue.
 * @param size The maximum size of the queue.
 */
 
void queue_init ( struct queue *q , int size ) {

    int k;
    
    /* Make size a power of two. */
    size -= 1;
    for ( k = 1 ; k < sizeof(int)*8 ; k *= 2 );
        size |= size >> k;
    size += 1;
    
    /* Allocate the task list if needed. */
    if ( q->inds == NULL || q->size < size ) {
        if ( q->inds != NULL )
            free( (int *)q->inds );
        q->size = size;
        if ( ( q->inds = (int *)malloc( sizeof(int) * size ) ) == NULL )
            error( "Failed to allocate queue inds." );
        }
    q->size = size;
        
    /* Fill the list with -1. */
    for ( k = 0 ; k < size ; k++ )
        q->inds[k] = -1;
        
    /* Init counters. */
    q->count = 0;
    q->first = 0;
    q->last = 0;
    
    }


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


/* The queue data structure. */
struct queue {

    /* Index of the first task ID. */
    volatile unsigned int first;
    
    /* Index of the last task ID. */
    volatile unsigned int last;
    
    /* Task indices. */
    volatile int *inds;
    
    /* Number of tasks waiting in the queue. */
    volatile int count;

    /* Maximum number of tasks in queue. */
    int size;

    };


/* Function prototypes. */
int queue_get ( struct queue *q );
void queue_put ( struct queue *q , int tid );
void queue_init ( struct queue *q , int size );
void queue_free ( struct queue *q );

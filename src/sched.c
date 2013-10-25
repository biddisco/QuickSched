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
#include "cycle.h"
#include "atomic.h"
#include "error.h"
#include "lock.h"
#include "task.h"
#include "res.h"
#include "sched.h"
#include "queue.h"

/* Max macro. */
#define max(a,b) ( ((a) >= (b)) ? (a) : (b) )


/**
 * @brief Fetch the data pointer of a task.
 *
 * @param s Pointer to the #sched.
 * @param t Pointer to the #task.
 */
 
void *sched_getdata( struct sched *s , struct task *t ) {

    return &s->data[ t->data ];
    
    }


/**
 * @brief Put the given task in the best possible queue.
 *
 * @param s Pointer to the #sched.
 * @param t Pointer to the #task.
 */
 
void sched_enqueue ( struct sched *s , struct task *t ) {

    int j, qid, scores[ s->nr_queues ];
    
    /* If this is a virtual task, just do its unlocks and leave. */
    if ( t->flags & task_flag_virtual ) {
        
        /* This task is done before it started. */
        t->tic = getticks();
        sched_done( s , t );
        
        }
        
    /* Otherwise, find a home (queue) for it. */
    else {
    
        /* Init the scores for each queue. */
        for ( j = 0 ; j < s->nr_queues ; j++ )
            scores[j] = 0;

        /* Loop over the locks and uses, and get their owners. */
        for ( j = 0 ; j < t->nr_locks ; j++ )
            scores[ s->res[ t->locks[j] ].owner ] += 1;
        for ( j = 0 ; j < t->nr_uses ; j++ )
            scores[ s->res[ t->uses[j] ].owner ] += 1;

        /* Find the queue with the highest score. */
        qid = 0;
        for ( j = 1 ; j < s->nr_queues ; j++ )
            if ( scores[j] > scores[qid] ||
                 ( scores[j] == scores[qid] && s->queues[j].count < s->queues[qid].count ) )
                qid = j;
                
        /* Put the unlocked task in that queue. */
        queue_put( &s->queues[qid] , s , t - s->tasks );
        
        }
    
    }


/**
 * @brief Tell the #sched that a task has completed.
 *
 * @param s Pointer to the #sched.
 * @param t Pointer to the completed #task.
 */
 
void sched_done ( struct sched *s , struct task *t ) {

    int k;
    struct task *t2;
    
    /* Release this task's locks. */
    for ( k = 0 ; k < t->nr_locks ; k++ )
        sched_unlockres( s , t->locks[k] );
    
    /* Loop over the task's unlocks... */
    for ( k = 0 ; k < t->nr_unlocks ; k++ ) {
    
        /* Get a grip on the unlocked task. */
        t2 = &s->tasks[ t->unlocks[k] ];

        /* Is the unlocked task ready to run? */
        if ( atomic_dec( &t2->wait ) == 1 && !( t2->flags & task_flag_skip ) )
            sched_enqueue( s , t2 );
            
        }
        
    /* Set the task stats. */
    t->toc = getticks();
    t->cost = t->toc - t->tic;

    }
    
    
/**
 * @brief Lock a resource and hold its parents.
 *
 * @param s Pointer to the #sched.
 * @param rid The ID of the resource to lock.
 *
 * @return @c 1 if the resource could be locked, @c 0 otherwise.
 */
 
int sched_lockres ( struct sched *s , int rid ) {

    int finger, finger2;
    
    /* Try to lock the root-level resource. */
    if ( s->res[rid].hold || lock_trylock( &s->res[rid].lock ) )
        return 0;
        
    /* Did the resource get held in the meantime? */
    if ( s->res[rid].hold ) {
        lock_unlock_blind( &s->res[rid].lock );
        return 0;
        }
        
    /* Follow parents and increase their hold counter, but fail
       if any are locked. */
    for ( finger = s->res[rid].parent ; finger != sched_res_none ; finger = s->res[finger].parent ) {
        if ( lock_trylock( &s->res[finger].lock ) )
            break;
        atomic_inc( &s->res[finger].hold );
        lock_unlock_blind( &s->res[finger].lock );
        }
        
    /* Did we fail on the way up? */
    if ( finger != sched_res_none ) {
    
        /* Unlock the resource. */
        lock_unlock_blind( &s->res[rid].lock );
    
        /* Go back up the tree and undo the holds. */
        for ( finger2 = s->res[rid].parent ; finger2 != finger ; finger2 = s->res[finger2].parent )
            atomic_dec( &s->res[finger2].hold );
            
        /* Fail. */
        return 0;
    
        }
        
    /* Otherwise, all went well. */
    else
        return 1;

    }


/**
 * @brief Unlock a resource and un-hold its parents.
 *
 * @param s Pointer to the #sched.
 * @param rid The ID of the resource to lock.
 */
 
void sched_unlockres ( struct sched *s , int rid ) {

    int finger;
        
    /* Unlock the resource. */
    lock_unlock_blind( &s->res[rid].lock );

    /* Go back up the tree and undo the holds. */
    for ( finger = s->res[rid].parent ; finger != sched_res_none ; finger = s->res[finger].parent )
        atomic_dec( &s->res[finger].hold );
            
    }
    
    
/**
 * @brief Try to get all the locks for a task.
 * 
 * @param s Pointer to the #sched.
 * @param tid The ID of the #task to lock.
 *
 * @return @c 1 if the resources could be locked, @c 0 otherwise.
 */
 
int sched_locktask ( struct sched *s , int tid ) {

    int k;
    struct task *t;

    /* Get a pointer on the task. */
    t = &s->tasks[tid];
        
    /* Try to lock all the task's locks. */
    for ( k = 0 ; k < t->nr_locks ; k++ )
        if ( sched_lockres( s , t->locks[k] ) == 0 )
            break;

    /* If I didn't get all the locks... */
    if ( k < t->nr_locks ) {

        /* Unroll the locks I got. */
        for ( k -= 1 ; k >= 0 ; k-- )
            sched_unlockres( s , t->locks[k] );

        /* Fail. */
        return 0;

        }
        
    /* Otherwise, all went well. */
    else
        return 1;
            
    }
    

/**
 * @brief Unlock the resources associated with a task.
 *
 * @param s Pointer to the #sched.
 * @param tid The ID of the #task to unlock.
 */
 
void sched_unlocktask ( struct sched *s , int tid ) {

    int k;
    struct task *t;

    /* Get a pointer on the task. */
    t = &s->tasks[tid];
        
    /* Unlock the used resources. */
    for ( k = 0 ; k < t->nr_locks ; k++ )
        sched_unlockres( s , t->locks[k] );

    }


/**
 * @brief Get a task from the #sched.
 *
 * @param s Pointer to the #sched.
 * @param qid The queue to use.
 *
 * @return A pointer to a task object.
 *
 * Note that the #sched has to have been prepared with #sched_prepare
 * before any tasks can be extracted. Adding dependencies or locks
 * will require the #sched to be re-prepared.
 */
 
struct task *sched_gettask ( struct sched *s , int qid ) {

    int k, maxq, tid;
    struct task *t;

    /* Check if the sched is ok. */
    if ( s->flags & sched_flag_dirty || !(s->flags & sched_flag_ready) )
        error( "Calling gettask with dirty or unprepared sched." );
        
    /* Check if the queue ID is ok. */
    if ( qid < 0 || qid >= s->nr_queues )
        error( "Invalid queue ID." );
        
    /* Main loop. */
    while ( s->waiting ) {
    
        /* Init some things. */
        maxq = qid;
    
        /* Try to get a task from my own queue. */
        if ( ( tid = queue_get( &s->queues[qid] , s ) ) < 0 ) {
            
            /* Otherwise, look for the largest queue. */
            maxq = 0;
            for ( k = 0 ; k < s->nr_queues ; k++ )
                if ( k != qid && s->queues[k].count > s->queues[maxq].count )
                    maxq = k;
            if ( ( tid = queue_get( &s->queues[maxq] , s ) ) < 0 )
                continue;
                
            }
            
        /* Bail if a valid task ID was returned. */
        if ( tid >= 0 ) {
        
            /* Get a pointer to the task. */
            t = &s->tasks[tid];
        
            /* Decrease the number of tasks in this space. */
            atomic_dec( &s->waiting );
            
            /* Own the resources. */
            for ( k = 0 ; k < t->nr_uses ; k++ )
                s->res[ t->uses[k] ].owner = qid;
            
            /* Set some stats data. */
            t->tic = getticks();
            t->qid = qid;
            
            /* Return the task. */
            return t;
            
            }
    
        }
        
    /* Return empty-handed. */
    return NULL;

    }


/**
 * @brief Sort the data according to the given indices.
 *
 * @param data The data to be sorted
 * @param ind The indices with respect to which the data are sorted.
 * @param N The number of entries
 * @param min Lowest index.
 * @param max highest index.
 *
 * This function calls itself recursively.
 */
 
void sched_sort ( int *restrict data , int *restrict ind , int N , int min , int max ) {

    int pivot = (min + max) / 2;
    int i = 0, j = N-1;
    int temp_i, temp_d;
    
    /* If N is small enough, just do insert sort. */
    if ( N < 16 ) {
    
        for ( i = 1 ; i < N ; i++ )
            if ( ind[i] < ind[i-1] ) {
                temp_i = ind[i];
                temp_d = data[i];
                for ( j = i ; j > 0 && ind[j-1] > temp_i ; j-- ) {
                    ind[j] = ind[j-1];
                    data[j] = data[j-1];
                    }
                ind[j] = temp_i;
                data[j] = temp_d;
                }
    
        }
        
    /* Otherwise, recurse with Quicksort. */
    else {
    
        /* One pass of quicksort. */
        while ( i < j ) {
            while ( i < N && ind[i] <= pivot )
                i++;
            while ( j >= 0 && ind[j] > pivot )
                j--;
            if ( i < j ) {
                temp_i = ind[i]; ind[i] = ind[j]; ind[j] = temp_i;
                temp_d = data[i]; data[i] = data[j]; data[j] = temp_d;
                }
            }
            
        /* Recurse in parallel? */
        if ( N > 100 ) {

            /* Recurse on the left? */
            if ( j > 0  && pivot > min ) {
                #pragma omp task untied
                sched_sort( data , ind , j+1 , min , pivot );
                }

            /* Recurse on the right? */
            if ( i < N && pivot+1 < max ) {
                #pragma omp task untied
                sched_sort( &data[i], &ind[i], N-i , pivot+1 , max );
                }

            }
        else {
            
            /* Recurse on the left? */
            if ( j > 0  && pivot > min )
                sched_sort( data , ind , j+1 , min , pivot );

            /* Recurse on the right? */
            if ( i < N && pivot+1 < max )
                sched_sort( &data[i], &ind[i], N-i , pivot+1 , max );
            
            }
            
        }
        
    }


/**
 * @brief Prepare a #sched for execution.
 * 
 * @param s Pointer to the #sched.
 */
 
void sched_prepare ( struct sched *s ) {

    int j, k, count;
    struct task *t, *tasks;

    /* Lock the sched. */
    lock_lock( &s->lock );
    
    /* Get a pointer to the tasks, set the count. */
    tasks = s->tasks;
    count = s->count;
    
    /* If the sched is dirty... */
    if ( s->flags & sched_flag_dirty ) {
    
        /* Do the sorts in parallel, if possible. */
        #pragma omp parallel
        {
    
            /* Sort the unlocks. */
            #pragma omp single nowait
            sched_sort( s->deps , s->deps_key , s->count_deps , 0 , count - 1 );

            /* Sort the locks. */
            #pragma omp single nowait
            sched_sort( s->locks , s->locks_key , s->count_locks , 0 , count - 1 );

            /* Sort the uses. */
            #pragma omp single nowait
            sched_sort( s->uses , s->uses_key , s->count_uses , 0 , count - 1 );
            
        }
        
        /* Run throught the tasks and link the locks and unlocks. */
        tasks[0].unlocks = s->deps;
        tasks[0].locks = s->locks;
        tasks[0].uses = s->uses;
        for ( k = 1 ; k < count ; k++ ) {
            tasks[k].unlocks = &tasks[k-1].unlocks[ tasks[k-1].nr_unlocks ];
            tasks[k].locks = &tasks[k-1].locks[ tasks[k-1].nr_locks ];
            tasks[k].uses = &tasks[k-1].uses[ tasks[k-1].nr_uses ];
            }
        
        /* All cleaned-up now! */
        s->flags &= ~sched_flag_dirty;
    
        }
        
    /* Init the queues. */
    for ( k = 0 ; k < s->nr_queues ; k++ )
        queue_init( &s->queues[k] , count );
    
    /* Run through the tasks and set the waits... */
    for ( k = 0 ; k < count ; k++ ) {
        t = &tasks[k];
        if ( !( t->flags & task_flag_skip ) )
            for ( j = 0 ; j < t->nr_unlocks ; j++ )
                tasks[ t->unlocks[j] ].wait += 1;
        }
        
    /* Sort the tasks topologically. */
    int *tid = (int *)malloc( sizeof(int) * count );
    for ( j = 0 , k = 0 ; k < count ; k++ )
        if ( tasks[k].wait == 0 ) {
            tid[j] = k;
            j += 1;
            }
    int ready = j;
    for ( k = 0 ; k < j ; k++ ) {
        t = &tasks[ tid[k] ];
        for ( int kk = 0 ; kk < t->nr_unlocks ; kk++ )
            if ( ( tasks[ t->unlocks[kk] ].wait -= 1 ) == 0 ) {
                tid[j] = t->unlocks[kk];
                j += 1;
                }
        }
    if ( k < count )
        error( "Circular dependencies detected." );
        
    /* Run through the topologically sorted tasks backwards and
       set their weights, re-setting the waits while we're at it. */
    for ( k = count-1 ; k >= 0 ; k-- ) {
        int maxweight = 0;
        t = &tasks[ tid[k] ];
        for ( j = 0 ; j < t->nr_unlocks ; j++ ) {
            tasks[ t->unlocks[j] ].wait += 1;
            if ( tasks[ t->unlocks[j] ].weight > maxweight )
                maxweight = tasks[ t->unlocks[j] ].weight;
            }
        t->weight = t->cost + maxweight;
        }

    /* Run through the tasks and enqueue the non-waiting ones. */
    for ( k = 0 ; k < ready ; k++ ) {
        t = &tasks[tid[k]];
        if ( t->wait == 0 && !( t->flags & task_flag_skip ) )
            sched_enqueue( s , t );
        }
        
    /* Clean up. */
    free( tid );
        
    /* Set the number of waiting tasks. */
    s->waiting = count;
        
    /* Set the ready flag. */
    s->flags |= sched_flag_ready;

    /* Unlock the sched. */
    lock_unlock_blind( &s->lock );

    }


/**
 * @brief Add a new resource to the #sched.
 *
 * @param s Pointer to the #sched.
 * @param parent ID of the parent resource or #sched_res_none if none.
 *
 * @return The ID of the new shared resource.
 */
 
int sched_addres ( struct sched *s , int parent ) {

    void *temp1, *temp2;
    int id;

    /* Lock the sched. */
    lock_lock( &s->lock );
    
    /* Do the deps need to be re-allocated? */
    if ( s->count_res == s->size_res ) {
    
        /* Scale the res list size. */
        s->size_res *= max( s->size_res + sched_inc , s->size_res * sched_stretch );
        
        /* Allocate a new task list. */
        if ( ( temp1 = malloc( sizeof(lock_type) * s->size_res ) ) == NULL ||
             ( temp2 = malloc( sizeof(int) * s->size_res ) ) == NULL )
            error( "Failed to allocate new res lists." );
            
        /* Copy the res and owners over to the new list. */
        memcpy( temp1 , (void *)s->res , sizeof(struct res) * s->count_res );
        
        /* Free the old res lists. */
        free( (int *)s->res );
        
        /* Set the new res lists. */
        s->res = (struct res *)temp1;
    
        }
        
    /* Increase the res counter. */
    id = s->count_res;
    s->count_res += 1;
    
    /* Init the resource. */
    lock_init( &s->res[ id ].lock );
    s->res[ id ].hold = 0;
    s->res[ id ].owner = -1;
    s->res[ id ].parent = parent;
        
    /* Unlock the sched. */
    lock_unlock_blind( &s->lock );

    /* Return the res ID. */
    return id;

    }


/**
 * @brief Add a resource requirement to a task.
 *
 * @param s Pointer to the #sched.
 * @param t ID of the task.
 * @param res ID of the resource.
 */
 
void sched_addlock ( struct sched *s , int t , int res ) {

    void *temp1, *temp2;

    /* Lock the sched. */
    lock_lock( &s->lock );
    
    /* Do the deps need to be re-allocated? */
    if ( s->count_locks == s->size_locks ) {
    
        /* Scale the locks list size. */
        s->size_locks *= max( s->size_locks + sched_inc , s->size_locks * sched_stretch );
        
        /* Allocate a new task list. */
        if ( ( temp1 = malloc( sizeof(int) * s->size_locks ) ) == NULL ||
             ( temp2 = malloc( sizeof(int) * s->size_locks ) ) == NULL )
            error( "Failed to allocate new locks lists." );
            
        /* Copy the locks and keys over to the new list. */
        memcpy( temp1 , s->locks , sizeof(int) * s->count_locks );
        memcpy( temp2 , s->locks_key , sizeof(int) * s->count_locks );
        
        /* Free the old locks lists. */
        free( s->locks );
        free( s->locks_key );
        
        /* Set the new locks lists. */
        s->locks = (int *)temp1;
        s->locks_key = (int *)temp2;
    
        }
        
    /* Add the new dependency. */
    s->locks[ s->count_locks ] = res;
    s->locks_key[ s->count_locks ] = t;
    s->tasks[t].nr_locks += 1;
    
    /* Increase the locks counter. */
    s->count_locks += 1;
    
    /* The sched is now dirty. */
    s->flags |= sched_flag_dirty;
    
    /* Unlock the sched. */
    lock_unlock_blind( &s->lock );

    }


/**
 * @brief Add a resource use to a task.
 *
 * @param s Pointer to the #sched.
 * @param t ID of the task.
 * @param res ID of the resource.
 */
 
void sched_adduse ( struct sched *s , int t , int res ) {

    void *temp1, *temp2;

    /* Lock the sched. */
    lock_lock( &s->lock );
    
    /* Do the deps need to be re-allocated? */
    if ( s->count_uses == s->size_uses ) {
    
        /* Scale the uses list size. */
        s->size_uses *= max( s->size_uses + sched_inc , s->size_uses * sched_stretch );
        
        /* Allocate a new task list. */
        if ( ( temp1 = malloc( sizeof(int) * s->size_uses ) ) == NULL ||
             ( temp2 = malloc( sizeof(int) * s->size_uses ) ) == NULL )
            error( "Failed to allocate new uses lists." );
            
        /* Copy the uses and keys over to the new list. */
        memcpy( temp1 , s->uses , sizeof(int) * s->count_uses );
        memcpy( temp2 , s->uses_key , sizeof(int) * s->count_uses );
        
        /* Free the old uses lists. */
        free( s->uses );
        free( s->uses_key );
        
        /* Set the new uses lists. */
        s->uses = (int *)temp1;
        s->uses_key = (int *)temp2;
    
        }
        
    /* Add the new dependency. */
    s->uses[ s->count_uses ] = res;
    s->uses_key[ s->count_uses ] = t;
    s->tasks[t].nr_uses += 1;
    
    /* Increase the uses counter. */
    s->count_uses += 1;
    
    /* The sched is now dirty. */
    s->flags |= sched_flag_dirty;
    
    /* Unlock the sched. */
    lock_unlock_blind( &s->lock );

    }


/**
 * @brief Add a task dependency.
 *
 * @param s Pointer to the #sched.
 * @param ta ID of the unlocking task.
 * @param tb ID of the unlocked task.
 *
 * A dependency is added such that @c tb depends on @c ta.
 */
 
void sched_addunlock ( struct sched *s , int ta , int tb ) {

    void *temp1, *temp2;

    /* Lock the sched. */
    lock_lock( &s->lock );
    
    /* Do the deps need to be re-allocated? */
    if ( s->count_deps == s->size_deps ) {
    
        /* Scale the deps list size. */
        s->size_deps *= max( s->size_deps + sched_inc , s->size_deps * sched_stretch );
        
        /* Allocate a new task list. */
        if ( ( temp1 = malloc( sizeof(int) * s->size_deps ) ) == NULL ||
             ( temp2 = malloc( sizeof(int) * s->size_deps ) ) == NULL )
            error( "Failed to allocate new deps lists." );
            
        /* Copy the deps and keys over to the new list. */
        memcpy( temp1 , s->deps , sizeof(int) * s->count_deps );
        memcpy( temp2 , s->deps_key , sizeof(int) * s->count_deps );
        
        /* Free the old deps lists. */
        free( s->deps );
        free( s->deps_key );
        
        /* Set the new deps lists. */
        s->deps = (int *)temp1;
        s->deps_key = (int *)temp2;
    
        }
        
    /* Add the new dependency. */
    s->deps[ s->count_deps ] = tb;
    s->deps_key[ s->count_deps ] = ta;
    s->tasks[ta].nr_unlocks += 1;
    
    /* Increase the deps counter. */
    s->count_deps += 1;
    
    /* The sched is now dirty. */
    s->flags |= sched_flag_dirty;
    
    /* Unlock the sched. */
    lock_unlock_blind( &s->lock );

    }


/**
 * @brief Add a new task to the #sched.
 *
 * @param s Pointer to the #sched
 * @param type Task type.
 * @param subtype Task subtype.
 * @param flags Task flags.
 * @param data Pointer to the task data.
 * @param data_size Size, in bytes, of the task data.
 * @param cost Approximate cost for this task.
 */
 
int sched_newtask ( struct sched *s , int type , int subtype , unsigned int flags , void *data , int data_size , int cost ) {

    void *temp;
    struct task *t;
    int id, data_size2;

    /* Lock the sched. */
    lock_lock( &s->lock );
    
    /* Do the tasks need to be re-allocated? */
    if ( s->count == s->size ) {
    
        /* Scale the task list size. */
        s->size *= max( s->size + sched_inc , s->size * sched_stretch );
        
        /* Allocate a new task list. */
        if ( ( temp = malloc( sizeof(struct task) * s->size ) ) == NULL )
            error( "Failed to allocate new task list." );
            
        /* Copy the tasks over to the new list. */
        memcpy( temp , s->tasks , sizeof(struct task) * s->count );
        
        /* Free the old task list. */
        free( s->tasks );
        
        /* Set the new task list. */
        s->tasks = (struct task *)temp;
    
        }
        
    /* Round-up the data size. */
    data_size2 = ( data_size + (sched_data_round-1) ) & ~(sched_data_round-1);
        
    /* Do the task data need to be re-allocated? */
    if ( s->count_data + data_size2 > s->size_data ) {
    
        /* Scale the task list size. */
        s->size_data *= max( s->size + data_size2 , s->size_data * sched_stretch );
        
        /* Allocate a new task list. */
        if ( ( temp = malloc( s->size_data ) ) == NULL )
            error( "Failed to allocate new task list." );
            
        /* Copy the tasks over to the new list. */
        memcpy( temp , s->data , s->count_data );
        
        /* Free the old task list. */
        free( s->data );
        
        /* Set the new task list. */
        s->data = temp;
    
        }
        
    /* Store the new task ID. */
    id = s->count;
        
    /* Init the new task. */
    t = &s->tasks[ id ];
    t->type = type;
    t->subtype = subtype;
    t->flags = flags;
    t->cost = cost;
    t->wait = 0;
    t->nr_conflicts = 0;
    t->nr_unlocks = 0;
    t->nr_locks = 0;
    t->nr_uses = 0;
    
    /* Add a relative pointer to the data. */
    memcpy( &s->data[ s->count_data ] , data , data_size );
    t->data = &s->data[ s->count_data ] - s->data;
    s->count_data += data_size2;
    
    /* Increase the task counter. */
    s->count += 1;
    
    /* Unlock the sched. */
    lock_unlock_blind( &s->lock );
    
    /* Return the task ID. */
    return id;

    }
    
    
/**
 * @brief Clean up a #sched, free all associated memory.
 *
 * @param s Pointer to the #sched.
 */
 
void sched_free ( struct sched *s ) {

    int k;

    /* Clear all the buffers if allocated. */
    if ( s->tasks != NULL ) { free( s->tasks ); s->tasks = NULL; }
    if ( s->deps != NULL ) { free( s->deps ); s->deps = NULL; }
    if ( s->deps_key != NULL ) { free( s->deps_key ); s->deps_key = NULL; }
    if ( s->locks != NULL ) { free( s->locks ); s->locks = NULL; }
    if ( s->locks_key != NULL ) { free( s->locks_key ); s->locks_key = NULL; }
    if ( s->uses != NULL ) { free( s->uses ); s->uses = NULL; }
    if ( s->uses_key != NULL ) { free( s->uses_key ); s->uses_key = NULL; }
    if ( s->res != NULL ) { free( (void *)s->res ); s->res = NULL; }
    if ( s->data != NULL ) { free( s->data ); s->data = NULL; }
    
    /* Loop over the queues and free them too. */
    for ( k = 0 ; k < s->nr_queues ; k++ )
        queue_free( &s->queues[k] );
    free( s->queues );
    s->queues = NULL;
        
    /* Clear the flags. */
    s->flags = sched_flag_none;

    }


/**
 * @brief Initialize the given #sched object.
 *
 * @param s Pointer to a #sched object.
 * @param nr_queues The number of queues in the #sched.
 * @param size The initial number of tasks in the queue.
 *
 * Initializes the given #sched with the given number of queues.
 * The initial size is not a fixed maximum, i.e. the #sched
 * will re-alloate its buffers if more tasks are added.
 */
 
void sched_init ( struct sched *s , int nr_queues , int size ) {
    
    /* Set the flags to begin with. */
    s->flags = sched_flag_none;
    
    /* Allocate and clear the queues (will init when sched is
       finalized. */
    if ( ( s->queues = (struct queue *)malloc( sizeof(struct queue) * nr_queues ) ) == NULL )
        error( "Failed to allocate memory for queues." );
    bzero( s->queues , sizeof(struct queue) * nr_queues );
    s->nr_queues = nr_queues;
    
    /* Allocate the task list. */
    if ( ( s->tasks = (struct task *)malloc( sizeof(struct task) * size ) ) == NULL )
        error( "Failed to allocate memory for tasks." );
    s->size = size;
    s->count = 0;
    
    /* Allocate the initial deps. */
    s->size_deps = sched_init_depspertask * size;
    if ( ( s->deps = (int *)malloc( sizeof(int) * s->size_deps ) ) == NULL ||
         ( s->deps_key = (int *)malloc( sizeof(int) * s->size_deps ) ) == NULL )
        error( "Failed to allocate memory for deps." );
    s->count_deps = 0;

    /* Allocate the initial locks. */
    s->size_locks = sched_init_lockspertask * size;
    if ( ( s->locks = (int *)malloc( sizeof(int) * s->size_locks ) ) == NULL ||
         ( s->locks_key = (int *)malloc( sizeof(int) * s->size_locks ) ) == NULL )
        error( "Failed to allocate memory for locks." );
    s->count_locks = 0;
    
    /* Allocate the initial res. */
    s->size_res = sched_init_respertask * size;
    if ( ( s->res = (struct res *)malloc( sizeof(struct res) * s->size_res ) ) == NULL )
        error( "Failed to allocate memory for res." );
    s->count_res = 0;
    
    /* Allocate the initial uses. */
    s->size_uses = sched_init_usespertask * size;
    if ( ( s->uses = (int *)malloc( sizeof(int) * s->size_uses ) ) == NULL ||
         ( s->uses_key = (int *)malloc( sizeof(int) * s->size_uses ) ) == NULL )
        error( "Failed to allocate memory for uses." );
    s->count_uses = 0;
    
    /* Allocate the initial data. */
    s->size_data = sched_init_datapertask * size;
    if ( ( s->data = malloc( sizeof(int) * s->size_data ) ) == NULL )
        error( "Failed to allocate memory for data." );
    s->count_data = 0;
    
    /* Init the sched lock. */
    lock_init( &s->lock );

    }

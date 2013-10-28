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

/* Scheduler flags. */
#define qsched_flag_none                 0
#define qsched_flag_dirty                1
#define qsched_flag_ready                2
#define qsched_flag_yield                4
#define qsched_flag_pthread              8

/* Some sched-specific constants. */
#define qsched_stretch                   2
#define qsched_size_init                 1000
#define qsched_init_depspertask          2
#define qsched_init_lockspertask         2
#define qsched_init_usespertask          2
#define qsched_init_respertask           2
#define qsched_init_datapertask          8
#define qsched_data_round                16
#define qsched_res_none                  (-1)
#define qsched_task_none                 (-1)


/** Type used for task identifiers. */
typedef int qsched_task_t;


/** Type used for resource identifiers. */
typedef int qsched_res_t;


/** Type definition for the execution function in #qsched_run. */
typedef void (*qsched_funtype)( int , void * );


/* The sched data structre. */
struct qsched {

    /* Flags for this scheduler. */
    unsigned int flags;

    /* The list of tasks in this scheduler. */
    struct task *tasks;
    
    /* The dependency indices. */
    int *deps, *deps_key;
    
    /* The conflict/lock array. */
    int *locks, *locks_key;
    
    /* The conflict/use array. */
    int *uses, *uses_key;
    
    /* The shared resources. */
    // lock_type *res;
    // int *res_owner;
    struct res *res;
    
    /* The task data. */
    char *data;
    
    /* The size of the data buffer. */
    int size_data;
    
    /* The number of data bytes used. */
    int count_data;
    
    /* Number of tasks in sched. */
    int count, waiting;
    
    /* The queues associated with this scheduler. */
    struct queue *queues;
    
    /* Number of queues in the scheduler. */
    int nr_queues;
    
    /* Size of the task array. */
    int size;
    
    /* Size of the dependencies array. */
    int size_deps;
    
    /* Total number of dependencies. */
    int count_deps;
    
    /* Size of the locks array. */
    int size_locks;
    
    /* Total number of locks. */
    int count_locks;
    
    /* Size of the uses array. */
    int size_uses;
    
    /* Total number of uses. */
    int count_uses;
    
    /* Size of the res array. */
    int size_res;
    
    /* Total number of res. */
    int count_res;
    
    /* A lock for the sched itself. */
    lock_type lock;
    
    /* Pthread stuff for condition variable for yielding threads. */
    #ifdef HAVE_PTHREAD
        pthread_cond_t cond;
        pthread_mutex_t mutex;
    #endif
    
    };


/* Function prototypes. */
/* Internal functions. */
void qsched_sort ( int *restrict data , int *restrict ind , int N , int min , int max );
void qsched_sort_rec ( int *restrict data , int *restrict ind , int N , int min , int max );
struct task *qsched_gettask ( struct qsched *s , int qid );
void qsched_done ( struct qsched *s , struct task *t );
void *qsched_getdata( struct qsched *s , struct task *t );
int qsched_lockres ( struct qsched *s , int rid );
void qsched_unlockres ( struct qsched *s , int rid );
int qsched_locktask ( struct qsched *s , int tid );
void qsched_unlocktask ( struct qsched *s , int tid );
void qsched_prepare ( struct qsched *s );

/* External functions. */
void qsched_init ( struct qsched *s , int nr_queues , int flags );
qsched_res_t qsched_addres ( struct qsched *s , qsched_res_t parent );
void qsched_addlock ( struct qsched *s , qsched_task_t t , qsched_res_t res );
void qsched_addunlock ( struct qsched *s , qsched_task_t ta , qsched_task_t tb );
qsched_task_t qsched_addtask ( struct qsched *s , int type , unsigned int flags , void *data , int data_size , int cost );
void qsched_adduse ( struct qsched *s , qsched_task_t t , qsched_res_t res );
void qsched_free ( struct qsched *s );
void qsched_run ( struct qsched *s , int nr_threads , qsched_funtype fun );
void qsched_reset ( struct qsched *s );

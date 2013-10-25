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
#define sched_flag_none                 0
#define sched_flag_dirty                1
#define sched_flag_ready                2

/* Some sched-specific constants. */
#define sched_stretch                   1.41
#define sched_inc                       10
#define sched_init_depspertask          2
#define sched_init_lockspertask         2
#define sched_init_usespertask          2
#define sched_init_respertask           2
#define sched_init_datapertask          2
#define sched_data_round                16
#define sched_res_none                  (-1)


/* The sched data structre. */
struct sched {

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

    };


/* Function prototypes. */
void sched_init ( struct sched *s , int nr_queues , int size );
void sched_sort ( int *restrict data , int *restrict ind , int N , int min , int max );
void sched_sort_rec ( int *restrict data , int *restrict ind , int N , int min , int max );
void sched_prepare ( struct sched *s );
int sched_addres ( struct sched *s , int parent );
void sched_addlock ( struct sched *s , int t , int res );
void sched_addunlock ( struct sched *s , int ta , int tb );
int sched_newtask ( struct sched *s , int type , int subtype , unsigned int flags , void *data , int data_size , int cost );
struct task *sched_gettask ( struct sched *s , int qid );
void sched_adduse ( struct sched *s , int t , int res );
void sched_done ( struct sched *s , struct task *t );
void *sched_getdata( struct sched *s , struct task *t );
void sched_free ( struct sched *s );
int sched_lockres ( struct sched *s , int rid );
void sched_unlockres ( struct sched *s , int rid );
int sched_locktask ( struct sched *s , int tid );
void sched_unlocktask ( struct sched *s , int tid );


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <R.h>

#define EPSILON 0.01
#define ROUND_CNST 1000

/*-----------------------------------------------------*/
/* ------------ priority queue pqops ------------------*/
#ifndef CS161_PROJLIB
#define CS161_PROJLIB

typedef unsigned long long PriType;

extern PriType newPriority(int CPUusage, int nice);
extern int ageCPUusage(int CPUusage, int nice);
extern int compare(PriType a, PriType b);
extern void initialize(void);
extern void finalize(void);

extern PriType strtopri(char *s);
extern char *pritostr(PriType p);
extern void die(char *s);
extern void warn(char *s);
extern void readCommand(void);
extern void printTop(int jobid, int duration, int nice, int heapsize);

extern char cmdletter;
extern int cmdtick, cmdjobid, cmdduration, cmdnice;
extern int tick;

#define CPU_PER_TICK 100
#define TICKS_UNTIL_AGING 100

#endif

#define FREE(x)  free(x) ; x = NULL /* core if free'ed pointer is used */
#define LEFT(x)  (2*x)          /* left child of a node */
#define RIGHT(x) ((2*x)+1)      /* right child of a node */
#define PARENT(x) (x/2)         /* parent of a node */
#define SWAP(t,x,y) tmp = x ; x = y ; y = tmp   /* swap to variables */

#define MSGSIZE 128             /* max size of a debug message */
#define FAILED   -1             /* define -1 to indicate failure of some
                                   operations */
#define EMPTY     1             /* define 1 to indicate that the heap is empty 
                                 */

/* just in case NULL is not defined */
#ifndef NULL
#define NULL    0
#endif

/* number of observed probes*/
int num_of_probes;

/* number of samples*/
int num_of_samples;

/* constant used in the stop condition*/
float beta_value;

char messages[MSGSIZE];

typedef float priority;

/* define the heap object*/
typedef struct node {
    priority p;
    unsigned int id;
    short sentinel;
} node;

/* define a structure representing an individual node in the heap, and
 * make it a valid type for convenience */
typedef float *sum_ptr;
typedef struct breakpoint {
    priority p;
    unsigned int id;
    int size;
    unsigned int prev;
    unsigned int next;
    sum_ptr sum;
} breakpoint;

/* create a global node tmp, for swaping purposes */
node tmp;

/* for convience in function declarations, typedef a pointer to a node
 * as its own type, node_ptr */
typedef node *node_ptr;


/* define a structure representing the heap, and make it a valid type
 * for convenience */
typedef struct binary_heap {
    int heap_size;
    int max_elems;
    node_ptr elements;
} binary_heap;

/* function prototypes for functions which operate on a binary heap */
extern void heapify(binary_heap * a, int i);
extern node_ptr heap_max(binary_heap * a);
extern node heap_extract_max(binary_heap * a);
extern void heap_insert(binary_heap * a, node key);
extern void heap_delete(binary_heap * a, int i);
extern void heap_increase_key(binary_heap * a, int i, priority p);
extern void heap_initalize(binary_heap * a, int nodes);
extern void heap_finalize(binary_heap * a);

/* function prototypes for functions which operate on a node */
extern int node_find(binary_heap a, unsigned int id);
extern node node_create(unsigned int id, priority p);
extern breakpoint breakpoint_create(unsigned int id, priority p, float sum[],
                                    int size, unsigned int prev,
                                    unsigned int next);

/* function prototypes for helper functions */
extern int compare_priority(node i, node j);
extern void print_error(char *msg);


void heapify(binary_heap * a, int i)
{
    register int l, r, largest;

    l = LEFT(i);
    r = RIGHT(i);

    /*
       check the left child 
     */
    largest = ((l <= a->heap_size &&
                compare_priority(a->elements[l], a->elements[i])) ? l : i);

    /*
       check the right child 
     */
    if (r <= a->heap_size &&
        compare_priority(a->elements[r], a->elements[largest]))
        largest = r;

    if (largest != i) {

        /*
           swap nodes largest and i, then heapify 
         */

        SWAP(node, a->elements[i], a->elements[largest]);
        heapify(a, largest);
    }
}

/* Function to return the max (first) node of a heap */

node_ptr heap_max(binary_heap * a)
{
    return ((a->heap_size <= 0) ? NULL : &(a->elements[1]));
}

/* Function to remove the max node from the heap and return it.  The
 * running time is O(lg(n)) since it performs only a costant amount of
 * work on top of the O(lg(n)) of heapify(). Adapted from Introduction
 * to Algorithms (Cormen, Leiserson, Rivest 1990) page 150 */

node heap_extract_max(binary_heap * a)
{
    node max;

    max.sentinel = 1;

    /*
       if there are elements in the heap, make the last item in the heap the
       first one, shorten the heap by one and call heapify(). 
     */

    if (a->heap_size >= 1) {
        max = a->elements[1];
        a->elements[1] = a->elements[(a->heap_size)--];
        heapify(a, 1);
    }

    return max;
}

/* Function to insert an element into the heap, worst case running
 * time is O(lg(n)) on an n element heap, since the path traced from
 * the new leaf to the root has at most length lg(n). This occurs when
 * the new leaf should be the root node.  Adapted from Introduction to
 * Algorithms (Cormen, Leiserson, Rivest 1990) page 150 */

void heap_insert(binary_heap * a, node key)
{
    register int i;

    /*
       if the heap already has the max number of elements we do not allow more 
       elements to be added 
     */
    if (a->heap_size >= a->max_elems) {
        print_error("Heap capacity exceeded, new element not added.");
        return;
    }

    /*
       increase the heap size to accomidate the new node, and set the inital
       position of this node to be the last node in the heap 
     */
    i = ++(a->heap_size);

    /*
       traverse the parth from the leaf to the root to find the a proper
       place for the new element 
     */
    while (i > 1 && compare_priority(key, a->elements[PARENT(i)])) {
        a->elements[i] = a->elements[PARENT(i)];
        i = PARENT(i);
    }

    /*
       insert the element at the position that was determined 
     */
    a->elements[i] = key;
}

/* Function to delete a node from the heap. Adapted from Introduction
 * to Algorithms (Cormen, Leiserson, Rivest 1990) page 151 Exercise
 * 7.5-5 */

void heap_delete(binary_heap * a, int i)
{
    node deleted;

    /*
       return with an error if the input is invalid, ie trying to delete
       elements that are outside of the heap bounds, 1 to heap_size 
     */
    if (i > a->heap_size || i < 1) {
        sprintf(messages, "heap_delete(): %d, no such element.", i);
        print_error(messages);
        return;
    }

    /*
       switch the item to be deleted with the last item, and then shorten the
       heap by one 
     */
    deleted = a->elements[i];
    a->elements[i] = a->elements[(a->heap_size)--];

    heapify(a, i);
}


/* Function to increase the key value of a node from in the
 * heap. Adapted from Introduction to Algorithms (Cormen, Leiserson,
 * Rivest 1990) page 151 Exercise 7.5-4 */
void heap_increase_key(binary_heap * a, int i, priority p)
{

    /*
       return with an error if the input is invalid, ie trying to increase
       elements that are outside of the heap bounds, 1 to heap_size 
     */
    if (i > a->heap_size || i < 1) {
        sprintf(messages, "heap_increase_key(): %d, no such element.", i);
        print_error(messages);
        return;
    }

    /*
       change and propagate 
     */
    a->elements[i].p = p;
    heapify(a, i);
}

/* function to initalize a given binary heap */
void heap_initalize(binary_heap * a, int nodes)
{

    /*
       We initalize heap_size to zero, since a newly created heap contains no
       elements. 
     */
    a->heap_size = 0;

    /*
       we set the max elems to the requested number of nodes, and the allocate 
       enough space for this + 1 number of nodes, since the heap is always
       numbered from 1, but array/pointer accesses are always from 0. 
     */
    a->max_elems = nodes;
    a->elements = (node_ptr) malloc(sizeof (node) * ((a->max_elems) + 1));

    /*
       mark the zero'th element of the heap a to be empty, just in case it is
       every accessed 
     */
    a->elements[0].sentinel = 1;
}


/* function to clean up after we are done with the heap */
void heap_finalize(binary_heap * a)
{
    FREE(a->elements);
}

/* function to create a node */
node node_create(unsigned int id, priority p)
{
    node n;
    n.id = id;
    n.p = p;
    n.sentinel = 0;
    return n;
}

/* function to create a breakpoint */
breakpoint breakpoint_create(unsigned int id, priority p, float sum[],
                             int size, unsigned int prev, unsigned int next)
{
    breakpoint b;
    b.id = id;
    b.p = p;
    b.size = size;
    b.prev = prev;
    b.next = next;
    b.sum = (sum_ptr) malloc(sizeof (float) * (num_of_samples));
    int i;
    for (i = 0; i < num_of_samples; i++) {
        b.sum[i] = (float) sum[i];

    }
    return b;
}

/* function to compare the priority of two given nodes, this is a
 * wrapper for the given compare routine, since in all heap
 * comparisions, we are only interested in greater than or less than
 * operations */
int compare_priority(node i, node j)
{
    if (i.p > j.p) {
        return 1;
    } else {
        if (i.p < j.p) {
            return 0;
        } else {
            if (i.id < j.id)
                return 1;
        }
    }
    return 0;
}

/* function to find if a node is in the heap, O(n) worst case, since
 * we will have to consider every element in a failed search */
int node_find(binary_heap a, unsigned int id)
{
    register int i;

    for (i = 1; i <= a.heap_size; i++)
        if (id == a.elements[i].id)
            return i;
    return FAILED;
}

/* function to print an error message */
void print_error(char *msg)
{
    Rprintf("# ERROR: %s\n", msg);
}

/* ------------ end priority queue pqops ------------------*/
/*---------------------------------------------------------*/

/* function prototype to compute the index entry of the matrix*/
int compute_index(int row, int col, int num_of_col);


int compute_index(int row, int col, int num_of_col)
{
    int index = (row * num_of_col) + col;
    return (index);
}

/* function prototypes for segmentation */

/* main function */
void vegaMC(float **data, int *markers_start, int *start, int *end, int *size,
            float *mean, int *n, float *be, float *std, int *n_reg,
            int *n_samples, float *weight, float *weight_sum);

/* this function computes the initial trivial segmentation*/
/* DA MODIFICARE*/
void init_trivial_segmentation(breakpoint * brks, float **data,
                               binary_heap * b, float weight[]);

/* this function updates the priority for the i-th breakpoint*/
float update_priority(breakpoint * brks, int i, float weight[]);


/*
data= the data matrix NxM (N samples and M markers)
markers_start = start/end marker
start/end/size/mean = output matrix pointers
n = number of probes
be = beta_value
std = data imput standard deviation
n_reg = number of comuted regions
n_samples = number of observed samples
weight = weight associated to each sample
weight_sum = the sum of the weight
*/


void vegaMC(float **data, int *markers_start, int *start, int *end, int *size,
            float *mean, int *n, float *be, float *std, int *n_reg,
            int *n_samples, float *weight, float *weight_sum)
{
    int i, j, k;
    float lambda, lambda_gradient, stop_lambda_gradient, tmp_lambda,
        tmp_sd = 0;
    num_of_probes = *n;
    num_of_samples = *n_samples;
    beta_value = *be;
    breakpoint brks[num_of_probes + 1], brk_del;
    binary_heap b;
    int node_index;
    node max, tmp_node;

    /*
       compute the stop condition based on the standard deviation of each
       sample, the weigth associated to each sample and the constant beta_value 
     */
    for (j = 0; j < num_of_samples; j++) {
        tmp_sd += *(std + j) * (*(weight + j));
    }

    
    stop_lambda_gradient = (tmp_sd / (*weight_sum)) * beta_value;
    
    /*
       Init"ialize the binary heap 
     */
    heap_initalize(&b, num_of_probes - 1);

    /*
       Compute the trivial segmenatation 
     */
    init_trivial_segmentation(brks, data, &b, weight);
    lambda = (heap_max(&b)->p) - EPSILON;
    lambda_gradient = fabs(lambda);


    while ((b.heap_size > 0) && (lambda_gradient <= stop_lambda_gradient)) {
        //Rprintf("b.heap_size %d \n", b.heap_size); || b.heap_size > 2
        
        while ((b.heap_size > 0) && (heap_max(&b)->p > lambda)) {

            max = heap_extract_max(&b);

            if ((max.p > (brks + max.id)->p) && ((brks + max.id)->p < lambda)) {
                max.p = (brks + max.id)->p;
                heap_insert(&b, max);

            } else {

                brk_del = *(brks + max.id);

                (brks + brk_del.prev)->next = brk_del.next;
                (brks + brk_del.next)->prev = brk_del.prev;

                for (i = 0; i < num_of_samples; i++) {
                    (brks + brk_del.prev)->sum[i] =
                        (brks + brk_del.prev)->sum[i] + brk_del.sum[i];

                }

                (brks + brk_del.prev)->size =
                    (brks + brk_del.prev)->size + brk_del.size;

                tmp_lambda = update_priority(brks, brk_del.prev, weight);

                if (brk_del.prev > 0) {
                    node_index = node_find(b, brk_del.prev);
                    tmp_node = b.elements[node_index];
                    heap_delete(&b, node_index);
                    tmp_node.p = tmp_lambda;
                    heap_insert(&b, tmp_node);
                }
                (brks + brk_del.prev)->p = tmp_lambda;




                if (brk_del.next < num_of_probes) {
                    tmp_lambda = update_priority(brks, brk_del.next, weight);
                    node_index = node_find(b, brk_del.next);
                    tmp_node = b.elements[node_index];
                    heap_delete(&b, node_index);
                    tmp_node.p = tmp_lambda;

                    heap_insert(&b, tmp_node);
                }
                (brks + brk_del.next)->p = tmp_lambda;
            }
        }


        /*
           Update the current lambda value 
         */
        if ((b.heap_size > 0)) {
            lambda_gradient = (float) lambda - (heap_max(&b)->p);
            lambda = (heap_max(&b)->p) - EPSILON;
        }
    }


    for (i = 0, j = 0; i < num_of_probes; j++) {
        brk_del = *(brks + i);

        start[j] = markers_start[brk_del.id];
        end[j] = markers_start[brk_del.next - 1];
        size[j] = brk_del.size;
        tmp_sd = 0;
        for (k = 0; k < num_of_samples; k++) {
            tmp_sd +=
                (((float) brk_del.sum[k]) / ((float) brk_del.size)) *
                (float) *(weight + k);
        }
        mean[j] = tmp_sd / (*weight_sum);
        if (i == 1) {
        }
        i = brk_del.next;

    }
    *n_reg = j;
    heap_finalize(&b);
}






float update_priority(breakpoint * brks, int i, float weight[])
{
    if (i != 0) {
        float lambda, first_term, second_term, prev_mean, curr_mean;
        int prev_size, curr_size, prev_id;
        float tmp;
        breakpoint brk_curr, brk_prev;
        prev_id = (brks + i)->prev;

        brk_curr = *(brks + i);
        brk_prev = *(brks + prev_id);
        prev_size = brk_prev.size;
        curr_size = brk_curr.size;

        second_term = 0;

        for (i = 0; i < num_of_samples; i++) {
            prev_mean = (float) (brk_prev.sum[i]) / (prev_size);
            curr_mean = (float) (brk_curr.sum[i]) / (curr_size);
            tmp = (float) prev_mean - curr_mean;

            if (tmp < 0)
                tmp = -1 * tmp;

            second_term += (float) pow(tmp, 2) * weight[i];

        }

        first_term = (float) (prev_size * curr_size) / (prev_size + curr_size);
        lambda = first_term * second_term;
        return (-lambda);
    } else {
        return (float) -1;
    }
}


void init_trivial_segmentation(breakpoint * brks, float **data,
                               binary_heap * b, float weight[])
{

    int i, j, index;
    float tmp_lambda, data_prev, data_curr, tmp_sum[num_of_samples];
    float a, c;


    /*
       The breakpoint 0 has no prev 
     */
    for (i = 0; i < num_of_samples; i++) {
        index = compute_index(i, 0, num_of_probes);
        tmp_sum[i] = floor(data[i][0] * ROUND_CNST) / ROUND_CNST;
    }

    *(brks) = breakpoint_create(0, -1, tmp_sum, 1, -1, 1);

    for (i = 0; i < num_of_samples; i++) {
        tmp_sum[i] = -1;
    }
    /*
       The breakpoint N has no next/prev and no data it is just a closure 
     */
    *(brks + num_of_probes) =
        breakpoint_create(num_of_probes, -1, tmp_sum, -1, -1, -1);

    /*
       Initialize all other nodes in the trivial segmentation 
     */
    for (i = 1; i < num_of_probes; i++) {
        /*
           Compute the lambda for the breakpoints as 1/2*(sum_1:NSamples
           lambda) 
         */
        tmp_lambda = 0;
        for (j = 0; j < num_of_samples; j++) {
            index = compute_index(j, i - 1, num_of_probes);
            data_prev = data[j][i - 1];

            index = compute_index(j, i, num_of_probes);
            data_curr = data[j][i];

            tmp_sum[j] = data_curr;
            a = data_prev - data_curr;
            if (a < 0) {
                a = -1 * a;
            }
            c = pow(a, 2);
            tmp_lambda += c * weight[j];
        }
        tmp_lambda = 0.5 * tmp_lambda;
        *(brks + i) =
            breakpoint_create(i, -tmp_lambda, tmp_sum, 1, i - 1, i + 1);
        heap_insert(b, node_create(i, -tmp_lambda));
    }
}

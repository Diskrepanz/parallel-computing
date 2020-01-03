/*
 * CSE 5441 : Lab 2
 * Sheng Ding
 * ding.853@osu.edu
 * The Ohio State University
 */

#include "Sheng_Ding_omp_persistent.h"

/*
 * Finds the contact distance between two boxes.
 * grid - The input grid.
 * current_box_id	- The ID of the current box.
 * neighbor_box_id	- The ID of the neighbor box.
 * direction		- The direction of the neighbor box in respect with the
 *					  current box. Either vertical or horizontal.
 *
 * Returns the contact distance between the two boxes.
 */
int find_contact_distance(int current_box_id, int neighbor_box_id, int direction){
    /* Box properties and the contact distance. */
    int contact_distance = 0;
    int current_box_left;
    int current_box_right;
    int current_box_top;
    int current_box_bottom;
    int neighbor_box_left;
    int neighbor_box_right;
    int neighbor_box_top;
    int neighbor_box_bottom;

    /* Neighbors can be in horizontal and vertical direction. */
    switch (direction) {
        case HORIZONTAL:
            current_box_left = grid[current_box_id].upper_left_x;
            current_box_right = grid[current_box_id].upper_left_x + grid[current_box_id].width;

            neighbor_box_left = grid[neighbor_box_id].upper_left_x;
            neighbor_box_right = grid[neighbor_box_id].upper_left_x + grid[neighbor_box_id].width;

            contact_distance = min(current_box_right, neighbor_box_right) - max(current_box_left, neighbor_box_left);
            break;
        case VERTICAL:
            current_box_top = grid[current_box_id].upper_left_y;
            current_box_bottom = grid[current_box_id].upper_left_y + grid[current_box_id].height;

            neighbor_box_top = grid[neighbor_box_id].upper_left_y;
            neighbor_box_bottom = grid[neighbor_box_id].upper_left_y + grid[neighbor_box_id].height;

            contact_distance = min(current_box_bottom, neighbor_box_bottom) - max(current_box_top, neighbor_box_top);
            break;
        default:
            /* Error! */
            assert(0);
            break;
    }

    return (contact_distance);
}

/*
 * Computes the contact distances among each neighbor and the current box.
 * grid - The input grid.
 * number_of_boxes  - The total number of boxes in the grid.
 *
 * Returns SUCCESS in case of success.
 */
int compute_contact_distance(){
    /* Iterator to go through the neighbors of a box. */
    neighbor *iter;

    // Compute for all the grid.
    for (int current_box_id = 0; current_box_id < number_of_boxes; current_box_id++) {
        // Compute the contact distances of the top neighbours.
        for (iter = grid[current_box_id].top; iter < grid[current_box_id].top + neighbors[current_box_id][0]; iter++) {
            iter->contact_distance = find_contact_distance(current_box_id, iter->id, HORIZONTAL);
        }
        // Compute the contact distances of the bottom neighbours.
        for (iter = grid[current_box_id].bottom; iter < grid[current_box_id].bottom + neighbors[current_box_id][1]; iter++) {
            iter->contact_distance = find_contact_distance(current_box_id, iter->id, HORIZONTAL);
        }
        // Compute the contact distances of the left neighbors.
        for (iter = grid[current_box_id].left; iter < grid[current_box_id].left + neighbors[current_box_id][2]; iter++) {
            iter->contact_distance = find_contact_distance(current_box_id, iter->id, VERTICAL);
        }
        // Compute the contact distances of the right neighbors.
        for (iter = grid[current_box_id].right; iter < grid[current_box_id].right + neighbors[current_box_id][3]; iter++) {
            iter->contact_distance = find_contact_distance(current_box_id, iter->id, VERTICAL);
        }
    }

    return (SUCCESS);
}

/*
 * Commits the updated DSVs in the original boxes in the grid.
 * Also finds the maximum and minimum box DSV (temp).
 *
 * grid 			- The input grid.
 * updated_dsvs 	- Array that contains the updated DSVs.
 * number_of_boxes  - The total number of boxes in the grid.
 * max 				- The maximum DSV (temperature) on convergence.
 * min 				- The minimum DSV (temperature) on convergence.
 *
 * Returns SUCCESS in case of successful grid update.
 */
int commit_dsv_update() {
    max_dsv = min_dsv = updated_DSVs[0];

    for(int i = 0; i < number_of_boxes; i++) {
        /* Commit the updated DSVs*/
        grid[i].box_dsv = updated_DSVs[i];

        /* Find the max and min DSVs */
        if (grid[i].box_dsv > max_dsv) {
            max_dsv = grid[i].box_dsv;
        } else if(grid[i].box_dsv < min_dsv) {
            min_dsv = grid[i].box_dsv;
        }
    }
    return (SUCCESS);
}

/*
 * Checks whether the DSVs reached convergence or not.
 *
 * grid 			- The input grid.
 * number_of_boxes  - The total number of boxes in the grid.
 * epsilon 			- The given epsilon.
 * max 				- The maximum DSV (temperature) on convergence.
 * min 				- The minimum DSV (temperature) on convergence.
 *
 * Returns SUCCESS in case of convergence reached, otherwise ERROR.
 */
int check_for_convergence() {
    if ((max_dsv - min_dsv) <= (max_dsv * epsilon)) {
        return (SUCCESS);
    }
    return (ERROR);
}

/*
 * All the new threads call this function and iterate the computation
 * until convergence reached.  After all threads computes the updated
 * DSVs, thread number 0 commits it and updates the maximum and minimum
 * DSVs.
 *
 * message - Thread specific message that provides information on the
 * 	     grid partitioning..
 */
void* converge_DSV(void *message) {

    /* Iterator to go through the neighbors of a box. */
    neighbor *iter;

    /* Used as a temporary variable to store the
    total DSV to aid in computing the weighed average. */
    float total_DSV = 0;

    /* Average adjacent temperature for each current box. */
    float average_DSV;
    int result;
    thread_data *mydata = (thread_data*) message;

    do {
        /*
         * Go through allowed boxes in the grid and compute updated DSV.
         * convergence loop.
         */
        for (int i = mydata->start; i <= mydata->end; i++) {
            total_DSV = 0;

            /* Check the top neighbors. */
            if (neighbors[i][0] == 0) {
                total_DSV += grid[i].box_dsv * (float)grid[i].width;
            } else {
                for (iter = grid[i].top; iter < grid[i].top + neighbors[i][0]; iter++) {
                    total_DSV += grid[iter->id].box_dsv * (float)(iter->contact_distance);
                }
            }

            /* Check the bottom neighbors. */
            if (neighbors[i][1] == 0) {
                total_DSV += grid[i].box_dsv * (float)grid[i].width;
            } else {
                for (iter = grid[i].bottom; iter != grid[i].bottom + neighbors[i][1]; iter++) {
                    total_DSV += grid[iter->id].box_dsv * (float)(iter->contact_distance);
                }
            }

            /* Check the left neighbors. */
            if (neighbors[i][2] == 0){
                total_DSV += grid[i].box_dsv * (float)grid[i].height;
            } else {
                for (iter = grid[i].left; iter != grid[i].left + neighbors[i][2]; iter++) {
                    total_DSV += grid[iter->id].box_dsv * (float)(iter->contact_distance);
                }
            }

            /* Check the right neighbors. */
            if (neighbors[i][3] == 0) {
                total_DSV += grid[i].box_dsv * (float)grid[i].height;
            } else {
                for (iter = grid[i].right; iter != grid[i].right + neighbors[i][3]; iter++) {
                    total_DSV += grid[iter->id].box_dsv * (float)(iter->contact_distance);
                }
            }

            /* Compute the average DSV. */
            average_DSV = (total_DSV / (float)grid[i].perimeter);

            /* Compute the updated values based on the average adjacent DSV. */
            if (average_DSV >= grid[i].box_dsv) {
                updated_DSVs[i] = grid[i].box_dsv +
                                  ((average_DSV - grid[i].box_dsv) * affect_rate);
            } else if (average_DSV < grid[i].box_dsv) {
                updated_DSVs[i] = grid[i].box_dsv -
                                  ((grid[i].box_dsv - average_DSV) * affect_rate);
            }
        }

        // Wait till all the threads reach here.
        #pragma omp barrier

        // Let thread 0 commit the result, update max and
        // min, and update iteration number.
        if (mydata->thread_id == 0) {
            result = commit_dsv_update();
            assert(result == SUCCESS);
            convergence_iteration = convergence_iteration + 1;
        }
        // Let thread 0 completes its update.
        #pragma omp barrier

        result = check_for_convergence();
    } while (result != SUCCESS);
}

/*
 * Function that does the stencil computation.
 * The iterative computations are repeated until convergence.
 *
 * Returns SUCCESS in case of successful computation.
 */
int do_stencil_computation() {

    int result;

    convergence_iteration = 0;
    thread_data message[number_of_threads]; //// Message to be passed to each thread.
    int start_index = 0;
    int number_of_partition = number_of_boxes/number_of_threads;

    // Partition the grid for each Thread.
    for (int i = 0; i < number_of_threads; i++) {
        message[i].thread_id = i;
        message[i].start = start_index;
        if (i == (number_of_threads - 1)) {
            message[i].end = number_of_boxes - 1;
        } else {
            message[i].end = start_index + number_of_partition - 1;
        }
        start_index = start_index + number_of_partition;
    }

    // create a barrier object with a count of number of threads
    #pragma omp parallel num_threads(number_of_threads)
    {
        //Get the number of thread and send its data.
        int thread_id = omp_get_thread_num();
        //See the actual number of threads created.
        if (thread_id == 0) {
            // Let the master thread update the
            // actual number of threads created.
            actual_number_of_threads = omp_get_num_threads();
        }
        converge_DSV((void *) &message[thread_id]);
    } // implicit barrier point.

    return (SUCCESS);
}

int main(int argc, char* argv[]){

    int result;

    /* The timers. */
    time_t time1, time2;
    clock_t clock1, clock2;
    struct timespec start = {0,0}, end = {0,0};

    /* Check whether wrong number of arguments. */
    //assert(argc == 4);

    /* Read affect_rate and epsilon. */
    affect_rate = atof(argv[1]);
    epsilon = atof(argv[2]);
    number_of_threads = atof(argv[3]);

    /* Read the number of grid and create grid. */
    scanf("%d",&number_of_boxes);

    grid = (box*)malloc(sizeof(box)*number_of_boxes);
    updated_DSVs = (float*)malloc(sizeof(float)*number_of_boxes);
    neighbors = (int**)malloc(sizeof(int*)*number_of_boxes);
    for(int i = 0;i< number_of_boxes;i++){
        neighbors[i] = (int*)malloc(sizeof(int)*4);
    }

    /*
     * Reads the input and prepare the data structures.
     *
     * grid 			- The input grid.
     * number_of_boxes  - The total number of boxes in the grid.
     *
     *  Returns SUCCESS in case of successful reading otherwise ERROR.
    */
    int current_box_id;
    int num_grid_rows;
    int num_grid_cols;
    int last_line;

    scanf("%d",&num_grid_rows);
    scanf("%d",&num_grid_cols);

    for (int i = 0; i < number_of_boxes; i++) {

        scanf("%d",&current_box_id);
        scanf("%d",&grid[current_box_id].upper_left_y);
        scanf("%d",&grid[current_box_id].upper_left_x);
        scanf("%d",&grid[current_box_id].height);
        scanf("%d",&grid[current_box_id].width);

        grid[current_box_id].perimeter = 2 * (grid[current_box_id].width + grid[current_box_id].height);

        //define neighbor vector
        grid[current_box_id].top = (neighbor*)malloc(sizeof(neighbor)*number_of_boxes);
        grid[current_box_id].bottom = (neighbor*)malloc(sizeof(neighbor)*number_of_boxes);
        grid[current_box_id].left = (neighbor*)malloc(sizeof(neighbor)*number_of_boxes);
        grid[current_box_id].right = (neighbor*)malloc(sizeof(neighbor)*number_of_boxes);

        //read neighbors
        neighbor temp;
        scanf("%d",&neighbors[current_box_id][0]);
        for (int j = 0; j < neighbors[current_box_id][0]; j++) {
            scanf("%d",&temp.id);
            *(grid[current_box_id].top + j) = temp;
        }
        scanf("%d",&neighbors[current_box_id][1]);
        for (int j = 0; j < neighbors[current_box_id][1]; j++) {
            scanf("%d",&temp.id);
            *(grid[current_box_id].bottom + j) = temp;
        }
        scanf("%d",&neighbors[current_box_id][2]);
        for (int j = 0; j < neighbors[current_box_id][2]; j++) {
            scanf("%d",&temp.id);
            *(grid[current_box_id].left + j) = temp;
        }
        scanf("%d",&neighbors[current_box_id][3]);
        for (int j = 0; j < neighbors[current_box_id][3]; j++) {
            scanf("%d",&temp.id);
            *(grid[current_box_id].right + j) = temp;
        }

        //read the box_dsv
        scanf("%f", &grid[current_box_id].box_dsv);
    }

    scanf("%d",&last_line);

    if (last_line != -1) {
        return (ERROR);
    }

    result =  compute_contact_distance();
    assert(result == SUCCESS);

    printf("\n*******************************************************************\n");

    /*Take the time stamps. */
    time(&time1);
    clock1 = clock();
    clock_gettime(CLOCK_REALTIME,&start);


    result = do_stencil_computation();
    assert(result == SUCCESS);

    /* Take the time stamps again. */
    time(&time2);
    clock2 = clock();
    clock_gettime(CLOCK_REALTIME,&end);
    float time_taken = ((float)(end.tv_sec - start.tv_sec)* 1000+ (float)(end.tv_nsec - start.tv_nsec)/1000000);

    printf("Threads requested: %d; Threads created: %d.\n", number_of_threads, actual_number_of_threads);
    printf("Dissipation converged in %d iterations.\n",convergence_iteration);
    printf("With max DSV = %f and min DSV = %f.\n", max_dsv, min_dsv);
    printf("Affect rate = %f; Epsilon: %f.\n",affect_rate, epsilon);
    printf("Elapsed convergence loop time (clock) : %d.\n", (int)(clock2-clock1));
    printf("Elapsed convergence loop time (time) : %d.\n", (int)difftime(time2, time1));
    printf("Elapsed convergence loop time (chrono) : %2f.\n", time_taken);
    printf("*******************************************************************\n");

    return 0;
}
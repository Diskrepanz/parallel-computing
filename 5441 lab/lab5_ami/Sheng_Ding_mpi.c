/*
 * CSE 5441 : Lab 5
 * Sheng Ding
 * ding.853@osu.edu
 * The Ohio State University
 */

#include "Sheng_Ding_mpi.h"

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
    * Reads the input and prepare the data structures.
    *
    * grid 			- The input grid.
    * number_of_boxes  - The total number of boxes in the grid.
    *
    *  Returns SUCCESS in case of successful reading otherwise ERROR.
   */
int read_input_file() {

    FILE  *fp = NULL;
    fp = fopen("testgrid_400_12206","r");

    /* Read the number of grid and create grid. */
    fscanf(fp,"%d",&number_of_boxes);

    grid = (box*)malloc(sizeof(box)*number_of_boxes);
    updated_DSVs = (float*)malloc(sizeof(float)*number_of_boxes);
    neighbors = (int**)malloc(sizeof(int*)*number_of_boxes);
    for(int i = 0;i< number_of_boxes;i++){
        neighbors[i] = (int*)malloc(sizeof(int)*4);
    }

    int current_box_id;
    int num_grid_rows;
    int num_grid_cols;
    int last_line;
    int result;

    fscanf(fp,"%d", &num_grid_rows);
    fscanf(fp,"%d", &num_grid_cols);

    for (int i = 0; i < number_of_boxes; i++) {

        fscanf(fp,"%d", &current_box_id);
        fscanf(fp,"%d", &grid[current_box_id].upper_left_y);
        fscanf(fp,"%d", &grid[current_box_id].upper_left_x);
        fscanf(fp,"%d", &grid[current_box_id].height);
        fscanf(fp,"%d", &grid[current_box_id].width);

        grid[current_box_id].perimeter = 2 * (grid[current_box_id].width + grid[current_box_id].height);

        //define neighbor vector
        grid[current_box_id].top = (neighbor *) malloc(sizeof(neighbor) * number_of_boxes);
        grid[current_box_id].bottom = (neighbor *) malloc(sizeof(neighbor) * number_of_boxes);
        grid[current_box_id].left = (neighbor *) malloc(sizeof(neighbor) * number_of_boxes);
        grid[current_box_id].right = (neighbor *) malloc(sizeof(neighbor) * number_of_boxes);

        //read neighbors
        neighbor temp;
        fscanf(fp,"%d", &neighbors[current_box_id][0]);
        for (int j = 0; j < neighbors[current_box_id][0]; j++) {
            fscanf(fp,"%d", &temp.id);
            *(grid[current_box_id].top + j) = temp;
        }
        fscanf(fp,"%d", &neighbors[current_box_id][1]);
        for (int j = 0; j < neighbors[current_box_id][1]; j++) {
            fscanf(fp,"%d", &temp.id);
            *(grid[current_box_id].bottom + j) = temp;
        }
        fscanf(fp,"%d", &neighbors[current_box_id][2]);
        for (int j = 0; j < neighbors[current_box_id][2]; j++) {
            fscanf(fp,"%d", &temp.id);
            *(grid[current_box_id].left + j) = temp;
        }
        fscanf(fp,"%d", &neighbors[current_box_id][3]);
        for (int j = 0; j < neighbors[current_box_id][3]; j++) {
            fscanf(fp,"%d", &temp.id);
            *(grid[current_box_id].right + j) = temp;
        }

        //read the box_dsv
        fscanf(fp,"%f", &grid[current_box_id].box_dsv);
    }

    fscanf(fp,"%d", &last_line);

    if (last_line != -1) {
        return (ERROR);
    }

    result = compute_contact_distance();
    assert(result == SUCCESS);

    return(SUCCESS);
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
    int iter;

    /* Used as a temporary variable to store the
    total DSV to aid in computing the weighed average. */
    float total_DSV = 0;

    /* Average adjacent temperature for each current box. */
    float average_DSV;
    int result;
    thread_data *mydata = (thread_data*) message;

        /*
         * Go through allowed boxes in the grid and compute updated DSV.
         * convergence loop.
         */
        for (int i = mydata->start; i <= mydata->end; i++) {
            total_DSV = 0;

            /* Check the top neighbors. */
            if (neighborstop[i] == 0) {
                total_DSV += grid_box_dsv[i] * (float)grid_width[i];
            } else {
                for (iter = 0; iter < neighborstop[i]; iter++) {
                    total_DSV += grid_box_dsv[grid_top_id[i * number_of_boxes + iter]] * (float)(grid_top_contact_distance[i * number_of_boxes + iter]);
                }
            }

            /* Check the bottom neighbors. */
            if (neighborsbottom[i] == 0) {
                total_DSV += grid_box_dsv[i] * (float)grid_width[i];
            } else {
                for (iter = 0; iter < neighborsbottom[i]; iter++) {
                    total_DSV += grid_box_dsv[grid_bottom_id[i * number_of_boxes + iter]] * (float)(grid_bottom_contact_distance[i * number_of_boxes + iter]);
                }
            }

            /* Check the left neighbors. */
            if (neighborsleft[i] == 0){
                total_DSV += grid_box_dsv[i] * (float)grid_height[i];
            } else{
                for (iter = 0; iter < neighborsleft[i]; iter++) {
                    total_DSV += grid_box_dsv[grid_left_id[i * number_of_boxes + iter]] * (float)(grid_left_contact_distance[i * number_of_boxes + iter]);
                }
            }

            /* Check the right neighbors. */
            if (neighborsright[i] == 0) {
                total_DSV += grid_box_dsv[i] * (float)grid_height[i];
            } else {
                for (iter = 0; iter < neighborsright[i]; iter++) {
                    total_DSV += grid_box_dsv[grid_right_id[i * number_of_boxes + iter]] * (float)(grid_right_contact_distance[i * number_of_boxes + iter]);
                }
            }

            /* Compute the average DSV. */
            average_DSV = (total_DSV / (float)grid_perimeter[i]);

            /* Compute the updated values based on the average adjacent DSV. */
            if (average_DSV >= grid_box_dsv[i]) {
                updated_DSVs[i] = grid_box_dsv[i] + (average_DSV - grid_box_dsv[i]) * affect_rate;
            } else if (average_DSV < grid_box_dsv[i]) {
                updated_DSVs[i] = grid_box_dsv[i] - (grid_box_dsv[i] - average_DSV) * affect_rate;
            }
        }

        // Wait till all the threads reach here.
        #pragma omp barrier
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
    thread_data message[size]; // Message to be passed to each thread.
    int start_index = (rank-1) * length_of_rank[1];
    int number_of_partition = length_of_rank[rank]/size;

    // Partition the grid for each Thread.
    for (int i = 0; i < size; i++) {
        message[i].thread_id = i;
        message[i].start = start_index;
        if (i == (size - 1)) {
            if(rank == 4){
                message[i].end = number_of_boxes - 1;
            }else {
                message[i].end = rank * length_of_rank[1] - 1;
            }
        } else {
            message[i].end = start_index + number_of_partition - 1;
        }
        start_index = start_index + number_of_partition;
    }

    // create a barrier object with a count of number of threads
    #pragma omp parallel num_threads(size)
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

int main(int argc, char* argv[]) {

    int result;
    int done;

    /* The timers. */
    time_t time1, time2;
    clock_t clock1, clock2;
    struct timespec start = {0, 0}, end = {0, 0};

    /* Read affect_rate and epsilon. */
    affect_rate = atof(argv[1]);
    epsilon = atof(argv[2]);

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        /* Read the input file. */
        result = read_input_file();
        assert(result == SUCCESS);

        grid_left_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_left_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_right_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_right_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_top_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_top_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_bottom_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_bottom_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);

        grid_box_dsv = (float *) malloc(sizeof(float) * number_of_boxes);
        grid_width = (int *) malloc(sizeof(int) * number_of_boxes);
        grid_height = (int *) malloc(sizeof(int) * number_of_boxes);
        grid_perimeter = (int *) malloc(sizeof(int) * number_of_boxes);

        neighborstop = (int *) malloc(sizeof(int) * number_of_boxes);
        neighborsbottom = (int *) malloc(sizeof(int) * number_of_boxes);
        neighborsleft = (int *) malloc(sizeof(int) * number_of_boxes);
        neighborsright = (int *) malloc(sizeof(int) * number_of_boxes);

        for (int i = 0; i < number_of_boxes; i++) {
            for (int j = 0; j < neighbors[i][0]; j++) {
                grid_top_id[i * number_of_boxes + j] = (*(grid[i].top + j)).id;
                grid_top_contact_distance[i * number_of_boxes + j] = (*(grid[i].top + j)).contact_distance;
            }

            for (int j = 0; j < neighbors[i][1]; j++) {
                grid_bottom_id[i * number_of_boxes + j] = (*(grid[i].bottom + j)).id;
                grid_bottom_contact_distance[i * number_of_boxes + j] = (*(grid[i].bottom + j)).contact_distance;
            }

            for (int j = 0; j < neighbors[i][2]; j++) {
                grid_left_id[i * number_of_boxes + j] = (*(grid[i].left + j)).id;
                grid_left_contact_distance[i * number_of_boxes + j] = (*(grid[i].left + j)).contact_distance;
            }

            for (int j = 0; j < neighbors[i][3]; j++) {
                grid_right_id[i * number_of_boxes + j] = (*(grid[i].right + j)).id;
                grid_right_contact_distance[i * number_of_boxes + j] = (*(grid[i].right + j)).contact_distance;
            }

            grid_box_dsv[i] = grid[i].box_dsv;
            grid_width[i] = grid[i].width;
            grid_height[i] = grid[i].height;
            grid_perimeter[i] = grid[i].perimeter;

            neighborstop[i] = neighbors[i][0];
            neighborsbottom[i] = neighbors[i][1];
            neighborsleft[i] = neighbors[i][2];
            neighborsright[i] = neighbors[i][3];
        }

        int set = number_of_boxes/4;

        for (int i = 0; i < 5; i++) {
            if(i == 0){
                length_of_rank[i] = 0;
            }
            if(i != 4){
                length_of_rank[i] = set;
            } 
            else {
                length_of_rank[i] = number_of_boxes - 3 * set;
            }
        }
    }

    MPI_Bcast(&number_of_boxes, 1, MPI_INT, 0, MPI_COMM_WORLD);  //bcast all data needed to all ranks
    MPI_Bcast(length_of_rank, 5, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank != 0) {

        updated_DSVs = (float*)malloc(sizeof(float)*number_of_boxes);

        grid_left_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_left_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_right_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_right_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_top_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_top_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_bottom_id = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);
        grid_bottom_contact_distance = (int *) malloc(sizeof(int) * number_of_boxes * number_of_boxes);

        grid_box_dsv = (float *) malloc(sizeof(float) * number_of_boxes);
        grid_width = (int *) malloc(sizeof(int) * number_of_boxes);
        grid_height = (int *) malloc(sizeof(int) * number_of_boxes);
        grid_perimeter = (int *) malloc(sizeof(int) * number_of_boxes);

        neighborstop = (int *) malloc(sizeof(int) * number_of_boxes);
        neighborsbottom = (int *) malloc(sizeof(int) * number_of_boxes);
        neighborsleft = (int *) malloc(sizeof(int) * number_of_boxes);
        neighborsright = (int *) malloc(sizeof(int) * number_of_boxes);
    }

    MPI_Bcast(grid_top_id, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);  //bcast all data needed to all ranks
    MPI_Bcast(grid_top_contact_distance, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_bottom_id, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_bottom_contact_distance, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_left_id, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_left_contact_distance, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_right_id, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_right_contact_distance, number_of_boxes * number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(grid_box_dsv, number_of_boxes, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_width, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_height, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid_perimeter, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(neighborstop, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(neighborsbottom, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(neighborsleft, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(neighborsright, number_of_boxes, MPI_INT, 0, MPI_COMM_WORLD);

    /*Take the time stamps. */
    time(&time1);
    clock1 = clock();
    clock_gettime(CLOCK_REALTIME, &start);

    do {
        if (rank == 1) {
            done = do_stencil_computation();
            assert(done == SUCCESS);
            MPI_Send(updated_DSVs, length_of_rank[1], MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else if(rank == 2){
            done = do_stencil_computation();
            assert(done == SUCCESS);
            MPI_Send(&updated_DSVs[length_of_rank[1]], length_of_rank[2], MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else if(rank == 3){
            done = do_stencil_computation();
            assert(done == SUCCESS);
            MPI_Send(&updated_DSVs[2 * length_of_rank[1]], length_of_rank[3], MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else if(rank == 4){
            done = do_stencil_computation();
            assert(done == SUCCESS);
            MPI_Send(&updated_DSVs[3 * length_of_rank[1]], length_of_rank[4], MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else {
            MPI_Recv(updated_DSVs, length_of_rank[1], MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&updated_DSVs[length_of_rank[1]], length_of_rank[2], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&updated_DSVs[2 * length_of_rank[1]], length_of_rank[3], MPI_INT, 3, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&updated_DSVs[3 * length_of_rank[1]], length_of_rank[4], MPI_INT, 4, 0, MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank ==0){
            result = commit_dsv_update();
            assert(result == SUCCESS);

            for(int i = 0; i < number_of_boxes; i++) {
                grid_box_dsv[i] = grid[i].box_dsv;//update grid_box_dsv
            }

            convergence_iteration = convergence_iteration + 1;
            result = check_for_convergence();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&convergence_iteration, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(grid_box_dsv, number_of_boxes, MPI_FLOAT, 0, MPI_COMM_WORLD);//new update grid_box_dsv bcast to all ranks
        MPI_Barrier(MPI_COMM_WORLD);

    }while(result != SUCCESS);

    assert(result == SUCCESS);

    /* Take the time stamps again. */
    time(&time2);
    clock2 = clock();
    clock_gettime(CLOCK_REALTIME,&end);
    float time_taken = ((float)(end.tv_sec - start.tv_sec)* 1000+ (float)(end.tv_nsec - start.tv_nsec)/1000000);

    if(rank ==0) {
        printf("\n*******************************************************************\n");
        printf("Dissipation converged in %d iterations.\n", convergence_iteration);
        printf("With max DSV = %f and min DSV = %f.\n", max_dsv, min_dsv);
        printf("Affect rate = %f; Epsilon: %f.\n", affect_rate, epsilon);
        printf("Elapsed convergence loop time (clock) : %d.\n", (int) (clock2 - clock1));
        printf("Elapsed convergence loop time (time) : %d.\n", (int) difftime(time2, time1));
        printf("Elapsed convergence loop time (chrono) : %2f.\n", time_taken);
        printf("*******************************************************************\n");
    }

    MPI_Finalize();
    return 0;
}
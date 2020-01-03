/*
 * CSE 5441 : Lab 4
 * Sheng Ding
 * ding.853@osu.edu
 * The Ohio State University
 */

#include "Sheng_Ding_cuda.h"

/*
 * Finds the contact distance between two boxes.
 * grid - The input grid.
 * current_box_id   - The ID of the current box.
 * neighbor_box_id  - The ID of the neighbor box. 
 * direction        - The direction of the neighbor box in respect with the
 *                    current box. Either vertical or horizontal.
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
 * grid             - The input grid.
 * number_of_boxes  - The total number of boxes in the grid.
 *
 *  Returns SUCCESS in case of successful reading otherwise ERROR.
 */
int read_input_file(){

    int current_box_id;
    int num_grid_rows;
    int num_grid_cols;
    int last_line;
    int result; 

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

    /*
     * Compute the contanct distances of
     * neighbour nodes for later use.
     */
    result =  compute_contact_distance();
    assert(result == SUCCESS);

    return(SUCCESS);
}

/*
 * Commits the updated DSVs in the original boxes in the grid.
 * Also finds the maximum and minimum box DSV (temp).
 *
 * grid             - The input grid.
 * updated_Dsvs     - Array that contains the updated DSVs.
 * number_of_boxes  - The total number of boxes in the grid.
 * max              - The maximum DSV (temperature) on convergence.
 * min              - The minimum DSV (temperature) on convergence.
 *
 * Returns SUCCESS in case of successful grid update.
 */
int commit_dsv_update() {
    max_dsv = min_dsv = updated_DSVs[0];

    for(int i = 0; i < number_of_boxes; i++) {
        /* Commit the updated DSVs*/
        grid_box_dsv[i] = updated_DSVs[i];
        /* Find the max and min DSVs */
        if (grid_box_dsv[i] > max_dsv) {
            max_dsv = grid_box_dsv[i];
        } else if(grid_box_dsv[i] < min_dsv) {
            min_dsv = grid_box_dsv[i];
        }
    }
    return (SUCCESS);
}

/*
 * Checks whether the DSVs reached convergence or not.
 *
 * grid             - The input grid.
 * number_of_boxes  - The total number of boxes in the grid.
 * epsilon          - The given epsilon.
 * max              - The maximum DSV (temperature) on convergence.
 * min              - The minimum DSV (temperature) on convergence.
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
 * Global the function and iterate the computation
 * until convergence reached.  After devices computes the updated
 * DSVs, copy the maximum and minimum DSVs back to Host.
 */
__global__ void compute_updated_DSV(
	int *grid_left_id,int *grid_left_contact_distance,
	int *grid_right_id,int *grid_right_contact_distance,
	int *grid_top_id,int *grid_top_contact_distance,
	int *grid_bottom_id,int *grid_bottom_contact_distance,
	float *grid_box_dsv,int *grid_height,int *grid_width,int *grid_perimeter,
	int *neighborsleft,int *neighborsright,int *neighborstop,int *neighborsbottom,
	int number_of_boxes, float affect_rate,float *update){

    /* Iterator to go through the neighbors of a box. */
    int iter;

    /* Used as a temporary variable to store the
    total DSV to aid in computing the weighed average. */
    float total_DSV;

    /* Average adjacent temperature for each current box. */
    float average_DSV;

    /*
     * Go through allowed boxes in the grid and compute updated DSV.
     */
    int index = blockIdx.x * blockDim.x;

    for(int i = index; i < (blockIdx.x+1) * blockDim.x && i<number_of_boxes; i++) {

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
            update[i] = grid_box_dsv[i] + (average_DSV - grid_box_dsv[i]) * affect_rate;
        } else if (average_DSV < grid_box_dsv[i]) {
            update[i] = grid_box_dsv[i] - (grid_box_dsv[i] - average_DSV) * affect_rate;
        }
    }

}

/*
 * Function that does the stencil computation.
 * The iterative computations are repeated until convergence.
 *
 * Returns SUCCESS in case of successful computation.
 */
int do_stencil_computation(){

    int result = 0;

    int *d_gli;//pointer for device grid
    int *d_glc;
    int *d_gri;
    int *d_grc;
    int *d_gti;
    int *d_gtc;
    int *d_gbi;
    int *d_gbc;

    float *d_gbd;//pointer for device updated dsv
    int *d_gh;//pointer for device box height
    int *d_gw;//pointer for device box width
    int *d_gp;//pointer for device box perimeter

    int *n_l;//pointer for device neighbors
    int *n_r;
    int *n_t;
    int *n_b;

    float *updsv;

    convergence_iteration = 0;

    //把所有需要的数据全部拉出来做多个一维数组!!!很重要!!!相当重要!!!
    int grid_left_id[number_of_boxes * number_of_boxes];
    int grid_left_contact_distance[number_of_boxes * number_of_boxes];
    int grid_right_id[number_of_boxes * number_of_boxes];
    int grid_right_contact_distance[number_of_boxes * number_of_boxes];
    int grid_top_id[number_of_boxes * number_of_boxes];
    int grid_top_contact_distance[number_of_boxes * number_of_boxes];
    int grid_bottom_id[number_of_boxes * number_of_boxes];
    int grid_bottom_contact_distance[number_of_boxes * number_of_boxes];

    grid_box_dsv = (float*)malloc(sizeof(float)*number_of_boxes);
    int grid_width[number_of_boxes];
    int grid_height[number_of_boxes];
    int grid_perimeter[number_of_boxes];

    int neighborstop[number_of_boxes];
    int neighborsbottom[number_of_boxes];
    int neighborsleft[number_of_boxes];
    int neighborsright[number_of_boxes];

    for(int i = 0; i < number_of_boxes; i++){

        for(int j = 0; j < neighbors[i][0]; j++){
            grid_top_id[i * number_of_boxes + j] = (*(grid[i].top + j)).id;
            grid_top_contact_distance[i * number_of_boxes+ j] = (*(grid[i].top + j)).contact_distance;
        }

        for(int j = 0; j < neighbors[i][1]; j++){
            grid_bottom_id[i * number_of_boxes+ j] = (*(grid[i].bottom + j)).id;
            grid_bottom_contact_distance[i * number_of_boxes + j] = (*(grid[i].bottom + j)).contact_distance;
        }

        for(int j = 0; j < neighbors[i][2]; j++){
            grid_left_id[i * number_of_boxes+ j] = (*(grid[i].left + j)).id;
            grid_left_contact_distance[i * number_of_boxes+ j] = (*(grid[i].left + j)).contact_distance;
        }

        for(int j = 0; j < neighbors[i][3]; j++){
            grid_right_id[i * number_of_boxes+ j] = (*(grid[i].right + j)).id;
            grid_right_contact_distance[i * number_of_boxes+ j] = (*(grid[i].right + j)).contact_distance;
        }

        grid_box_dsv[i]=grid[i].box_dsv;
        grid_width[i]=grid[i].width;
        grid_height[i]=grid[i].height;
        grid_perimeter[i]=grid[i].perimeter;

        neighborstop[i] = neighbors[i][0];
        neighborsbottom[i] = neighbors[i][1];
        neighborsleft[i] = neighbors[i][2];
        neighborsright[i] = neighbors[i][3];
    }

    number_of_threads = 200;
    int number_of_blocks = number_of_boxes / number_of_threads;

    // allocate host and device memory
    size_t memSizei;
    size_t memSizef;
    size_t memSizen;

    memSizei = number_of_boxes * sizeof(int);
    memSizef = number_of_boxes * sizeof(float);
    memSizen = number_of_boxes * number_of_boxes * sizeof(int);

    //memory allocate for all variables
    cudaMalloc( (void**) &d_gti, memSizen);
    cudaMalloc( (void**) &d_gtc, memSizen);
    cudaMalloc( (void**) &d_gbi, memSizen);
    cudaMalloc( (void**) &d_gbc, memSizen);
    cudaMalloc( (void**) &d_gli, memSizen);
    cudaMalloc( (void**) &d_glc, memSizen);
    cudaMalloc( (void**) &d_gri, memSizen);        
    cudaMalloc( (void**) &d_grc, memSizen);

    cudaMalloc( (void**) &d_gbd, memSizef);
    cudaMalloc( (void**) &d_gh, memSizei);
    cudaMalloc( (void**) &d_gw, memSizei);
    cudaMalloc( (void**) &d_gp, memSizei);

    cudaMalloc( (void**) &n_t, memSizei);
    cudaMalloc( (void**) &n_b, memSizei);
    cudaMalloc( (void**) &n_l, memSizei);
    cudaMalloc( (void**) &n_r, memSizei);

    cudaMalloc( (void**) &updsv, memSizef);

    //memory copy to device
    cudaMemcpy(d_gti, grid_top_id, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gtc, grid_top_contact_distance, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gbi, grid_bottom_id, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gbc, grid_bottom_contact_distance, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gli, grid_left_id, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_glc, grid_left_contact_distance, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gri, grid_right_id, memSizen, cudaMemcpyHostToDevice);
    cudaMemcpy(d_grc, grid_right_contact_distance, memSizen, cudaMemcpyHostToDevice);

    cudaMemcpy(d_gh, grid_height, memSizei, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gw, grid_width, memSizei, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gp, grid_perimeter,memSizei, cudaMemcpyHostToDevice);

    cudaMemcpy(n_t, neighborstop, memSizei, cudaMemcpyHostToDevice);
    cudaMemcpy(n_b, neighborsbottom, memSizei, cudaMemcpyHostToDevice);
    cudaMemcpy(n_l, neighborsleft, memSizei, cudaMemcpyHostToDevice);
    cudaMemcpy(n_r, neighborsright, memSizei, cudaMemcpyHostToDevice);

    //  launch kernel
    dim3 ThreadsPerBlock(number_of_threads);
    dim3 BlocksPerGrid(number_of_blocks);

    do{
        cudaMemcpy(d_gbd, grid_box_dsv, memSizef, cudaMemcpyHostToDevice);

        compute_updated_DSV<<<BlocksPerGrid, ThreadsPerBlock>>>(d_gli,d_glc,d_gri,d_grc,
        	d_gti,d_gtc,d_gbi,d_gbc,d_gbd,d_gh,d_gw,d_gp,n_l,n_r,n_t,n_b,number_of_boxes,affect_rate,updsv);

        // get results
        cudaMemcpy(updated_DSVs, updsv, memSizef, cudaMemcpyDeviceToHost);

        /* Update the grid. */
        result = commit_dsv_update();

        assert(result == SUCCESS);

        convergence_iteration = convergence_iteration + 1;
        /* Check if the convergence reached. */
        result = check_for_convergence();
    }while(result != SUCCESS);

    //release GPU cache
    cudaFree(d_gli);
    cudaFree(d_glc);
    cudaFree(d_gri);
    cudaFree(d_grc);
    cudaFree(d_gti);   
    cudaFree(d_gtc);   
    cudaFree(d_gbi);   
    cudaFree(d_gbc);   

    cudaFree(d_gbd);
    cudaFree(d_gh);
    cudaFree(d_gw);
    cudaFree(d_gp);

    cudaFree(n_l);   
    cudaFree(n_r);   
    cudaFree(n_t);
    cudaFree(n_b);

    cudaFree(updsv);

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

    /* Read the number of grid and create grid. */
    scanf("%d",&number_of_boxes);

    grid = (box*)malloc(sizeof(box)*number_of_boxes);
    updated_DSVs = (float*)malloc(sizeof(float)*number_of_boxes);
    neighbors = (int**)malloc(sizeof(int*)*number_of_boxes);
    for(int i = 0;i < number_of_boxes;i++){
        neighbors[i] = (int*)malloc(sizeof(int)*4);
    }

    /* Read the input file. */
    result = read_input_file();

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

    printf("Threads requested: %d.\n", number_of_threads);
    printf("Dissipation converged in %d iterations.\n",convergence_iteration);
    printf("With max DSV = %f and min DSV = %f.\n", max_dsv, min_dsv);
    printf("Affect rate = %f; Epsilon: %f.\n",affect_rate, epsilon);
    printf("Elapsed convergence loop time (clock) : %d.\n", (int)(clock2-clock1));
    printf("Elapsed convergence loop time (time) : %d.\n", (int)difftime(time2, time1));
    printf("Elapsed convergence loop time (chrono) : %2f.\n", time_taken);
    printf("*******************************************************************\n");

    return 0;
}


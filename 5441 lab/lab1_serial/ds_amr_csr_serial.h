/*
 * CSE 5441 : Lab 1
 * Sheng Ding
 * ding.853@osu.edu
 * The Ohio State University
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>

#define SUCCESS 1
#define ERROR 0
#define HORIZONTAL 1
#define VERTICAL 0

#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)

/*
 * Structure to store neighbor’s information.
 */
typedef struct neighbor {
    int id; /* Stores the neighbor’s ID. */
    int contact_distance; /* Stores the contact distance between two neighbors. */
} neighbor;

/*
 * Structure to store properties of a box.
 */
typedef struct box {
    /*
     * position current box on underlying
     * co-ordinated grid.
     */
    int upper_left_y;
    int upper_left_x;

    /* Height, width and perimeter of the current box. */
    int height;
    int width;
    int perimeter;

    /* Current box DSV (temperature). */
    float box_dsv;

    /* The neighbors of the current box. */
    neighbor *left;
    neighbor *right;
    neighbor *top;
    neighbor *bottom;
} box;

int commit_dsv_update(box grid[], float updated_DSVs[], int number_of_boxes, float *max, float *min);
int check_for_convergence(float epsilon, float max, float min);

#ifndef STRUCTURE_H
#define STRUCTURE_H

#define D2Q9_DIM 2
#define D2Q9_Q 9

typedef struct
{
    int dimensions;
    int velocity_number;
    float weights[D2Q9_Q]; /* 9 weights (1 per velocity) */

    /* Matrix [Q][DIM]: 9 rows, 2 cols */
    int velocities_by_dir[D2Q9_Q][D2Q9_DIM];

    /* Matrix [DIM][Q]: 2 rows, 9 cols */
    int velocities_by_dim[D2Q9_DIM][D2Q9_Q];

    int opposite[D2Q9_Q];

} Structure;

/*
D2Q9 lattice
* 6 2 5
* 3 0 1
* 7 4 8
*/

extern const Structure D2Q9;

#endif

/*

#define D3Q27_DIM 3
#define D3Q27_Q   27

typedef struct {
    int dimensions;
    int velocity_number;
    float weights[D3Q27_Q];
    int velocities_by_dir[D3Q27_Q][D3Q27_DIM];
    int velocities_by_dim[D3Q27_DIM][D3Q27_Q];
    int opposite[D3Q27_Q];
} Structure;


*/
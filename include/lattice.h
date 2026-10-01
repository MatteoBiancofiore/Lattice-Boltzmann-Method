#ifndef LATTICE_H
#define LATTICE_H

#include "node.h"
#include "structure.h"

typedef struct
{
    int width;                  /* Dimensione orizzontale della griglia (Nx) */
    int height;                 /* Dimensione verticale della griglia (Ny) */
    Node *grid;                 /* Matrice 2D di celle [height][width], array di Node* */
    const Structure *structure; /* Puntatore alla struttura (es. &D2Q9) */

    /* Parametri fisici*/
    float tau;        /* Tempo di rilassamento collisionale */
    float omega;      /* Frequenza di rilassamento: 1.0f / tau */
    int current_step; /* Iterazione temporale corrente */
} Lattice;

Lattice *lattice_init(int width, int height, const Structure *structure, float tau);

void lattice_free(Lattice *lattice);

void lattice_collide(Lattice *lattice); // wrapper per equilibrium_collision

void lattice_stream(Lattice *lattice); // wrapper per streaming

void lattice_save_csv(Lattice *lattice, char *filename);

// void lattice_apply_boundaries(Lattice *lattice, float u_lid); // wrapper per

#endif
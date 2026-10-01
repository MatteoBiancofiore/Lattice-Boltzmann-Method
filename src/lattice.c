#include "../include/lattice.h"
#include <stdio.h>
#include <stdlib.h>

#define IDX(x, y, width) ((y) * (width) + (x))

Lattice *lattice_init(int width, int height, const Structure *structure,
                      float tau)
{
    Lattice *lat = (Lattice *)malloc(sizeof(Lattice));
    if (lat == NULL)
    {
        printf("Malloc error for Lattice");
        exit(1);
    }

    lat->width = width;
    lat->height = height;
    lat->structure = structure;
    lat->tau = tau;
    lat->omega = 1.0f / tau;
    lat->current_step = 0;

    // Singola allocazione 1D contigua per tutti i nodi
    lat->grid = (Node *)malloc(width * height * sizeof(Node));
    if (lat->grid == NULL)
    {
        free(lat);
        printf("Malloc error for Grid");
        exit(1);
    }

    // Inizializzazione di ciascun nodo direttamente in-place
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int pos[2] = {x, y};
            int bnd[2] = {0, 0};

            /* Rilevamento bordi */
            if (x == 0)
                bnd[0] = -1; // Sinistra
            else if (x == width - 1)
                bnd[0] = 1; // Destra

            if (y == 0)
                bnd[1] = -1; // Sotto
            else if (y == height - 1)
                bnd[1] = 1; // Sopra

            int idx = IDX(x, y, width);

            // Inizializza direttamente la memoria nella posizione idx
            node_init(&lat->grid[idx], pos, bnd, false);

            // Inizializza le distribuzioni f_i al valore di equilibrio
            init_eq(&lat->grid[idx], structure, 0);
        }
    }

    return lat;
}

void lattice_destroy(Lattice *lat)
{
    if (lat == NULL)
        return;

    // Si libera l'unico blocco contiguo della griglia
    if (lat->grid != NULL)
    {
        free(lat->grid);
    }
    free(lat);
}

void lattice_collide(Lattice *lattice)
{
    for (int i = 0; i < lattice->width * lattice->height; i++)
    {
        Node *n = &lattice->grid[i];
        if (!n->is_obstacle)
        {
            node_update_macro(n, lattice->structure);
            equilibrium_collision(n, lattice->structure, lattice->omega);
        }
    }
}

void lattice_stream(Lattice *lattice)
{
    for (int i = 0; i < lattice->width * lattice->height; i++)
    {
        Node *n = &lattice->grid[i];
        streaming(n, lattice, lattice->structure);
    }
}

void lattice_step(Lattice *lattice, float lid_velocity)
{
    // 1. Collisione su tutti i nodi fluidi
    lattice_collide(lattice);

    // 2. Propagazione e rimbalzo
    lattice_stream(lattice);

    // 3. Applicazione condizioni al contorno dinamiche (es. profilo velocità inlet)
    for (int i = 0; i < lattice->width * lattice->height; i++)
    {
        // set_boundary_velocity(&lattice->grid[i], lattice, lid_velocity);
    }

    lattice->current_step++;

    if (lattice->current_step % 100 == 0)
    {
        char *filename;
        snprintf(filename, sizeof(filename), "../out/output_%06.csv", lattice->current_step);
        lattice_save_csv(lattice, filename);
    }
}

// save parameters on file
void lattice_save_csv(Lattice *lattice, char *filename)
{
    FILE *fp = fopen(filename, "w");

    if (fp == NULL)
    {
        printf("Error while opening %s", filename);
        exit(1);
    }

    // Intestazione

    fprintf(fp, "x,y,u_x,u_y,rho,is_obstacle\n");

    for (int j = 0; j < lattice->height; j++)
    {
        for (int i = 0; i < lattice->width; i++)
        {
            int idx = IDX(i, j, lattice->width);

            Node *n = &lattice->grid[idx];

            fprintf(fp, "%d,%d,%.6f,%.6f,%.6f,%d\n", i, j, n->macroscopic_velocity[0], n->macroscopic_velocity[1], n->rho, n->is_obstacle);
        }
    }
    fclose(fp);
}
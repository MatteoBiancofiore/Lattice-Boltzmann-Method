#ifndef NODE_H
#define NODE_H

#include <stdbool.h>
#include <math.h>
#include "structure.h"

typedef struct
{
  // macroscopic density.
  float rho;

  // particle distribution function.
  float distribution[D2Q9_Q];
  float new_distribution[D2Q9_Q];

  // macroscopic velocity of the cell.
  float macroscopic_velocity[D2Q9_DIM];
  int boundary[D2Q9_DIM];
  int position[D2Q9_DIM];

  bool is_obstacle;
} Node;

void node_init(Node *node, const int position[D2Q9_DIM], const int boundary[D2Q9_DIM], bool is_obstacle);

void node_update_macro(Node *node, const Structure *structure);

void equilibrium_collision(Node *node, const Structure *structure, float omega);

void init_eq(Node *node, const Structure *structure, const int problem_type);

// utils
float scalarProduct_2(const float a[], const float b[]);
float scalarProduct_3(const float a[], const float b[], const float c[]);

#endif // !NODE_H

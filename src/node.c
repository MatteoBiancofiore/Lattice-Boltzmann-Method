#include "../include/node.h"

// Utils
float scalar_product_2(const float a[], const float b[])
{
  return a[0] * b[0] + a[1] * b[1];
}

float scalar_product_3(const float a[], const float b[], const float c[])
{
  return a[0] * b[0] * c[0] + a[1] * b[1] * c[1];
}

// Node
void node_init(Node *node, const int position[2], const int boundary[2],
               bool is_obstacle)
{
  node->rho = 1.0f;

  for (int i = 0; i < D2Q9_Q; i++)
  {
    node->distribution[i] = 0.0f;
    node->new_distribution[i] = 0.0f;
  }

  for (int i = 0; i < D2Q9_DIM; i++)
  {
    node->macroscopic_velocity[i] = 0.0f;
    node->position[i] = position[i];
    node->boundary[i] = boundary[i];
  }

  node->is_obstacle = is_obstacle;
}

void node_update_macro(Node *node, const Structure *structure)
{
  if (node->is_obstacle)
    return;

  float rho = 0.0f;

  /* rho = sum_i f_i */
  for (int i = 0; i < structure->velocity_number; i++)
  {
    rho += node->distribution[i];
  }

  node->rho = rho;

  /* u = (1/rho) * sum_i f_i * c_i */
  for (int d = 0; d < structure->dimensions; d++)
  {
    float velocity = 0.0f;

    for (int i = 0; i < structure->velocity_number; i++)
    {
      velocity += node->distribution[i] * structure->velocities_by_dim[i][d];
    }

    node->macroscopic_velocity[d] = velocity / rho;
  }
}

// Node collision with BGK: f_i_new = f_i + omega * (f_i^eq - f_i)
void equilibrium_collision(Node *node, const Structure *structure,
                           float omega)
{
  if (node->is_obstacle)
    return;

  float f_eq[D2Q9_Q];

  float u_square =
      node->macroscopic_velocity[0] * node->macroscopic_velocity[0] +
      node->macroscopic_velocity[1] * node->macroscopic_velocity[1];

  for (int i = 0; i < structure->velocity_number; i++)
  {
    float cu =
        structure->velocities_by_dir[i][0] * node->macroscopic_velocity[0] +
        structure->velocities_by_dir[i][1] * node->macroscopic_velocity[1];

    f_eq[i] = structure->weights[i] * node->rho *
              (1.0f + 3.0f * cu + 4.5f * cu * cu - 1.5f * u_square);
  }

  // BGK collision
  for (int i = 0; i < structure->velocity_number; i++)
  {
    node->new_distribution[i] =
        node->distribution[i] + omega * (f_eq[i] - node->distribution[i]);
  }
}

// Initialize the distribution f to its equilibrium value.
// Equilibrium distribution:
// f_i^eq = omega_i * rho * [1 + 3(c_i · u) + (9/2)(c_i · u)^2 - (3/2)|u|^2]
void init_eq(Node *node, const Structure *structure, const int problem_type)
{
  if (node->is_obstacle && problem_type == 2)
  {
    node->macroscopic_velocity[0] = NAN;
    node->macroscopic_velocity[1] = NAN;
    return;
  }

  const float velocity_squared =
      scalar_product_2(node->macroscopic_velocity, node->macroscopic_velocity);

  const float three_halves_velocity_squared = 1.5f * velocity_squared;

  for (int i = 0; i < structure->velocity_number; i++)
  {
    const float three_times_velocity_projection =
        3.0f * scalar_product_2(structure->velocities_by_dir[i],
                                node->macroscopic_velocity);

    node->distribution[i] = structure->weights[i] * node->rho *
                            (1.0f + three_times_velocity_projection +
                             0.5f * three_times_velocity_projection *
                                 three_times_velocity_projection -
                             three_halves_velocity_squared);

    node->new_distribution[i] = node->distribution[i];
  }
}

// Stream the post-collision distribution to adjacent nodes
void streaming(Node *node, Lattice *lattice, const Structure *structure)
{
  if (node->is_obstacle)
    return;

  // The rest particle does not move
  node->distribution[0] = node->new_distribution[0];

  for (int i = 1; i < structure->velocity_number; i++)
  {
    // Calculate the position of the destination node
    int new_position[D2Q9_DIM];

    for (int d = 0; d < structure->dimensions; d++)
    {
      new_position[d] = node->position[d] + structure->velocities_by_dir[i][d];
    }

    // Check that the destination is inside the lattice
    if (new_position[0] >= 0 && new_position[0] < lattice->width &&
        new_position[1] >= 0 && new_position[1] < lattice->height)
    {
      Node *destination =
          &lattice->grid[IDX(new_position[0], new_position[1], lattice->width)];

      // Do not stream into obstacle nodes
      if (!destination->is_obstacle)
      {
        destination->distribution[i] = node->new_distribution[i];
      }
    }
    else
    {

      // Andrebbe implementata la logica di bounce-back, se c'è un ostacolo rimbalza indietro ()
        }
  }
}

// Set the macroscopic velocity at the domain boundaries
void set_boundary_velocity(Node *node, const Lattice *lattice,
                           const float lid_velocity)
{
  if (node->is_obstacle)
    return;

  const int x_len = lattice->width;
  const int y_len = lattice->height;
  const int problem_type = 2;

  // Set velocity to zero at all walls.
  if (node->position[0] == 0 || node->position[1] == 0 ||
      node->position[0] == x_len - 1 || node->position[1] == y_len - 1)
  {
    node->macroscopic_velocity[0] = 0.0f;
    node->macroscopic_velocity[1] = 0.0f;
  }

  switch (problem_type)
  {
  case 1:
    // Lid driven cavity: the top wall moves with velocity.
    if (node->position[1] == 0)
    {
      node->macroscopic_velocity[0] = lid_velocity;
    }

    break;

  case 2:
    // Channel flow: impose a parabolic velocity profile at the inlet.
    if (node->position[0] == 0)
    {
      const float half_dimension = (float)y_len / 2.0f;

      const float normalized_position =
          ((float)node->position[1] / half_dimension) - 1.0f;

      const float parabolic_profile =
          1.0f - normalized_position * normalized_position;

      node->macroscopic_velocity[0] = lid_velocity * parabolic_profile;
    }

    break;

  default:
    break;
  }
}

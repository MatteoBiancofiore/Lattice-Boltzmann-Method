#include "../include/node.h"
#include "../include/structure.h"

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
void node_init(Node *node, const int position[2], const int boundary[2], bool is_obstacle)
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
void equilibrium_collision(Node *node, const Structure *structure, float omega)
{
  if (node->is_obstacle)
    return;

  float f_eq[D2Q9_Q];

  float u_square = node->macroscopic_velocity[0] * node->macroscopic_velocity[0] +
                   node->macroscopic_velocity[1] * node->macroscopic_velocity[1];

  for (int i = 0; i < structure->velocity_number; i++)
  {
    float cu = structure->velocities_by_dir[i][0] *
                   node->macroscopic_velocity[0] +
               structure->velocities_by_dir[i][1] *
                   node->macroscopic_velocity[1];

    f_eq[i] = structure->weights[i] * node->rho * (1.0f + 3.0f * cu + 4.5f * cu * cu - 1.5f * u_square);
  }

  // BGK collision
  for (int i = 0; i < structure->velocity_number; i++)
  {
    node->new_distribution[i] = node->distribution[i] + omega * (f_eq[i] - node->distribution[i]);
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

  // |u|^2
  const float velocity_squared =
      scalar_product_2(
          node->macroscopic_velocity,
          node->macroscopic_velocity);

  // (3/2) |u|^2
  const float three_halves_velocity_squared =
      1.5f * velocity_squared;

  for (int i = 0; i < structure->velocity_number; i++)
  {
    // 3(c_i · u)
    const float three_times_velocity_projection =
        3.0f * scalar_product_2(
                   structure->velocities_by_dir[i],
                   node->macroscopic_velocity);

    node->distribution[i] =
        structure->weights[i] *
        node->rho *
        (1.0f + three_times_velocity_projection + 0.5f * three_times_velocity_projection * three_times_velocity_projection - three_halves_velocity_squared);
  }
}

#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void apply_velocity_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_u,
    const double *boundary_conditions_u_values,
    const size_t *boundary_conditions_u_ids,
    double *velocities
);

void apply_pressure_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_p,
    const double *boundary_conditions_p_values,
    const size_t *boundary_conditions_p_ids,
    double *pressures
);

#endif
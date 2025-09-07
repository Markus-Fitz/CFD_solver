#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void apply_velocity_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_u,
    double *velocities
);

void apply_pressure_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_p,
    double *pressures
);

#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "solver.h"

void apply_velocity_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_u,
    const double *boundary_conditions_u_values,
    const size_t *boundary_conditions_u_ids,
    double *velocities
){
    for(size_t i = 0; i < 3*num_cells; i++){
        if(boundary_conditions_u[i] == 1){
            velocities[i] = boundary_conditions_u_values[i];
        }
        if(boundary_conditions_u[i] == 2){
            velocities[i] = velocities[boundary_conditions_u_ids[i]];
        }
    }
}

void apply_pressure_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_p,
    const double *boundary_conditions_p_values,
    const size_t *boundary_conditions_p_ids,
    double *pressures
){
    for(size_t i = 0; i < num_cells; i++){
        if(boundary_conditions_p[i] == 1){
            pressures[i] = boundary_conditions_p_values[i];
        }
        if(boundary_conditions_p[i] == 2){
            pressures[i] = pressures[boundary_conditions_p_ids[i]];
        }
    }
}
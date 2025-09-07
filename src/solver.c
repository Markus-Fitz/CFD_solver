#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "solver.h"

void apply_velocity_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_u,
    double *velocities
){
    for(size_t i = 0; i < 3*num_cells; i++){
        //BCs
    }
}

void apply_pressure_BCs(
    size_t num_cells,
    const uint8_t *boundary_conditions_p,
    double *pressures
){
    for(size_t i = 0; i < num_cells; i++){
        //BCs
    }
}
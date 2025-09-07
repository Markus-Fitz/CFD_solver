#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "output.h"
#include "solver.h"

int main() {

    size_t nx = 31, ny = 51, nz = 1;    // size of the field of cell centers
    size_t num_cells = nx * ny * nz; // number of cells
    double dx = 0.1, dy = 0.1, dz = 0.1;
    size_t num_points = (nx + 1)*(ny + 1)*(nz + 1); // number of mesh node points

    double dt = 1;  // (maximum) timestep in seconds; possibly adaptive timestep later on by examining cell flow velocities?
    int num_steps = 50; // number of timesteps

    double nu = 0.0001; // viscosity
    double rho = 1.0;   // density

    // coordinate array for cell centers
    double *cell_centers = malloc(num_cells * 3 * sizeof(double));
    if (cell_centers == NULL) return 1;

    // array to denote neighbouring cell IDs for each cell; at each cube face the cell has 4 neighbours (= 24 neighbours for whole cell) "worst case" (AMR)
    size_t *cell_interfaces = malloc(num_cells * 24 * sizeof(size_t));
    if (cell_interfaces == NULL) return 1;

    // initializing cell_interface array
    for(size_t i = 0; i < num_cells*24; i++){
        cell_interfaces[i] = 0;
    }

    // coordinate array for mesh points (n+1 in each direction, staggered by 0.5 cell legth)
    double *mesh_point_coordinates = malloc(num_points * 3 * sizeof(double));
    if (mesh_point_coordinates == NULL) return 1;

    // array to denote edge point IDs (in mesh point array) for each cell
    size_t *cell_points = malloc(num_cells * 8 * sizeof(size_t));
    if (cell_points == NULL) return 1;

    // array for (updated) pressure at each cell center
    double *pressure = malloc(num_cells * sizeof(double));
    if (pressure == NULL) return 1;

    // array for the three velocity components at each cell center
    double *velocities = malloc(num_cells * 3 * sizeof(double));
    if (velocities == NULL) return 1;

    // array for the three tentative velocity components at each cell center
    double *tent_velocities = malloc(num_cells * 3 * sizeof(double));
    if (tent_velocities == NULL) return 1;

    // TODO: combine u and p BCs into one array?
    // array for velocity boundary conditions:
    // 0: no BC
    // 1: const_velocity: velocity is set to constant value e.g. zero -> see velocity_BC_val array
    // 2: copy velocity from other cell face: velocity is set to velocity value of ohter cell ID in bc_id-array
    uint8_t *boundary_conditions_u = malloc(num_cells * 3 * sizeof(uint8_t));
    if(boundary_conditions_u == NULL) return 1;

    // initialize BCs
    for(size_t i = 0; i < num_cells*3; i++){
        boundary_conditions_u[i] = 0;
    }

    // array for pressure boundary conditions:
    // 0: no BC
    // 1: const pressure
    // 2: copy pressure from other cell center
    uint8_t *boundary_conditions_p = malloc(num_cells * sizeof(uint8_t));
    if(boundary_conditions_p == NULL) return 1;

    // initialize BCs
    for(size_t i = 0; i < num_cells; i++){
        boundary_conditions_p[i] = 0;
    }

    // generate initial pressure and velocity fields
    for(size_t z = 0; z < nz; z++){
        for(size_t y = 0; y < ny; y++) {
            for(size_t x = 0; x < nx; x++) {

                // coordinates of cell-center in question
                double coord_x = x*dx;
                double coord_y = y*dy;
                double coord_z = z*dz;

                // index of current cell: after nx cells: y gets iterated by one, after nx*ny cells, z gets iterated by one
                size_t id = x + nx*y + nx*ny*z;

                // write down coordinate of cell into array
                cell_centers[3*id] = coord_x;       // x-coordinate
                cell_centers[3*id + 1] = coord_y;   // y-coordinate
                cell_centers[3*id + 2] = coord_z;   // z-coordinate

                // lower left corner ID of the voxel
                size_t p0_ID = x + y*(nx+1) + z*(nx+1)*(ny+1);

                // for each ID list all eight points connecting cell; one row in y-dir equals ID-offset of nx; one row in z-dir equals ID-offset of nx*ny
                cell_points[8*id] = p0_ID;
                cell_points[8*id + 1] = p0_ID + 1;
                cell_points[8*id + 2] = p0_ID + (nx+1);
                cell_points[8*id + 3] = p0_ID + (nx+1) + 1;
                cell_points[8*id + 4] = p0_ID + (nx+1)*(ny+1);
                cell_points[8*id + 5] = p0_ID + (nx+1)*(ny+1) + 1;
                cell_points[8*id + 6] = p0_ID + (nx+1)*(ny+1) + (nx+1);
                cell_points[8*id + 7] = p0_ID + (nx+1)*(ny+1) + (nx+1) + 1;

                // generate cell interface array; ordering: x+ neighbours (x4), y+ neighbours (x4), z+ neighbours (x4), then x-, y-, z-
                if(x < nx - 1){
                    cell_interfaces[24*id] = id + 1;
                }else{
                    cell_interfaces[24*id] = 0;    // last cell in x-direction -> no neighbour in x+ direction
                }

                // y+ direction
                if(y < ny - 1){
                    cell_interfaces[24*id + 4] = id + nx;
                }else{
                    cell_interfaces[24*id + 4] = 0;
                }

                // z+ direction
                if(z < nz - 1){
                    cell_interfaces[24*id + 8] = id + nx*ny;
                }else{
                    cell_interfaces[24*id + 8] = 0;
                }

                // x- direction
                if(x > 0){
                    cell_interfaces[24*id + 12] = id-1;
                }else{
                    cell_interfaces[24*id + 12] = 0;
                }

                // y- direction
                if(y > 0){
                    cell_interfaces[24*id + 16] = id - nx;
                }else{
                    cell_interfaces[24*id + 16] = 0;
                }

                // z- direction
                if(z > 0){
                    cell_interfaces[24*id + 20] = id - nx*ny;
                }else{
                    cell_interfaces[24*id + 20] = 0;
                }


                // initial pressure field
                pressure[id] = 0;                 //(coord_x - 5)*(coord_x - 5) + (coord_y - 5)*(coord_y - 5) + (coord_z - 5)*(coord_z - 5);     // parabolic field centered around center
                
                // initial velocity field
                if(y == ny-1) velocities[3*id] = 1.0;  // velocities in x-dir; lid-driven cavity, u = 1.0 for top cells
                velocities[3*id+1] = 0;             // velocities in y-dir
                velocities[3*id+2] = 0;             // velocities in z-dir

            }
        }
    }

    // generate point mesh
    for(size_t z = 0; z < nz+1; z++){
        for(size_t y = 0; y < ny+1; y++){
            for(size_t x = 0; x < nx+1; x++){

                // coordinate of cell center in question
                double coord_x = x*dx;
                double coord_y = y*dy;
                double coord_z = z*dz;
                
                // point grid has one more coordinates in x-direction!
                size_t id = x + (nx+1)*y + (nx+1)*(ny+1)*z;

                // coordinates of mesh grid points - grid starts 0.5*dx/dy/dz "before" first cell center
                mesh_point_coordinates[3*id] = coord_x - 0.5*dx;
                mesh_point_coordinates[3*id + 1] = coord_y - 0.5*dy;
                mesh_point_coordinates[3*id + 2] = coord_z - 0.5*dz;

            }
        }
    }

    // main simulation loop
    for(int step = 0; step < num_steps; step++){
        
        double t = step*dt;

        // iterating update equations over all 
        for(size_t cell_id = 0; cell_id < num_cells; cell_id++){

            // boundary conditions on velocities
            
            // calculate tentative velocities


            // enforce boundary conditions on tentative velocities

            // calculate pressures and enforce pressure boundary conditions

            // write output into vtu file
            write_vtu("output.vtu", num_cells, num_points, mesh_point_coordinates, cell_points, pressure, velocities, t);
        }
    }

    free(cell_centers);
    free(mesh_point_coordinates);
    free(cell_points);
    free(pressure);
    free(velocities);

    printf("Wrote output.vtu — open in ParaView!\n");

    return 0;
}
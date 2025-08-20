#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "output.h"

int main() {

    size_t nx = 11, ny = 11, nz = 11;    // size of the field of cell centers
    size_t num_cells = nx * ny * nz; // number of cells
    double dx = 1.0, dy = 1.0, dz = 1.0;
    size_t num_points = (nx + 1)*(ny + 1)*(nz + 1);

    double *cell_centers = malloc(num_cells * 3 * sizeof(double));
    if (cell_centers == NULL) return 1;

    double *mesh_point_coordinates = malloc(num_points * 3 * sizeof(double));
    if (mesh_point_coordinates == NULL) return 1;

    size_t *cell_points = malloc(num_cells * 8 * sizeof(size_t));
    if (cell_points == NULL) return 1;

    double *pressure = malloc(num_cells * sizeof(double));
    if (pressure == NULL) return 1;

    double *velocities = malloc(num_cells * 3 * sizeof(double));
    if (velocities == NULL) return 1;

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

                // initial pressure field
                pressure[id] = (coord_x - 5)*(coord_x - 5) + (coord_y - 5)*(coord_y - 5) + (coord_z - 5)*(coord_z - 5);     // parabolic field centered around center
                
                // initial velocity field
                velocities[3*id] = coord_x - 5;         // velocities in x-dir
                velocities[3*id+1] = coord_y - 5;       // velocities in y-dir
                velocities[3*id+2] = coord_z - 5;       // velocities in z-dir
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

    // generate conections of grid nodes
    for(size_t z = 0; z < nz; z++){
        for(size_t y = 0; y < ny; y++){
            for(size_t x = 0; x < nx; x++){
                
                size_t id = x + nx*y + nx*ny*z;

                // for each id list all eight points connecting cell; one row in y-dir equals ID-offset of nx; one row in z-dir equals ID-offset of nx*ny
                cell_points[8*id] = id;
                cell_points[8*id + 1] = id + 1;
                cell_points[8*id + 2] = id + nx + 1;
                cell_points[8*id + 3] = id + nx;
                cell_points[8*id + 4] = id + nx*ny;
                cell_points[8*id + 5] = id + nx*ny + 1;
                cell_points[8*id + 6] = id + nx*ny + nx + 1;
                cell_points[8*id + 7] = id + nx*ny + nx;

            }
        }
    }

    write_vtu("output.vtu", num_cells, num_points, mesh_point_coordinates, cell_points, pressure, velocities);
    free(cell_centers);
    free(mesh_point_coordinates);
    free(cell_points);
    free(pressure);
    free(velocities);

    printf("Wrote output.vtu — open in ParaView!\n");

    return 0;
}
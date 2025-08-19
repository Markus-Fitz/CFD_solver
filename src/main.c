#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "output.h"

int main() {

    size_t nx = 11, ny = 11;    // size of the field of cell centers
    size_t num_cells = nx * ny; // number of cells
    double dx = 1.0, dy = 1.0;
    int dim = 2;
    size_t num_points = (nx + 1)*(ny + 1);

    double *cell_centers = malloc(nx * ny * 2 * sizeof(double));
    if (cell_centers == NULL) return 1;

    double *mesh_point_coordinates = malloc((nx + 1)*(ny + 1) * 2 * sizeof(double));
    if (mesh_point_coordinates == NULL) return 1;

    size_t *cell_points = malloc(nx * ny * 4 * sizeof(size_t));
    if (cell_points == NULL) return 1;

    double *pressure = malloc(nx * ny * sizeof(double));
    if (pressure == NULL) return 1;

    double *velocities = malloc(nx * ny * 2 * sizeof(double));
    if (velocities == NULL) return 1;

    // generate debug dummy fields for pressure and velocity
    for(size_t j = 0; j < ny; j++) {
        for(size_t i = 0; i < nx; i++){

            double x = i*dx;
            double y = j*dy;

            size_t id = j*nx + i;   // index with "stride" for 2D arrays

            cell_centers[2*id] = x;  // x-coordinate
            cell_centers[2*id + 1] = y;  // y-coordinate

            pressure[j*ny + i] = (i*dx - 5)*(i*dx - 5) + (j*dy - 5)*(j*dy - 5); // parabolic field centered around middle

            velocities[2*id] = i*dx - 5.5; // velocities pointing outward of center in x-dir
            velocities[2*id+1] = j*dy - 5.5;    //velocities in y-dir
        }
    }

    // generate point mesh
    for(size_t j = 0; j < ny+1; j++){
        for(size_t i = 0; i < nx+1; i++){

            double x = i*dx;
            double y = j*dy;

            size_t id = j*(nx+1) + i;   // point grid has one more coordinates in x-direction!

            mesh_point_coordinates[2*id] = x - 0.5*dx;  // grid starts 0.5*dx "before" first cell center
            mesh_point_coordinates[2*id + 1] = y - 0.5*dy;  // grid starts 0.5*dy "before" first cell center

        }
    }

    // generate conections
    for(size_t j = 0; j < ny; j++){
        for(size_t i = 0; i < nx; i++){
            
            size_t id = j*nx + i;

            cell_points[4*id] = j * (nx + 1) + i;
            cell_points[4*id + 1] = j * (nx + 1) + i + 1;
            cell_points[4*id + 2] = (j + 1) * (nx + 1) + i + 1;
            cell_points[4*id + 3] = (j + 1) * (nx + 1) + i;

        }
    }

    write_vtu("output.vtu", dim, num_cells, num_points, mesh_point_coordinates, cell_points, pressure, velocities);
    free(cell_centers);
    free(mesh_point_coordinates);
    free(cell_points);
    free(pressure);
    free(velocities);

    printf("Wrote output.vtu — open in ParaView!\n");

    return 0;
}
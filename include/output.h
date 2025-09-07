#ifndef OUTPUT_H
#define OUTPUT_H

#include <stddef.h> // for size_t

void write_vtu(
    const char *filename,       // name of output file
    size_t num_cells,           // number of cells for extraction of values from arrays
    size_t num_points,          // number of points making up the edges of the cells -> the coordinates for the mesh
    const double *mesh_point_coordinates,   // array of mesh point coordinates making up the corners of the cells
    const size_t *cell_points,  // array of indices which for each cell-ID point to the point-IDs making up the cells
    const double *pressure,     // array of pressures at the cell centers indexed by the cell-IDs
    const double *velocities,    // array of flow-velocities at the cell centers indexed by the cell-IDs
    double time
);

#endif
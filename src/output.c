#include <stdio.h>
#include "output.h"

void write_vtu(
    const char *filename, 
    size_t num_cells,
    size_t num_points,
    const double *mesh_point_coordinates, 
    const size_t *cell_points, 
    const double *pressure, 
    const double *velocities
    double time) {
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        return;
    }

    fprintf(fp,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        "  <UnstructuredGrid>\n"
        "    <Piece NumberOfPoints=\"%zu\" NumberOfCells=\"%zu\">\n",
        num_points, num_cells);
    
    // Write points
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(size_t i = 0; i < num_points; i++){
        fprintf(fp, "%f %f %f \n", mesh_point_coordinates[i*3], mesh_point_coordinates[i*3 + 1], mesh_point_coordinates[i*3 + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // Write cells (voxels)
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for(size_t i = 0; i < num_cells; i++){
        fprintf(fp, "%zu %zu %zu %zu %zu %zu %zu %zu\n", cell_points[i*8], cell_points[i*8 + 1], cell_points[i*8 + 2], cell_points[i*8 + 3], cell_points[i*8 + 4], cell_points[i*8 + 5], cell_points[i*8 + 6], cell_points[i*8 + 7]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (size_t c = 1; c <= num_cells; c++) {
        fprintf(fp, "%zu\n", c*8);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (size_t c = 0; c < num_cells; c++) {
        fprintf(fp, "11\n"); // VTK_VOXEL = 11
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // Write pressure as CELL_DATA
    fprintf(fp, "      <CellData Scalars=\"pressure\">\n");
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">\n");
    for(size_t i = 0; i < num_cells; i++){
        fprintf(fp, "%f\n", pressure[i]);
    }
    fprintf(fp, "        </DataArray>\n");

    // Write Velocities as CELL_DATA
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (size_t c = 0; c < num_cells; c++) {
        fprintf(fp, "      %f %f %f\n", velocities[3*c], velocities[3*c + 1], velocities[3*c + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </CellData>\n");

    fprintf(fp,
        "    </Piece>\n"
        "  </UnstructuredGrid>\n"
        "</VTKFile>\n");

    fclose(fp);
}

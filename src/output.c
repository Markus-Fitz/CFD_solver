#include <stdio.h>
#include "output.h"

void write_vtu(
    const char *filename, 
    int dim,
    size_t num_cells,
    size_t num_points,
    const double *mesh_point_coordinates, 
    const size_t *connections, 
    const double *pressure, 
    const double *velocities) {
    
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
        if(dim == 2) fprintf(fp, "%f %f %f \n", mesh_point_coordinates[i*2], mesh_point_coordinates[i*2 + 1], 0.0);
        if(dim == 3) fprintf(fp, "%f %f %f \n", mesh_point_coordinates[i*3], mesh_point_coordinates[i*3 + 1], mesh_point_coordinates[i*3 + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // Write cells (quads / hexagons)
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

    for(size_t i = 0; i < num_cells; i++){
        if(dim == 2) fprintf(fp, "%d %d %d %d\n", connections[i*4], connections[i*4 + 1], connections[i*4 + 2], connections[i*4 + 3]);
        if(dim == 3) fprintf(fp, "%d %d %d %d\n", connections[i*8], connections[i*8 + 1], connections[i*8 + 2], connections[i*8 + 3], connections[i*8 + 4], connections[i*8 + 5], connections[i*8 + 6], connections[i*8 + 7]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (size_t c = 1; c <= num_cells; c++) {
        if(dim == 2) fprintf(fp, "%zu\n", c*4);
        if(dim == 3) fprintf(fp, "%zu\n", c*8);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (size_t c = 0; c < num_cells; c++) {
        if(dim == 2) fprintf(fp, "9\n");  // VTK_QUAD = 9
        if(dim == 3) fprintf(fp, "12\n"); // VTK_HEX = 12
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
    fprintf(fp, "      </CellData>\n");

    fprintf(fp,
        "    </Piece>\n"
        "  </UnstructuredGrid>\n"
        "</VTKFile>\n");

    fclose(fp);
}

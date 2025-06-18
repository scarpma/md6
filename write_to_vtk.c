#include "defs.h"

void writePointCloudToVTK(const char *filename, const vec *points, const vec *velocity, int numPoints) {
    FILE *fp = NULL;

    // Open the file for writing
    fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }

    // --- VTK Legacy File Header ---
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Point Cloud\n"); // Title
    fprintf(fp, "ASCII\n");      // Data type (ASCII for human readability)

    // --- Dataset Structure ---
    fprintf(fp, "DATASET POLYDATA\n"); // We're representing points, which can be thought of as a degenerate polydata

    // --- Points Section ---
    fprintf(fp, "POINTS %d float\n", numPoints); // Number of points and data type
    for (int i = 0; i < numPoints; i++) {
        fprintf(fp, "%.4e %.4e %.4e\n", points[i].x, points[i].y, points[i].z);
    }

    // // --- Vertices Section (to make each point a distinct vertex) ---
    // // This tells VTK to render each point individually.
    // // The format is: num_vertices num_indices_total
    // // For each vertex: 1 index_of_point
    // fprintf(fp, "VERTICES %d %d\n", numPoints, numPoints * 2); // numPoints vertices, each with 1 index + its own count
    // for (int i = 0; i < numPoints; i++) {
    //     fprintf(fp, "1 %d\n", i); // '1' indicates one point in this vertex, 'i' is the index of the point
    // }

    // --- POINT_DATA Section for Velocity ---
    fprintf(fp, "POINT_DATA %d\n", numPoints); // Point data section for the number of points

    // Add velocity data as an additional scalar field
    fprintf(fp, "SCALARS vel float 3\n"); // 3 components (x, y, z) for each velocity
    fprintf(fp, "LOOKUP_TABLE default\n"); // Lookup table name (default is fine)

    // Write the velocity data for each point (3 floats for each velocity)
    for (int i = 0; i < numPoints; i++) {
      fprintf(fp, "%.4e %.4e %.4e\n", velocity[i].x, velocity[i].y, velocity[i].z);
    }


    // Close the file
    fclose(fp);
    // printf("Successfully wrote point cloud to %s\n", filename);
}

void writePointCloudToVTKBinary(const char *filename, const vec *points, int numPoints) {
    FILE *fp = NULL;

    // Open the file for binary writing
    fp = fopen(filename, "wb");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }

    // --- VTK Legacy File Header ---
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Point Cloud\n"); // Title
    fprintf(fp, "BINARY\n");      // Data type (BINARY for faster processing)

    // --- Dataset Structure ---
    fprintf(fp, "DATASET POLYDATA\n"); // We're representing points, which can be thought of as a degenerate polydata

    // --- Points Section ---
    fprintf(fp, "POINTS %d float\n", numPoints); // Number of points and data type
    // Write the points in binary format
    // fwrite(points, sizeof(vec), numPoints, fp);
      // Write the points in binary format (3 floats for each point)
    for (int i = 0; i < numPoints; i++) {
        fwrite(&points[i].x, sizeof(float), 1, fp); // Write x-coordinate
        fwrite(&points[i].y, sizeof(float), 1, fp); // Write y-coordinate
        fwrite(&points[i].z, sizeof(float), 1, fp); // Write z-coordinate
    }

    // // --- Vertices Section (to make each point a distinct vertex) ---
    // // This tells VTK to render each point individually.
    // // The format is: num_vertices num_indices_total
    // // For each vertex: 1 index_of_point
    // fprintf(fp, "VERTICES %d %d\n", numPoints, numPoints * 2); // numPoints vertices, each with 1 index + its own count
    // // Write the vertex indices in binary format
    // for (int i = 0; i < numPoints; i++) {
    //     unsigned int index = i;
    //     fwrite(&index, sizeof(unsigned int), 1, fp);
    // }

    // Close the file
    fclose(fp);
    // printf("Successfully wrote point cloud to %s\n", filename);
}

#include "defs.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>


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
        fprintf(fp, "%lf %lf %lf\n", points[i].x, points[i].y, points[i].z);
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
      fprintf(fp, "%lf %lf %lf\n", velocity[i].x, velocity[i].y, velocity[i].z);
    }


    // Close the file
    fclose(fp);
    // printf("Successfully wrote point cloud to %s\n", filename);
}


// detect at runtime whether this machine is little-endian
static bool is_little_endian(void) {
    uint16_t one = 1;
    return *((uint8_t*)&one) == 1;
}

// byte-swap a 64-bit word
static uint64_t swap_uint64(uint64_t v) {
    return  ((v & 0x00000000000000FFULL) << 56) |
            ((v & 0x000000000000FF00ULL) << 40) |
            ((v & 0x0000000000FF0000ULL) << 24) |
            ((v & 0x00000000FF000000ULL) <<  8) |
            ((v & 0x000000FF00000000ULL) >>  8) |
            ((v & 0x0000FF0000000000ULL) >> 24) |
            ((v & 0x00FF000000000000ULL) >> 40) |
            ((v & 0xFF00000000000000ULL) >> 56);
}

// Byte-swap a 32-bit word
static uint32_t swap_uint32(uint32_t v) {
    return  ((v & 0x000000FF) << 24) |
            ((v & 0x0000FF00) << 8)  |
            ((v & 0x00FF0000) >> 8)  |
            ((v & 0xFF000000) >> 24);
}

// write one double in big-endian byte order
static void writeDoubleBE(FILE *fp, double val) {
    uint64_t bits;
    // copy the bit‐pattern of the double into a uint64_t
    memcpy(&bits, &val, sizeof(bits));
    // if host is little-endian, swap to network (big) order
    if (is_little_endian()) {
        bits = swap_uint64(bits);
    }
    // write the 8 bytes
    fwrite(&bits, sizeof(bits), 1, fp);
}

// Write one float in big-endian byte order
static void writeFloatBE(FILE *fp, float val) {
    uint32_t bits;
    // Copy the bit‐pattern of the float into a uint32_t
    memcpy(&bits, &val, sizeof(bits));
    // If host is little-endian, swap to network (big) order
    if (is_little_endian()) {
        bits = swap_uint32(bits);
    }
    // Write the 4 bytes
    fwrite(&bits, sizeof(bits), 1, fp);
}



void writePointCloudToVTKBinaryBigEndian(const char *filename,
                                         const vec *points,
                                         const vec *velocity,
                                         int numPoints)
{
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("Error opening file");
        return;
    }

    // --- VTK Legacy File Header ---
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Point Cloud\n");
    fprintf(fp, "BINARY\n");

    // point coordinates
    fprintf(fp, "DATASET POLYDATA\n");
    fprintf(fp, "POINTS %d float\n", numPoints);
    for (int i = 0; i < numPoints; i++) {
        writeFloatBE(fp, (float)points[i].x);
        writeFloatBE(fp, (float)points[i].y);
        writeFloatBE(fp, (float)points[i].z);
    }

    fprintf(fp, "POINT_DATA %d\n", numPoints); // Point data section for the number of points

    // velocity attribute
    fprintf(fp, "SCALARS vel float 3\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < numPoints; i++) {
        writeFloatBE(fp, (float)velocity[i].x);
        writeFloatBE(fp, (float)velocity[i].y);
        writeFloatBE(fp, (float)velocity[i].z);
    }

    fclose(fp);
}

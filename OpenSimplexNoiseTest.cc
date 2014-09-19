/*
 * OpenSimplex (Simplectic) Noise Test in C++
 * Copyright 2014 Arthur Tombs
 *
 * This file is intended to test the function of OpenSimplexNoise.hh.
 * It requires that you have development packages for libpng installed.
 * I am aware this is bad code; I wrote it in a hurry.
 *
 * Compile with:
 *   g++ -o OpenSimplexNoiseTest -O3 OpenSimplexNoiseTest.cc -lpng
 *
 */


#include <png.h>
#include <cmath>
#include <iostream>

#include "OpenSimplexNoise.hh"


const int WIDTH = 512;
const int HEIGHT = 512;
const double FEATURE_SIZE = 12.0;


int main (int argc, char **args) {

  OpenSimplexNoise noise;

  {
    FILE * fp;
    png_structp png_ptr;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;

    fp = fopen("noise.png", "wb");
    if (fp == NULL) {
      std::cerr << "Error: Failed to open file 'noise.png' for writing" << std::endl;
      goto finish;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
      std::cerr << "Error: Failed to allocate libpng write struct" << std::endl;
      goto finish;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
      std::cerr << "Error: Failed to allocate libpng info struct" << std::endl;
      goto finish;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      std::cerr << "Error: Failed to set libpng jump points" << std::endl;
      goto finish;
    }

    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, WIDTH, HEIGHT, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    for (int y = 0; y < HEIGHT; y++) {
      png_byte row [WIDTH * 3];
      for (int x = 0; x < WIDTH; x++) {
        double value = noise.eval((double)x / FEATURE_SIZE, (double)y / FEATURE_SIZE, 0.0);
        png_byte rgbval = (png_byte)std::floor((value * 0.5 + 0.5) * 255.0 + 0.5);
        row[x*3] = row[x*3+1] = row[x*3+2] = rgbval;
      }
      png_write_row(png_ptr, row);
    }

    png_write_end(png_ptr, NULL);

finish:

    if (fp != NULL) fclose(fp);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, NULL);

  }

  return 0;
}

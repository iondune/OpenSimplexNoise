/*
 * OpenSimplex (Simplectic) Noise Test in C++
 * by Arthur Tombs
 *
 * Modified 2014-12-04
 *
 * This file is intended to test the function of OpenSimplexNoise.hh
 * by creating example images of 2D, 3D and 4D noise.
 * It requires that you have development packages for libpng installed.
 *
 * Compile with:
 *   g++ -o OpenSimplexNoiseTest -O2 OpenSimplexNoiseTest.cc -lpng
 *
 * Additional optimization can be obtained with -Ofast (at the cost of accuracy)
 * and -msse4 (or the highest level of SSE your CPU supports).
 */


#include <png.h>
#include <cmath>
#include <iostream>

#include "OpenSimplexNoise.hh"


const int WIDTH = 512;
const int HEIGHT = 512;
const double FEATURE_SIZE = 24.0;

int write_image (const char * fname, double * data) {

  int retval = 0;

  FILE * fp;
  png_structp png_ptr;
  png_infop info_ptr = NULL;

  fp = fopen(fname, "wb");
  if (fp == NULL) {
    std::cerr << "Error: Failed to open file '" << fname << "' for writing" << std::endl;
    retval = 1;
    goto finish;
  }

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    std::cerr << "Error: Failed to allocate libpng write struct" << std::endl;
    retval = 2;
    goto finish;
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    std::cerr << "Error: Failed to allocate libpng info struct" << std::endl;
    retval = 3;
    goto finish;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    std::cerr << "Error: Failed to set libpng jump points" << std::endl;
    retval = 4;
    goto finish;
  }

  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, WIDTH, HEIGHT, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr, info_ptr);

  for (int y = 0; y < HEIGHT; y++) {
    png_byte row [WIDTH];
    for (int x = 0; x < WIDTH; x++) {
      png_byte rgbval = (png_byte)(std::min((data[x+y*WIDTH] + 1.0) * 128.0, 255.0));
      row[x] = rgbval;
    }
    png_write_row(png_ptr, row);
  }

  png_write_end(png_ptr, NULL);

finish:

  if (fp != NULL) fclose(fp);
  if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, NULL);

  return retval;

}

int main (void) {

  double * data = new double [WIDTH * HEIGHT];

  OSN::Noise<2, double> noise;

  for (int yi = 0; yi < HEIGHT; ++yi) {
    double y = yi / FEATURE_SIZE;
    for (int xi = 0; xi < WIDTH; ++xi) {
      double x = xi / FEATURE_SIZE;
      data[xi+yi*WIDTH] = noise.eval(x, y);
    }
  }

  write_image("noise2.png", data);

  return 0;
}

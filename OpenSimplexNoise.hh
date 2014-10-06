/*
 * OpenSimplex (Simplectic) Noise in C++
 * by Arthur Tombs
 *
 * Modified 2014-09-22
 *
 * This is a derivative work based on OpenSimplex by Kurt Spencer:
 *   https://gist.github.com/KdotJPG/b1270127455a94ac5d19
 *
 * Anyone is free to make use of this software in whatever way they want.
 * Attribution is appreciated, but not required.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef OPENSIMPLEXNOISE_HH
#define OPENSIMPLEXNOISE_HH

#include <cstdlib>
#include <cmath>

// (1 / sqrt(3 + 1) - 1) / 3
#define STRETCH_CONSTANT_3D (-1.0f / 6.0f)
// (sqrt(3 + 1) - 1) / 3
#define SQUISH_CONSTANT_3D (1.0f / 3.0f)

#define NORM_CONSTANT_3D 103.0f

class OpenSimplexNoise {

private:

  // Array of gradient values for 3D. Values are defined below the class definition.
  static const int gradients3D [72];

  // The default permutation order. Values are defined below the class definition.
  static const int permDefault [256];

  int perm [256];
  int permGradIndex3D [256];

  inline float extrapolate (long xsb, long ysb, long zsb, float dx, float dy, float dz) const {
    unsigned int index = permGradIndex3D[(perm[(perm[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF];
    return gradients3D[index] * dx +
           gradients3D[index + 1] * dy +
           gradients3D[index + 2] * dz;
  }

public:

  OpenSimplexNoise (const int * p = permDefault) {
    // Copy the supplied permutation array into this instance
    // TODO: Replace this with memcpy? Means including an extra header cstring
    for (int i = 0; i < 256; ++i) {
      perm[i] = p[i];
    }
    for (int i = 0; i < 256; ++i) {
      // NB: 72 is the number of elements in the gradients3D array
      permGradIndex3D[i] = (int)((perm[i] % (72 / 3)) * 3);
    }
  }

  // Initializes the class using a permutation array generated from a
  // 64-bit seed.
  // Generates a proper permutation (i.e. doesn't merely perform N
  // successive pair swaps on a base array)
  OpenSimplexNoise (int seed) {
    int source [256];
    for (int i = 0; i < 256; ++i) {
      source[i] = i;
    }
    // TODO: Use random number generator classes (C++11 only)
    srand(seed);
    for (int i = 255; i >= 0; --i) {
      int r = (rand() % (i+1));
      perm[i] = source[r];
      // NB: 72 is the number of elements in the gradients3D array
      permGradIndex3D[i] = (int)((perm[i] % (72 / 3)) * 3);
      source[r] = source[i];
    }
  }

  float eval (float x, float y, float z) const {

    // Place input coordinates on simplectic lattice.
    float stretchOffset = (x + y + z) * STRETCH_CONSTANT_3D;
    float xs = x + stretchOffset;
    float ys = y + stretchOffset;
    float zs = z + stretchOffset;

    // Floor to get simplectic lattice coordinates of rhombohedron
    // (stretched cube) super-cell.
    float xsbd = std::floor(xs);
    float ysbd = std::floor(ys);
    float zsbd = std::floor(zs);
    long xsb = (long)xsbd;
    long ysb = (long)ysbd;
    long zsb = (long)zsbd;

    // Skew out to get actual coordinates of rhombohedron origin.
    // These are needed later.
    float squishOffset = (xsbd + ysbd + zsbd) * SQUISH_CONSTANT_3D;
    float xb = xsbd + squishOffset;
    float yb = ysbd + squishOffset;
    float zb = zsbd + squishOffset;

    // Positions relative to origin point.
    float dx0 = x - xb;
    float dy0 = y - yb;
    float dz0 = z - zb;

    // Compute simplectic lattice coordinates relative to rhombohedral origin.
    float xins = xs - xsbd;
    float yins = ys - ysbd;
    float zins = zs - zsbd;

    // Sum together to get a value that determines which cell we are in.
    float inSum = xins + yins + zins;

    // These are given values inside the next block, and used afterwards.
    long xsv_ext0, ysv_ext0, zsv_ext0;
    long xsv_ext1, ysv_ext1, zsv_ext1;
    float dx_ext0, dy_ext0, dz_ext0;
    float dx_ext1, dy_ext1, dz_ext1;

    float value = 0.0f;

    if (inSum > 1.0f && inSum < 2.0f) {
      // The point is inside the octahedron (rectified 3-Simplex) inbetween.

      float aScore;
      unsigned char aPoint;
      bool aIsFurtherSide;
      float bScore;
      unsigned char bPoint;
      bool bIsFurtherSide;

      // Decide between point (1,0,0) and (0,1,1) as closest.
      float p1 = xins + yins;
      if (p1 <= 1.0f) {
        aScore = 1.0f - p1;
        aPoint = 4;
        aIsFurtherSide = false;
      } else {
        aScore = p1 - 1.0f;
        aPoint = 3;
        aIsFurtherSide = true;
      }

      // Decide between point (0,1,0) and (1,0,1) as closest.
      float p2 = xins + zins;
      if (p2 <= 1.0f) {
        bScore = 1.0f - p2;
        bPoint = 2;
        bIsFurtherSide = false;
      } else {
        bScore = p2 - 1.0f;
        bPoint = 5;
        bIsFurtherSide = true;
      }

      // The closest out of the two (0,0,1) and (1,1,0) will replace the
      // furthest out of the two decided above if closer.
      float p3 = yins + zins;
      if (p3 > 1.0f) {
        float score = p3 - 1.0f;
        if (aScore > bScore && bScore < score) {
          bScore = score;
          bPoint = 6;
          bIsFurtherSide = true;
        } else if (aScore <= bScore && aScore < score) {
          aScore = score;
          aPoint = 6;
          aIsFurtherSide = true;
        }
      } else {
        float score = 1.0f - p3;
        if (aScore > bScore && bScore < score) {
          bScore = score;
          bPoint = 1;
          bIsFurtherSide = false;
        } else if (aScore <= bScore && aScore < score) {
          aScore = score;
          aPoint = 1;
          aIsFurtherSide = false;
        }
      }

      // Where each of the two closest points are determines how the
      // extra two vertices are calculated.
      if (aIsFurtherSide == bIsFurtherSide) {
        if (aIsFurtherSide) {
          // Both closest points on (1,1,1) side.

          // One of the two extra points is (1,1,1)
          xsv_ext0 = xsb + 1;
          ysv_ext0 = ysb + 1;
          zsv_ext0 = zsb + 1;
          dx_ext0 = dx0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
          dy_ext0 = dy0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
          dz_ext0 = dz0 - 1.0f - (SQUISH_CONSTANT_3D * 3);

          // Other extra point is based on the shared axis.
          unsigned char c = aPoint & bPoint;
          if (c & 0x01) {
            xsv_ext1 = xsb + 2;
            ysv_ext1 = ysb;
            zsv_ext1 = zsb;
            dx_ext1 = dx0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
            dy_ext1 = dy0 - (SQUISH_CONSTANT_3D * 2);
            dz_ext1 = dz0 - (SQUISH_CONSTANT_3D * 2);
          } else if (c & 0x02) {
            xsv_ext1 = xsb;
            ysv_ext1 = ysb + 2;
            zsv_ext1 = zsb;
            dx_ext1 = dx0 - (SQUISH_CONSTANT_3D * 2);
            dy_ext1 = dy0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
            dz_ext1 = dz0 - (SQUISH_CONSTANT_3D * 2);
          } else {
            xsv_ext1 = xsb;
            ysv_ext1 = ysb;
            zsv_ext1 = zsb + 2;
            dx_ext1 = dx0 - (SQUISH_CONSTANT_3D * 2);
            dy_ext1 = dy0 - (SQUISH_CONSTANT_3D * 2);
            dz_ext1 = dz0 - 2 - (SQUISH_CONSTANT_3D * 2);
          }
        } else {
          // Both closest points are on the (0,0,0) side.

          // One of the two extra points is (0,0,0).
          xsv_ext0 = xsb;
          ysv_ext0 = ysb;
          zsv_ext0 = zsb;
          dx_ext0 = dx0;
          dy_ext0 = dy0;
          dz_ext0 = dz0;

          // The other extra point is based on the omitted axis.
          unsigned char c = aPoint | bPoint;
          if ((c & 0x01) == 0) {
            xsv_ext1 = xsb - 1;
            ysv_ext1 = ysb + 1;
            zsv_ext1 = zsb + 1;
            dx_ext1 = dx0 + 1.0f - SQUISH_CONSTANT_3D;
            dy_ext1 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
            dz_ext1 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
          } else if ((c & 0x02) == 0) {
            xsv_ext1 = xsb + 1;
            ysv_ext1 = ysb - 1;
            zsv_ext1 = zsb + 1;
            dx_ext1 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
            dy_ext1 = dy0 + 1.0f - SQUISH_CONSTANT_3D;
            dz_ext1 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
          } else {
            xsv_ext1 = xsb + 1;
            ysv_ext1 = ysb + 1;
            zsv_ext1 = zsb - 1;
            dx_ext1 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
            dy_ext1 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
            dz_ext1 = dz0 + 1.0f - SQUISH_CONSTANT_3D;
          }
        }
      } else {
        // TODO: Instrumentation suggests this branch is never taken in 2D
        // One point is on the (0,0,0) side, one point is on the (1,1,1) side.

        unsigned char c1, c2;
        if (aIsFurtherSide) {
          c1 = aPoint;
          c2 = bPoint;
        } else {
          c1 = bPoint;
          c2 = aPoint;
        }

        // One contribution is a permutation of (1,1,-1).
        if ((c1 & 0x01) == 0) {
          xsv_ext0 = xsb - 1;
          ysv_ext0 = ysb + 1;
          zsv_ext0 = zsb + 1;
          dx_ext0 = dx0 + 1.0f - SQUISH_CONSTANT_3D;
          dy_ext0 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
          dz_ext0 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
        } else if ((c1 & 0x02) == 0) {
          xsv_ext0 = xsb + 1;
          ysv_ext0 = ysb - 1;
          zsv_ext0 = zsb + 1;
          dx_ext0 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
          dy_ext0 = dy0 + 1.0f - SQUISH_CONSTANT_3D;
          dz_ext0 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
        } else {
          xsv_ext0 = xsb + 1;
          ysv_ext0 = ysb + 1;
          zsv_ext0 = zsb - 1;
          dx_ext0 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
          dy_ext0 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
          dz_ext0 = dz0 + 1.0f - SQUISH_CONSTANT_3D;
        }

        // One contribution is a permutation of (0,0,2).
        if (c2 & 0x01) {
          xsv_ext1 = xsb + 2;
          ysv_ext1 = ysb;
          zsv_ext1 = zsb;
          dx_ext1 = dx0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
          dy_ext1 = dy0 - (SQUISH_CONSTANT_3D * 2);
          dz_ext1 = dz0 - (SQUISH_CONSTANT_3D * 2);
        } else if (c2 & 0x02) {
          xsv_ext1 = xsb;
          ysv_ext1 = ysb + 2;
          zsv_ext1 = zsb;
          dx_ext1 = dx0 - (SQUISH_CONSTANT_3D * 2);
          dy_ext1 = dy0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
          dz_ext1 = dz0 - (SQUISH_CONSTANT_3D * 2);
        } else {
          xsv_ext1 = xsb;
          ysv_ext1 = ysb;
          zsv_ext1 = zsb + 2;
          dx_ext1 = dx0 - (SQUISH_CONSTANT_3D * 2);
          dy_ext1 = dy0 - (SQUISH_CONSTANT_3D * 2);
          dz_ext1 = dz0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
        }
      }

      // Contribution (0,0,1).
      float dx1 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
      float dy1 = dy0 - SQUISH_CONSTANT_3D;
      float dz1 = dz0 - SQUISH_CONSTANT_3D;
      float attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1);
      if (attn1 < 2.0f) {
        value = std::pow(2.0f-attn1, 4) * extrapolate(xsb + 1, ysb, zsb, dx1, dy1, dz1);
      }

      // Contribution (0,1,0).
      float dx2 = dx0 - SQUISH_CONSTANT_3D;
      float dy2 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
      float dz2 = dz1;
      float attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2);
      if (attn2 < 2.0f) {
        value += std::pow(2.0f-attn2, 4) * extrapolate(xsb, ysb + 1, zsb, dx2, dy2, dz2);
      }

      // Contribution (1,0,0).
      float dx3 = dx2;
      float dy3 = dy1;
      float dz3 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
      float attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3);
      if (attn3 < 2.0f) {
        value += std::pow(2.0f-attn3, 4) * extrapolate(xsb, ysb, zsb + 1, dx3, dy3, dz3);
      }

      // Contribution (1,1,0).
      float dx4 = dx0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
      float dy4 = dy0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
      float dz4 = dz0 - (SQUISH_CONSTANT_3D * 2);
      float attn4 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4);
      if (attn4 < 2.0f) {
        value += std::pow(2.0f-attn4, 4) * extrapolate(xsb + 1, ysb + 1, zsb + 0, dx4, dy4, dz4);
      }

      // Contribution (1,0,1).
      float dx5 = dx4;
      float dy5 = dy0 - (SQUISH_CONSTANT_3D * 2);
      float dz5 = dz0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
      float attn5 = (dx5 * dx5) + (dy5 * dy5) + (dz5 * dz5);
      if (attn5 < 2.0f) {
        value += std::pow(2.0f-attn5, 4) * extrapolate(xsb + 1, ysb, zsb + 1, dx5, dy5, dz5);
      }

      // Contribution (0,1,1).
      float dx6 = dx0 - (SQUISH_CONSTANT_3D * 2);
      float dy6 = dy4;
      float dz6 = dz5;
      float attn6 = (dx6 * dx6) + (dy6 * dy6) + (dz6 * dz6);
      if (attn6 < 2.0f) {
        value += std::pow(2.0f-attn6, 4) * extrapolate(xsb, ysb + 1, zsb + 1, dx6, dy6, dz6);
      }
    } else if (inSum <= 1.0f) {
      // The point is inside the tetrahedron (3-Simplex) at (0,0,0)

      // Determine which of (0,0,1), (0,1,0), (1,0,0) are closest.
      unsigned char aPoint = 1;
      float aScore = xins;
      unsigned char bPoint = 2;
      float bScore = yins;
      if (aScore < bScore && zins > aScore) {
        aScore = zins;
        aPoint = 4;
      } else if (aScore >= bScore && zins > bScore) {
        bScore = zins;
        bPoint = 4;
      }

      // Determine the two lattice points not part of the tetrahedron that may contribute.
      // This depends on the closest two tetrahedral vertices, including (0,0,0).
      float wins = 1.0f - inSum;
      if (wins > aScore || wins > bScore) {
        // (0,0,0) is one of the closest two tetrahedral vertices.

        // The other closest vertex is the closer of a and b.
        unsigned char c = ((bScore > aScore) ? bPoint : aPoint);

        if (c != 1) {
          xsv_ext0 = xsb - 1;
          xsv_ext1 = xsb;
          dx_ext0 = dx0 + 1.0f;
          dx_ext1 = dx0;
        } else {
          xsv_ext0 = xsv_ext1 = xsb + 1;
          dx_ext0 = dx_ext1 = dx0 - 1.0f;
        }

        if (c != 2) {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy_ext1 = dy0;
          if (c == 1) {
            ysv_ext0 -= 1;
            dy_ext0 += 1.0f;
          } else {
            ysv_ext1 -= 1;
            dy_ext1 += 1.0f;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy_ext1 = dy0 - 1.0f;
        }

        if (c != 4) {
          zsv_ext0 = zsb;
          zsv_ext1 = zsb - 1;
          dz_ext0 = dz0;
          dz_ext1 = dz0 + 1.0f;
        } else {
          zsv_ext0 = zsv_ext1 = zsb + 1;
          dz_ext0 = dz_ext1 = dz0 - 1.0f;
        }
      } else {
        // (0,0,0) is not one of the closest two tetrahedral vertices.

        // The two extra vertices are determined by the closest two.
        unsigned char c = (aPoint | bPoint);

        if (c & 0x01) {
          xsv_ext0 = xsv_ext1 = xsb + 1;
          dx_ext0 = dx0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
          dx_ext1 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
        } else {
          xsv_ext0 = xsb;
          xsv_ext1 = xsb - 1;
          dx_ext0 = dx0 - (SQUISH_CONSTANT_3D * 2);
          dx_ext1 = dx0 + 1.0f - SQUISH_CONSTANT_3D;
        }

        if (c & 0x02) {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
          dy_ext1 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
        } else {
          ysv_ext0 = ysb;
          ysv_ext1 = ysb - 1;
          dy_ext0 = dy0 - (SQUISH_CONSTANT_3D * 2);
          dy_ext1 = dy0 + 1.0f - SQUISH_CONSTANT_3D;
        }

        if (c & 0x04) {
          zsv_ext0 = zsv_ext1 = zsb + 1;
          dz_ext0 = dz0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
          dz_ext1 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
        } else {
          zsv_ext0 = zsb;
          zsv_ext1 = zsb - 1;
          dz_ext0 = dz0 - (SQUISH_CONSTANT_3D * 2);
          dz_ext1 = dz0 + 1.0f - SQUISH_CONSTANT_3D;
        }
      }

      // Contribution (0,0,0)
      float attn0 = (dx0 * dx0) + (dy0 * dy0) + (dz0 * dz0);
      if (attn0 < 2.0f) {
        value = std::pow(2.0f-attn0, 4) * extrapolate(xsb, ysb, zsb, dx0, dy0, dz0);
      }

      // Contribution (0,0,1)
      float dx1 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
      float dy1 = dy0 - SQUISH_CONSTANT_3D;
      float dz1 = dz0 - SQUISH_CONSTANT_3D;
      float attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1);
      if (attn1 < 2.0f) {
        value += std::pow(2.0f-attn1, 4) * extrapolate(xsb + 1, ysb, zsb, dx1, dy1, dz1);
      }

      // Contribution (0,1,0)
      float dx2 = dx0 - SQUISH_CONSTANT_3D;
      float dy2 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
      float dz2 = dz1;
      float attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2);
      if (attn2 < 2.0f) {
        value += std::pow(2.0f-attn2, 4) * extrapolate(xsb, ysb + 1, zsb, dx2, dy2, dz2);
      }

      // Contribution (1,0,0)
      float dx3 = dx2;
      float dy3 = dy1;
      float dz3 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
      float attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3);
      if (attn3 < 2.0f) {
        value += std::pow(2.0f-attn3, 4) * extrapolate(xsb, ysb, zsb + 1, dx3, dy3, dz3);
      }
    } else {
      // The point is inside the tetrahedron (3-Simplex) at (1,1,1)

      // Determine which two tetrahedral vertices are the closest
      // out of (1,1,0), (1,0,1), and (0,1,1), but not (1,1,1).
      unsigned char aPoint = 6;
      float aScore = xins;
      unsigned char bPoint = 5;
      float bScore = yins;
      if (aScore <= bScore && zins < bScore) {
        bScore = zins;
        bPoint = 3;
      } else if (aScore > bScore && zins < aScore) {
        aScore = zins;
        aPoint = 3;
      }

      // Determine the two lattice points not part of the tetrahedron that may contribute.
      // This depends on the closest two tetrahedral vertices, including (1,1,1).
      float wins = 3.0f - inSum;
      if (wins < aScore || wins < bScore) {
        // (1,1,1) is one of the closest two tetrahedral vertices.

        // The other closest vertex is the closest of a and b.
        unsigned char c = ((bScore < aScore) ? bPoint : aPoint);

        if (c & 0x01) {
          xsv_ext0 = xsb + 2;
          xsv_ext1 = xsb + 1;
          dx_ext0 = dx0 - 2.0f - (SQUISH_CONSTANT_3D * 3);
          dx_ext1 = dx0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
        } else {
          xsv_ext0 = xsv_ext1 = xsb;
          dx_ext0 = dx_ext1 = dx0 - (SQUISH_CONSTANT_3D * 3);
        }

        if (c & 0x02) {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy_ext1 = dy0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
          if (c & 0x01) {
            ysv_ext1 += 1;
            dy_ext1 -= 1.0f;
          } else {
            ysv_ext0 += 1;
            dy_ext0 -= 1.0f;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy_ext1 = dy0 - (SQUISH_CONSTANT_3D * 3);
        }

        if (c & 0x04) {
          zsv_ext0 = zsb + 1;
          zsv_ext1 = zsb + 2;
          dz_ext0 = dz0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
          dz_ext1 = dz0 - 2.0f - (SQUISH_CONSTANT_3D * 3);
        } else {
          zsv_ext0 = zsv_ext1 = zsb;
          dz_ext0 = dz_ext1 = dz0 - (SQUISH_CONSTANT_3D * 3);
        }
      } else {
        // (1,1,1) is not one of the closest two tetrahedral vertices.

        // The two extra vertices are determined by the closest two.
        unsigned char c = aPoint & bPoint;

        if (c & 0x01) {
          xsv_ext0 = xsb + 1;
          xsv_ext1 = xsb + 2;
          dx_ext0 = dx0 - 1.0f - SQUISH_CONSTANT_3D;
          dx_ext1 = dx0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
        } else {
          xsv_ext0 = xsv_ext1 = xsb;
          dx_ext0 = dx0 - SQUISH_CONSTANT_3D;
          dx_ext1 = dx0 - (SQUISH_CONSTANT_3D * 2);
        }

        if (c & 0x02) {
          ysv_ext0 = ysb + 1;
          ysv_ext1 = ysb + 2;
          dy_ext0 = dy0 - 1.0f - SQUISH_CONSTANT_3D;
          dy_ext1 = dy0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
        } else {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy0 - SQUISH_CONSTANT_3D;
          dy_ext1 = dy0 - (SQUISH_CONSTANT_3D * 2);
        }

        if (c & 0x04) {
          zsv_ext0 = zsb + 1;
          zsv_ext1 = zsb + 2;
          dz_ext0 = dz0 - 1.0f - SQUISH_CONSTANT_3D;
          dz_ext1 = dz0 - 2.0f - (SQUISH_CONSTANT_3D * 2);
        } else {
          zsv_ext0 = zsv_ext1 = zsb;
          dz_ext0 = dz0 - SQUISH_CONSTANT_3D;
          dz_ext1 = dz0 - (SQUISH_CONSTANT_3D * 2);
        }
      }

      // Contribution (1,1,0)
      float dx3 = dx0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
      float dy3 = dy0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
      float dz3 = dz0 - (SQUISH_CONSTANT_3D * 2);
      float attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3);
      if (attn3 < 2.0f) {
        value = std::pow(2.0f-attn3, 4) * extrapolate(xsb + 1, ysb + 1, zsb, dx3, dy3, dz3);
      }

      // Contribution (1,0,1)
      float dx2 = dx3;
      float dy2 = dy0 - (SQUISH_CONSTANT_3D * 2);
      float dz2 = dz0 - 1.0f - (SQUISH_CONSTANT_3D * 2);
      float attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2);
      if (attn2 < 2.0f) {
        value += std::pow(2.0f-attn2, 4) * extrapolate(xsb + 1, ysb, zsb + 1, dx2, dy2, dz2);
      }

      // Contribution (0,1,1)
      float dx1 = dx0 - (SQUISH_CONSTANT_3D * 2);
      float dy1 = dy3;
      float dz1 = dz2;
      float attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1);
      if (attn1 < 2.0f) {
        value += std::pow(2.0f-attn1, 4) * extrapolate(xsb, ysb + 1, zsb + 1, dx1, dy1, dz1);
      }

      // Contribution (1,1,1)
      dx0 = dx0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
      dy0 = dy0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
      dz0 = dz0 - 1.0f - (SQUISH_CONSTANT_3D * 3);
      float attn0 = (dx0 * dx0) + (dy0 * dy0) + (dz0 * dz0);
      if (attn0 < 2.0f) {
        value += std::pow(2.0f-attn0, 4) * extrapolate(xsb + 1, ysb + 1, zsb + 1, dx0, dy0, dz0);
      }
    }

    // First extra vertex.
    float attn_ext0 = (dx_ext0 * dx_ext0) + (dy_ext0 * dy_ext0) + (dz_ext0 * dz_ext0);
    if (attn_ext0 < 2.0f) {
      value += std::pow(2.0f-attn_ext0, 4) * extrapolate(xsv_ext0, ysv_ext0, zsv_ext0, dx_ext0, dy_ext0, dz_ext0);
    }

    // Second extra vertex.
    float attn_ext1 = (dx_ext1 * dx_ext1) + (dy_ext1 * dy_ext1) + (dz_ext1 * dz_ext1);
    if (attn_ext1 < 2.0f) {
      value += std::pow(2.0f-attn_ext1, 4) * extrapolate(xsv_ext1, ysv_ext1, zsv_ext1, dx_ext1, dy_ext1, dz_ext1);
    }

    return (value / NORM_CONSTANT_3D);
  }

};

// Array of gradient values for 3D. They approximate the directions to the
// vertices of a rhombicuboctahedron from its center, skewed so that the
// triangular and square facets can be inscribed in circles of the same radius.
// New gradient set 2014-10-06.
const int OpenSimplexNoise::gradients3D [72] = {
  -11, 4, 4,  -4, 11, 4,  -4, 4, 11,   11, 4, 4,   4, 11, 4,   4, 4, 11,
  -11,-4, 4,  -4,-11, 4,  -4,-4, 11,   11,-4, 4,   4,-11, 4,   4,-4, 11,
  -11, 4,-4,  -4, 11,-4,  -4, 4,-11,   11, 4,-4,   4, 11,-4,   4, 4,-11,
  -11,-4,-4,  -4,-11,-4,  -4,-4,-11,   11,-4,-4,   4,-11,-4,   4,-4,-11
};

// The standard permutation order as used in Ken Perlin's "Improved Noise"
// (and basically every noise implementation on the Internet).
const int OpenSimplexNoise::permDefault [256] = {
  151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
  140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
  247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32,
   57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
   74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
   60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
   65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
  200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
   52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
  207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
  119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
  129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
  218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
   81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
  184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
  222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180
};


#undef STRETCH_CONSTANT_3D
#undef SQUISH_CONSTANT_3D
#undef NORM_CONSTANT


#endif

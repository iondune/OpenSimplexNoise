/*
 * OpenSimplex (Simplectic) Noise in C++
 * by Arthur Tombs
 *
 * Modified 2014-12-04
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

#include <cmath>

#if __cplusplus < 201103L
#pragma message("Info: Your compiler does not claim C++11 support. Some features may be unavailable.")
#else
#define OSN_USE_LONG_LONG
#define OSN_USE_STATIC_ASSERT
#endif

#ifndef OSN_USE_LONG_LONG
#pragma message("Info: Using legacy constructors. Define OSN_USE_LONG_LONG before including this header to force use of the 'long long' type.")
// cstdlib is required for the srand and rand functions
#include <cstdlib>
#endif

#ifdef OSN_USE_STATIC_ASSERT
#include <type_traits>
#endif

namespace OSN {

typedef unsigned char OSN_BYTE;


class NoiseBase {

protected:

  int perm [256];

  // Empty constructor to allow child classes to set up perm themselves.
  NoiseBase (void) {}

#ifdef OSN_USE_LONG_LONG
  // Perform one step of the Linear Congruential Generator algorithm.
  inline static void LCG_STEP (long long & x) {
    // Magic constants are attributed to Donald Knuth's MMIX implementation.
    const long long MULTIPLIER = 6364136223846793005LL;
    const long long INCREMENT  = 1442695040888963407LL;
    x = ((x * MULTIPLIER) + INCREMENT);
  }

  // Initializes the class using a permutation array generated from a 64-bit seed.
  // Generates a proper permutation (i.e. doesn't merely perform N successive
  // pair swaps on a base array).
  // Uses a simple 64-bit LCG.
  NoiseBase (long long seed) {
    int source [256];
    for (int i = 0; i < 256; ++i) {
      source[i] = i;
    }
    LCG_STEP(seed);
    LCG_STEP(seed);
    LCG_STEP(seed);
    for (int i = 255; i >= 0; --i) {
      LCG_STEP(seed);
      int r = (int)((seed + 31) % (i + 1));
      if (r < 0) r += (i + 1);
      perm[i] = source[r];
      source[r] = source[i];
    }
  }
#else
  // Initializes the class using a permutation array generated from a 32-bit seed.
  // Generates a proper permutation (i.e. doesn't merely perform N successive
  // pair swaps on a base array).
  NoiseBase (long seed) {
    int source [256];
    for (int i = 0; i < 256; ++i) {
      source[i] = i;
    }
    srand(seed);
    for (int i = 255; i >= 0; --i) {
      int r = (int)(rand() % (i + 1));
      perm[i] = source[r];
      source[r] = source[i];
    }
  }
#endif

  NoiseBase (const int * p) {
    // Copy the supplied permutation array into this instance
    for (int i = 0; i < 256; ++i) {
      perm[i] = p[i];
    }
  }

};


template <int N, typename T = double>
class Noise : public NoiseBase {
};

// 2D Implementation of the OpenSimplexNoise generator.
template <typename T>
class Noise <2, T> : public NoiseBase {
private:

  static const int gradients [];

  inline T extrapolate (long xsb, long ysb, T dx, T dy) const {
    unsigned int index = perm[(perm[xsb & 0xFF] + ysb) & 0xFF] & 0x0E;
    return gradients[index] * dx +
           gradients[index + 1] * dy;
  }

public:

#ifdef OSN_USE_LONG_LONG
  Noise (long long seed = 0LL) : NoiseBase (seed) {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif
  }
#else
  Noise (long seed = 0L) : NoiseBase (seed) {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif
  }
#endif
  Noise (const int * p) : NoiseBase (p) {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif
  }

  T eval (T x, T y) const {

    const T STRETCH_CONSTANT = (T)((1.0 / std::sqrt(2.0 + 1.0) - 1.0) * 0.5);
    const T SQUISH_CONSTANT  = (T)((std::sqrt(2.0 + 1.0) - 1.0) * 0.5);
    const T NORM_CONSTANT    = (T)47.0;

    long xsb, ysb, xsv_ext, ysv_ext;
    T dx0, dy0, dx_ext, dy_ext;
    T xins, yins;

    // Parameters for the four contributions
    T contr_m [4], contr_ext [4];

    {
      // Place input coordinates on a grid.
      T stretchOffset = (x + y) * STRETCH_CONSTANT;
      T xs = x + stretchOffset;
      T ys = y + stretchOffset;

      // Floor to get grid coordinates of rhombus super-cell origin.
      T xsbd = std::floor(xs);
      T ysbd = std::floor(ys);
      xsb = (long)xsbd;
      ysb = (long)ysbd;

      // Skew out to get actual coordinates of rhombohedron origin.
      T squishOffset = (xsbd + ysbd) * SQUISH_CONSTANT;
      T xb = xsbd + squishOffset;
      T yb = ysbd + squishOffset;

      // Positions relative to origin point.
      dx0 = x - xb;
      dy0 = y - yb;

      // Compute grid coordinates relative to rhomboidal origin.
      xins = xs - xsbd;
      yins = ys - ysbd;
    }

    // Contribution (1,0).
    {
      T dx1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
      T dy1 = dy0 - SQUISH_CONSTANT;
      T attn1 = (dx1 * dx1) + (dy1 * dy1);
      contr_m[0]   = attn1;
      contr_ext[0] = extrapolate(xsb + 1, ysb, dx1, dy1);
    }

    // Contribution (0,1).
    {
      T dx2 = dx0 - SQUISH_CONSTANT;
      T dy2 = dy0 - (T)1.0 - SQUISH_CONSTANT;
      T attn2 = (dx2 * dx2) + (dy2 * dy2);
      contr_m[1]   = attn2;
      contr_ext[1] = extrapolate(xsb, ysb + 1, dx2, dy2);
    }

    if ((xins + yins) <= (T)1.0) {
      // Inside the triangle (2-Simplex) at (0,0).
      T zins = (T)1.0 - (xins + yins);
      if (zins > xins || zins > yins) {
        // (0,0) is one of the closest two triangular vertices.
        if (xins > yins) {
          xsv_ext = xsb + 1;
          ysv_ext = ysb - 1;
          dx_ext = dx0 - (T)1.0;
          dy_ext = dy0 + (T)1.0;
        } else {
          xsv_ext = xsb - 1;
          ysv_ext = ysb + 1;
          dx_ext = dx0 + (T)1.0;
          dy_ext = dy0 - (T)1.0;
        }
      } else {
        // (1,0) and (0,1) are the closest two vertices.
        xsv_ext = xsb + 1;
        ysv_ext = ysb + 1;
        dx_ext = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        dy_ext = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      }
    } else {
      // Inside the triangle (2-Simplex) at (1,1).
      T zins = (T)2.0 - (xins + yins);
      if (zins < xins || zins < yins) {
        // (0,0) is one of the closest two triangular vertices.
        if (xins > yins) {
          xsv_ext = xsb + 2;
          ysv_ext = ysb;
          dx_ext = dx0 - (T)2.0 - (SQUISH_CONSTANT * 2);
          dy_ext = dy0 - (SQUISH_CONSTANT * 2);
        } else {
          xsv_ext = xsb;
          ysv_ext = ysb + 2;
          dx_ext = dx0 - (SQUISH_CONSTANT * 2);
          dy_ext = dy0 - (T)2.0 - (SQUISH_CONSTANT * 2);
        }
      } else {
        // (1,0) and (0,1) are the closest two vertices.
        xsv_ext = xsb;
        ysv_ext = ysb;
        dx_ext = dx0;
        dy_ext = dy0;
      }
      xsb += 1;
      ysb += 1;
      dx0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      dy0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
    }

    // Contribution (0,0) or (1,1).
    {
      T attn0 = (dx0 * dx0) + (dy0 * dy0);
      contr_m[2]   = attn0;
      contr_ext[2] = extrapolate(xsb, ysb, dx0, dy0);
    }

    // Extra vertex.
    {
      T attn_ext = (dx_ext * dx_ext) + (dy_ext * dy_ext);
      contr_m[3]   = attn_ext;
      contr_ext[3] = extrapolate(xsv_ext, ysv_ext, dx_ext, dy_ext);
    }

    T value = 0.0;
    for (int i=0; i<4; ++i) {
      value += std::pow(std::max((T)2.0 - contr_m[i], (T)0.0), 4) * contr_ext[i];
    }

    return (value / NORM_CONSTANT);
  }

};

// Array of gradient values for 2D. They approximate the directions to the
// vertices of a octagon from its center.
// Gradient set 2014-10-06.
template <typename T>
const int Noise<2, T>::gradients [] = {
   5, 2,   2, 5,  -5, 2,  -2, 5,
   5,-2,   2,-5,  -5,-2,  -2,-5
};


// 3D Implementation of the OpenSimplexNoise generator.
template <typename T>
class Noise <3, T> : public NoiseBase {
private:

  // Array of gradient values for 3D. Values are defined below the class definition.
  static const int gradients [72];

  // Because 72 is not a power of two, extrapolate cannot use a bitmask to index
  // into the perm array. Pre-calculate and store the indices instead.
  int permGradIndex [256];

  inline T extrapolate (long xsb, long ysb, long zsb, T dx, T dy, T dz) const {
    unsigned int index = permGradIndex[(perm[(perm[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF];
    return gradients[index] * dx +
           gradients[index + 1] * dy +
           gradients[index + 2] * dz;
  }

public:

#ifdef OSN_USE_LONG_LONG
  // Initializes the class using a permutation array generated from a 64-bit seed.
  // Generates a proper permutation (i.e. doesn't merely perform N successive
  // pair swaps on a base array).
  // Uses a simple 64-bit LCG.
  Noise (long long seed = 0LL) : NoiseBase () {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif

    int source [256];
    for (int i = 0; i < 256; ++i) {
      source[i] = i;
    }
    LCG_STEP(seed);
    LCG_STEP(seed);
    LCG_STEP(seed);
    for (int i = 255; i >= 0; --i) {
      LCG_STEP(seed);
      int r = (int)((seed + 31) % (i + 1));
      if (r < 0) r += (i + 1);
      perm[i] = source[r];
      permGradIndex[i] = (int)((perm[i] % (72 / 3)) * 3);
      source[r] = source[i];
    }
  }
#else
  // Initializes the class using a permutation array generated from a 32-bit seed.
  // Generates a proper permutation (i.e. doesn't merely perform N successive
  // pair swaps on a base array).
  Noise (long seed = 0L) : NoiseBase () {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif

    int source [256];
    for (int i = 0; i < 256; ++i) {
      source[i] = i;
    }
    srand(seed);
    for (int i = 255; i >= 0; --i) {
      int r = (int)(rand() % (i + 1));
      perm[i] = source[r];
      // NB: 72 is the number of elements of the gradients3D array
      permGradIndex[i] = (int)((perm[i] % (72 / 3)) * 3);
      source[r] = source[i];
    }
  }
#endif

  Noise (const int * p) : NoiseBase () {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif

    // Copy the supplied permutation array into this instance.
    for (int i = 0; i < 256; ++i) {
      perm[i] = p[i];
      permGradIndex[i] = (int)((perm[i] % (72 / 3)) * 3);
    }
  }


  T eval (T x, T y, T z) const {

    const T STRETCH_CONSTANT = (T)(-1.0 / 6.0); // (1 / sqrt(3 + 1) - 1) / 3
    const T SQUISH_CONSTANT  = (T)(1.0 / 3.0);  // (sqrt(3 + 1) - 1) / 3
    const T NORM_CONSTANT    = (T)103.0;

    long xsb, ysb, zsb;
    T dx0, dy0, dz0;
    T xins, yins, zins;

    {
      // Place input coordinates on simplectic lattice.
      T stretchOffset = (x + y + z) * STRETCH_CONSTANT;
      T xs = x + stretchOffset;
      T ys = y + stretchOffset;
      T zs = z + stretchOffset;

      // Floor to get simplectic lattice coordinates of rhombohedron
      // (stretched cube) super-cell.
      T xsbd = std::floor(xs);
      T ysbd = std::floor(ys);
      T zsbd = std::floor(zs);
      xsb = (long)xsbd;
      ysb = (long)ysbd;
      zsb = (long)zsbd;

      // Skew out to get actual coordinates of rhombohedron origin.
      T squishOffset = (xsbd + ysbd + zsbd) * SQUISH_CONSTANT;
      T xb = xsbd + squishOffset;
      T yb = ysbd + squishOffset;
      T zb = zsbd + squishOffset;

      // Positions relative to origin point.
      dx0 = x - xb;
      dy0 = y - yb;
      dz0 = z - zb;

      // Compute simplectic lattice coordinates relative to rhombohedral origin.
      xins = xs - xsbd;
      yins = ys - ysbd;
      zins = zs - zsbd;
    }

    // These are given values inside the next block, and used afterwards.
    long xsv_ext0, ysv_ext0, zsv_ext0;
    long xsv_ext1, ysv_ext1, zsv_ext1;
    T dx_ext0, dy_ext0, dz_ext0;
    T dx_ext1, dy_ext1, dz_ext1;

    T value = 0.0;

    // Sum together to get a value that determines which cell we are in.
    T inSum = xins + yins + zins;

    if (inSum > (T)1.0 && inSum < (T)2.0) {
      // The point is inside the octahedron (rectified 3-Simplex) inbetween.

      T aScore;
      OSN_BYTE aPoint;
      bool aIsFurtherSide;
      T bScore;
      OSN_BYTE bPoint;
      bool bIsFurtherSide;

      // Decide between point (1,0,0) and (0,1,1) as closest.
      T p1 = xins + yins;
      if (p1 <= (T)1.0) {
        aScore = (T)1.0 - p1;
        aPoint = 4;
        aIsFurtherSide = false;
      } else {
        aScore = p1 - (T)1.0;
        aPoint = 3;
        aIsFurtherSide = true;
      }

      // Decide between point (0,1,0) and (1,0,1) as closest.
      T p2 = xins + zins;
      if (p2 <= (T)1.0) {
        bScore = (T)1.0 - p2;
        bPoint = 2;
        bIsFurtherSide = false;
      } else {
        bScore = p2 - (T)1.0;
        bPoint = 5;
        bIsFurtherSide = true;
      }

      // The closest out of the two (0,0,1) and (1,1,0) will replace the
      // furthest out of the two decided above if closer.
      T p3 = yins + zins;
      if (p3 > (T)1.0) {
        T score = p3 - (T)1.0;
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
        T score = (T)1.0 - p3;
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
          dx_ext0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          dy_ext0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          dz_ext0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);

          // Other extra point is based on the shared axis.
          OSN_BYTE c = aPoint & bPoint;
          if (c & 0x01) {
            xsv_ext1 = xsb + 2;
            ysv_ext1 = ysb;
            zsv_ext1 = zsb;
            dx_ext1 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 2);
            dy_ext1 = dy0 - (SQUISH_CONSTANT * 2);
            dz_ext1 = dz0 - (SQUISH_CONSTANT * 2);
          } else if (c & 0x02) {
            xsv_ext1 = xsb;
            ysv_ext1 = ysb + 2;
            zsv_ext1 = zsb;
            dx_ext1 = dx0 - (SQUISH_CONSTANT * 2);
            dy_ext1 = dy0 - (T)2.0 - (SQUISH_CONSTANT * 2);
            dz_ext1 = dz0 - (SQUISH_CONSTANT * 2);
          } else {
            xsv_ext1 = xsb;
            ysv_ext1 = ysb;
            zsv_ext1 = zsb + 2;
            dx_ext1 = dx0 - (SQUISH_CONSTANT * 2);
            dy_ext1 = dy0 - (SQUISH_CONSTANT * 2);
            dz_ext1 = dz0 - (T)2.0 - (SQUISH_CONSTANT * 2);
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
          OSN_BYTE c = aPoint | bPoint;
          if ((c & 0x01) == 0) {
            xsv_ext1 = xsb - 1;
            ysv_ext1 = ysb + 1;
            zsv_ext1 = zsb + 1;
            dx_ext1 = dx0 + (T)1.0 - SQUISH_CONSTANT;
            dy_ext1 = dy0 - (T)1.0 - SQUISH_CONSTANT;
            dz_ext1 = dz0 - (T)1.0 - SQUISH_CONSTANT;
          } else if ((c & 0x02) == 0) {
            xsv_ext1 = xsb + 1;
            ysv_ext1 = ysb - 1;
            zsv_ext1 = zsb + 1;
            dx_ext1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
            dy_ext1 = dy0 + (T)1.0 - SQUISH_CONSTANT;
            dz_ext1 = dz0 - (T)1.0 - SQUISH_CONSTANT;
          } else {
            xsv_ext1 = xsb + 1;
            ysv_ext1 = ysb + 1;
            zsv_ext1 = zsb - 1;
            dx_ext1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
            dy_ext1 = dy0 - (T)1.0 - SQUISH_CONSTANT;
            dz_ext1 = dz0 + (T)1.0 - SQUISH_CONSTANT;
          }
        }
      } else {
        // One point is on the (0,0,0) side, one point is on the (1,1,1) side.

        OSN_BYTE c1, c2;
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
          dx_ext0 = dx0 + (T)1.0 - SQUISH_CONSTANT;
          dy_ext0 = dy0 - (T)1.0 - SQUISH_CONSTANT;
          dz_ext0 = dz0 - (T)1.0 - SQUISH_CONSTANT;
        } else if ((c1 & 0x02) == 0) {
          xsv_ext0 = xsb + 1;
          ysv_ext0 = ysb - 1;
          zsv_ext0 = zsb + 1;
          dx_ext0 = dx0 - (T)1.0 - SQUISH_CONSTANT;
          dy_ext0 = dy0 + (T)1.0 - SQUISH_CONSTANT;
          dz_ext0 = dz0 - (T)1.0 - SQUISH_CONSTANT;
        } else {
          xsv_ext0 = xsb + 1;
          ysv_ext0 = ysb + 1;
          zsv_ext0 = zsb - 1;
          dx_ext0 = dx0 - (T)1.0 - SQUISH_CONSTANT;
          dy_ext0 = dy0 - (T)1.0 - SQUISH_CONSTANT;
          dz_ext0 = dz0 + (T)1.0 - SQUISH_CONSTANT;
        }

        // One contribution is a permutation of (0,0,2).
        if (c2 & 0x01) {
          xsv_ext1 = xsb + 2;
          ysv_ext1 = ysb;
          zsv_ext1 = zsb;
          dx_ext1 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz0 - (SQUISH_CONSTANT * 2);
        } else if (c2 & 0x02) {
          xsv_ext1 = xsb;
          ysv_ext1 = ysb + 2;
          zsv_ext1 = zsb;
          dx_ext1 = dx0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy0 - (T)2.0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz0 - (SQUISH_CONSTANT * 2);
        } else {
          xsv_ext1 = xsb;
          ysv_ext1 = ysb;
          zsv_ext1 = zsb + 2;
          dx_ext1 = dx0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz0 - (T)2.0 - (SQUISH_CONSTANT * 2);
        }
      }

      // Contribution (0,0,1).
      T dx1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
      T dy1 = dy0 - SQUISH_CONSTANT;
      T dz1 = dz0 - SQUISH_CONSTANT;
      T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1);
      value = std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb, dx1, dy1, dz1);

      // Contribution (0,1,0).
      T dx2 = dx0 - SQUISH_CONSTANT;
      T dy2 = dy0 - (T)1.0 - SQUISH_CONSTANT;
      T dz2 = dz1;
      T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2);
      value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb, dx2, dy2, dz2);

      // Contribution (1,0,0).
      T dx3 = dx2;
      T dy3 = dy1;
      T dz3 = dz0 - (T)1.0 - SQUISH_CONSTANT;
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3);
      value += std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb, ysb, zsb + 1, dx3, dy3, dz3);

      // Contribution (1,1,0).
      T dx4 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dy4 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dz4 = dz0 - (SQUISH_CONSTANT * 2);
      T attn4 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4);
      value += std::pow(std::max((T)2.0 - attn4, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb + 0, dx4, dy4, dz4);

      // Contribution (1,0,1).
      T dx5 = dx4;
      T dy5 = dy0 - (SQUISH_CONSTANT * 2);
      T dz5 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T attn5 = (dx5 * dx5) + (dy5 * dy5) + (dz5 * dz5);
      value += std::pow(std::max((T)2.0 - attn5, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb + 1, dx5, dy5, dz5);

      // Contribution (0,1,1).
      T dx6 = dx0 - (SQUISH_CONSTANT * 2);
      T dy6 = dy4;
      T dz6 = dz5;
      T attn6 = (dx6 * dx6) + (dy6 * dy6) + (dz6 * dz6);
      value += std::pow(std::max((T)2.0 - attn6, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb + 1, dx6, dy6, dz6);

    } else if (inSum <= (T)1.0) {
      // The point is inside the tetrahedron (3-Simplex) at (0,0,0)

      // Determine which of (0,0,1), (0,1,0), (1,0,0) are closest.
      OSN_BYTE aPoint = 1;
      T aScore = xins;
      OSN_BYTE bPoint = 2;
      T bScore = yins;
      if (aScore < bScore && zins > aScore) {
        aScore = zins;
        aPoint = 4;
      } else if (aScore >= bScore && zins > bScore) {
        bScore = zins;
        bPoint = 4;
      }

      // Determine the two lattice points not part of the tetrahedron that may contribute.
      // This depends on the closest two tetrahedral vertices, including (0,0,0).
      T wins = (T)1.0 - inSum;
      if (wins > aScore || wins > bScore) {
        // (0,0,0) is one of the closest two tetrahedral vertices.

        // The other closest vertex is the closer of a and b.
        OSN_BYTE c = ((bScore > aScore) ? bPoint : aPoint);

        if (c != 1) {
          xsv_ext0 = xsb - 1;
          xsv_ext1 = xsb;
          dx_ext0 = dx0 + (T)1.0;
          dx_ext1 = dx0;
        } else {
          xsv_ext0 = xsv_ext1 = xsb + 1;
          dx_ext0 = dx_ext1 = dx0 - (T)1.0;
        }

        if (c != 2) {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy_ext1 = dy0;
          if (c == 1) {
            ysv_ext0 -= 1;
            dy_ext0 += (T)1.0;
          } else {
            ysv_ext1 -= 1;
            dy_ext1 += (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy_ext1 = dy0 - (T)1.0;
        }

        if (c != 4) {
          zsv_ext0 = zsb;
          zsv_ext1 = zsb - 1;
          dz_ext0 = dz0;
          dz_ext1 = dz0 + (T)1.0;
        } else {
          zsv_ext0 = zsv_ext1 = zsb + 1;
          dz_ext0 = dz_ext1 = dz0 - (T)1.0;
        }
      } else {
        // (0,0,0) is not one of the closest two tetrahedral vertices.

        // The two extra vertices are determined by the closest two.
        OSN_BYTE c = (aPoint | bPoint);

        if (c & 0x01) {
          xsv_ext0 = xsv_ext1 = xsb + 1;
          dx_ext0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dx_ext1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
        } else {
          xsv_ext0 = xsb;
          xsv_ext1 = xsb - 1;
          dx_ext0 = dx0 - (SQUISH_CONSTANT * 2);
          dx_ext1 = dx0 + (T)1.0 - SQUISH_CONSTANT;
        }

        if (c & 0x02) {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy0 - (T)1.0 - SQUISH_CONSTANT;
        } else {
          ysv_ext0 = ysb;
          ysv_ext1 = ysb - 1;
          dy_ext0 = dy0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy0 + (T)1.0 - SQUISH_CONSTANT;
        }

        if (c & 0x04) {
          zsv_ext0 = zsv_ext1 = zsb + 1;
          dz_ext0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz0 - (T)1.0 - SQUISH_CONSTANT;
        } else {
          zsv_ext0 = zsb;
          zsv_ext1 = zsb - 1;
          dz_ext0 = dz0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz0 + (T)1.0 - SQUISH_CONSTANT;
        }
      }

      // Contribution (0,0,0)
      T attn0 = (dx0 * dx0) + (dy0 * dy0) + (dz0 * dz0);
      value = std::pow(std::max((T)2.0 - attn0, (T)0.0), 4) * extrapolate(xsb, ysb, zsb, dx0, dy0, dz0);

      // Contribution (0,0,1)
      T dx1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
      T dy1 = dy0 - SQUISH_CONSTANT;
      T dz1 = dz0 - SQUISH_CONSTANT;
      T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1);
      value += std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb, dx1, dy1, dz1);

      // Contribution (0,1,0)
      T dx2 = dx0 - SQUISH_CONSTANT;
      T dy2 = dy0 - (T)1.0 - SQUISH_CONSTANT;
      T dz2 = dz1;
      T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2);
      value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb, dx2, dy2, dz2);

      // Contribution (1,0,0)
      T dx3 = dx2;
      T dy3 = dy1;
      T dz3 = dz0 - (T)1.0 - SQUISH_CONSTANT;
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3);
      value += std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb, ysb, zsb + 1, dx3, dy3, dz3);

    } else {
      // The point is inside the tetrahedron (3-Simplex) at (1,1,1)

      // Determine which two tetrahedral vertices are the closest
      // out of (1,1,0), (1,0,1), and (0,1,1), but not (1,1,1).
      OSN_BYTE aPoint = 6;
      T aScore = xins;
      OSN_BYTE bPoint = 5;
      T bScore = yins;
      if (aScore <= bScore && zins < bScore) {
        bScore = zins;
        bPoint = 3;
      } else if (aScore > bScore && zins < aScore) {
        aScore = zins;
        aPoint = 3;
      }

      // Determine the two lattice points not part of the tetrahedron that may contribute.
      // This depends on the closest two tetrahedral vertices, including (1,1,1).
      T wins = 3.0 - inSum;
      if (wins < aScore || wins < bScore) {
        // (1,1,1) is one of the closest two tetrahedral vertices.

        // The other closest vertex is the closest of a and b.
        OSN_BYTE c = ((bScore < aScore) ? bPoint : aPoint);

        if (c & 0x01) {
          xsv_ext0 = xsb + 2;
          xsv_ext1 = xsb + 1;
          dx_ext0 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 3);
          dx_ext1 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
        } else {
          xsv_ext0 = xsv_ext1 = xsb;
          dx_ext0 = dx_ext1 = dx0 - (SQUISH_CONSTANT * 3);
        }

        if (c & 0x02) {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy_ext1 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          if (c & 0x01) {
            ysv_ext1 += 1;
            dy_ext1 -= (T)1.0;
          } else {
            ysv_ext0 += 1;
            dy_ext0 -= (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy_ext1 = dy0 - (SQUISH_CONSTANT * 3);
        }

        if (c & 0x04) {
          zsv_ext0 = zsb + 1;
          zsv_ext1 = zsb + 2;
          dz_ext0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          dz_ext1 = dz0 - (T)2.0 - (SQUISH_CONSTANT * 3);
        } else {
          zsv_ext0 = zsv_ext1 = zsb;
          dz_ext0 = dz_ext1 = dz0 - (SQUISH_CONSTANT * 3);
        }
      } else {
        // (1,1,1) is not one of the closest two tetrahedral vertices.

        // The two extra vertices are determined by the closest two.
        OSN_BYTE c = aPoint & bPoint;

        if (c & 0x01) {
          xsv_ext0 = xsb + 1;
          xsv_ext1 = xsb + 2;
          dx_ext0 = dx0 - (T)1.0 - SQUISH_CONSTANT;
          dx_ext1 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 2);
        } else {
          xsv_ext0 = xsv_ext1 = xsb;
          dx_ext0 = dx0 - SQUISH_CONSTANT;
          dx_ext1 = dx0 - (SQUISH_CONSTANT * 2);
        }

        if (c & 0x02) {
          ysv_ext0 = ysb + 1;
          ysv_ext1 = ysb + 2;
          dy_ext0 = dy0 - (T)1.0 - SQUISH_CONSTANT;
          dy_ext1 = dy0 - (T)2.0 - (SQUISH_CONSTANT * 2);
        } else {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy0 - SQUISH_CONSTANT;
          dy_ext1 = dy0 - (SQUISH_CONSTANT * 2);
        }

        if (c & 0x04) {
          zsv_ext0 = zsb + 1;
          zsv_ext1 = zsb + 2;
          dz_ext0 = dz0 - (T)1.0 - SQUISH_CONSTANT;
          dz_ext1 = dz0 - (T)2.0 - (SQUISH_CONSTANT * 2);
        } else {
          zsv_ext0 = zsv_ext1 = zsb;
          dz_ext0 = dz0 - SQUISH_CONSTANT;
          dz_ext1 = dz0 - (SQUISH_CONSTANT * 2);
        }
      }

      // Contribution (1,1,0)
      T dx3 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dy3 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dz3 = dz0 - (SQUISH_CONSTANT * 2);
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3);
      value = std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb, dx3, dy3, dz3);

      // Contribution (1,0,1)
      T dx2 = dx3;
      T dy2 = dy0 - (SQUISH_CONSTANT * 2);
      T dz2 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2);
      value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb + 1, dx2, dy2, dz2);

      // Contribution (0,1,1)
      T dx1 = dx0 - (SQUISH_CONSTANT * 2);
      T dy1 = dy3;
      T dz1 = dz2;
      T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1);
      value += std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb + 1, dx1, dy1, dz1);

      // Contribution (1,1,1)
      dx0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      dy0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      dz0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T attn0 = (dx0 * dx0) + (dy0 * dy0) + (dz0 * dz0);
      value += std::pow(std::max((T)2.0 - attn0, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb + 1, dx0, dy0, dz0);
    }

    // First extra vertex.
    {
      T attn_ext0 = (dx_ext0 * dx_ext0) + (dy_ext0 * dy_ext0) + (dz_ext0 * dz_ext0);
      value += std::pow(std::max((T)2.0 - attn_ext0, (T)0.0), 4) * extrapolate(xsv_ext0, ysv_ext0, zsv_ext0, dx_ext0, dy_ext0, dz_ext0);
    }

    // Second extra vertex.
    {
      T attn_ext1 = (dx_ext1 * dx_ext1) + (dy_ext1 * dy_ext1) + (dz_ext1 * dz_ext1);
      value += std::pow(std::max((T)2.0 - attn_ext1, (T)0.0), 4) * extrapolate(xsv_ext1, ysv_ext1, zsv_ext1, dx_ext1, dy_ext1, dz_ext1);
    }

    return (value / NORM_CONSTANT);
  }

};


// Array of gradient values for 3D. They approximate the directions to the
// vertices of a rhombicuboctahedron from its center, skewed so that the
// triangular and square facets can be inscribed in circles of the same radius.
// New gradient set 2014-10-06.
template <typename T>
const int Noise<3, T>::gradients [] = {
  -11, 4, 4,  -4, 11, 4,  -4, 4, 11,   11, 4, 4,   4, 11, 4,   4, 4, 11,
  -11,-4, 4,  -4,-11, 4,  -4,-4, 11,   11,-4, 4,   4,-11, 4,   4,-4, 11,
  -11, 4,-4,  -4, 11,-4,  -4, 4,-11,   11, 4,-4,   4, 11,-4,   4, 4,-11,
  -11,-4,-4,  -4,-11,-4,  -4,-4,-11,   11,-4,-4,   4,-11,-4,   4,-4,-11
};


// 4D Implementation of the OpenSimplexNoise generator.
template <typename T>
class Noise <4, T> : public NoiseBase {
private:

  // Array of gradient values for 4D. Values are defined below the class definition.
  static const int gradients [];

  inline T extrapolate (long xsb, long ysb, long zsb, long wsb, T dx, T dy, T dz, T dw) const {
    unsigned int index = perm[(perm[(perm[(perm[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF] + wsb) & 0xFF] & 0xFC;
    return gradients[index] * dx +
           gradients[index + 1] * dy +
           gradients[index + 2] * dz +
           gradients[index + 3] * dw;
  }

public:

#ifdef OSN_USE_LONG_LONG
  Noise (long long seed = 0LL) : NoiseBase (seed) {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif
}
#else
  Noise (long seed = 0L) : NoiseBase (seed) {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif
}
#endif
  Noise (const int * p) : NoiseBase (p) {
#ifdef OSN_USE_STATIC_ASSERT
    static_assert(std::is_floating_point<T>::value, "OpenSimplexNoise can only be used with floating-point types");
#endif
}


  T eval (T x, T y, T z, T w) const {

    const T STRETCH_CONSTANT = (T)((1.0 / std::sqrt(4.0 + 1.0) - 1.0) * 0.25);
    const T SQUISH_CONSTANT  = (T)((std::sqrt(4.0 + 1.0) - 1.0) * 0.25);
    const T NORM_CONSTANT    = (T)30.0;

    T dx0, dy0, dz0, dw0;
    long xsb, ysb, zsb, wsb;
    T xins, yins, zins, wins;

    {
      // Place input coordinates on simplectic honeycomb.
      T stretchOffset = (x + y + z + w) * STRETCH_CONSTANT;
      T xs = x + stretchOffset;
      T ys = y + stretchOffset;
      T zs = z + stretchOffset;
      T ws = w + stretchOffset;

      // Floor to get simplectic honeycomb coordinates of rhombo-hypercube origin.
      T xsbd = std::floor(xs);
      T ysbd = std::floor(ys);
      T zsbd = std::floor(zs);
      T wsbd = std::floor(ws);
      xsb = (long)xsbd;
      ysb = (long)ysbd;
      zsb = (long)zsbd;
      wsb = (long)wsbd;

      // Skew out to get actual coordinates of stretched rhombo-hypercube origin.
      T squishOffset = (xsbd + ysbd + zsbd + wsbd) * SQUISH_CONSTANT;
      T xb = xsbd + squishOffset;
      T yb = ysbd + squishOffset;
      T zb = zsbd + squishOffset;
      T wb = wsbd + squishOffset;

      // Positions relative to origin point.
      dx0 = x - xb;
      dy0 = y - yb;
      dz0 = z - zb;
      dw0 = w - wb;

      // Compute simplectic honeycomb coordinates relative to rhombo-hypercube origin.
      xins = xs - xsbd;
      yins = ys - ysbd;
      zins = zs - zsbd;
      wins = ws - wsbd;
    }

    // These are given values inside the next block, and used afterwards.
    long xsv_ext0, ysv_ext0, zsv_ext0, wsv_ext0;
    long xsv_ext1, ysv_ext1, zsv_ext1, wsv_ext1;
    long xsv_ext2, ysv_ext2, zsv_ext2, wsv_ext2;
    T dx_ext0, dy_ext0, dz_ext0, dw_ext0;
    T dx_ext1, dy_ext1, dz_ext1, dw_ext1;
    T dx_ext2, dy_ext2, dz_ext2, dw_ext2;

    T value = 0.0;

    // Sum together to get a value that determines which cell we are in.
    T inSum = xins + yins + zins + wins;

    if (inSum <= (T)1.0) {
      // Inside a pentachoron (4-Simplex) at (0,0,0,0)

      // Determine which two of (0,0,0,1), (0,0,1,0), (0,1,0,0) and (1,0,0,0) are closest.
      OSN_BYTE aPoint = 0x01, bPoint = 0x02;
      T aScore = xins, bScore = yins;
      if (aScore >= bScore && zins > bScore) {
        bPoint = 0x04;
        bScore = zins;
      } else if (aScore < bScore && zins > aScore) {
        aPoint = 0x04;
        aScore = zins;
      }
      if (aScore >= bScore && wins > bScore) {
        bPoint = 0x08;
        bScore = wins;
      } else if (aScore < bScore && wins > aScore) {
        aPoint = 0x08;
        aScore = wins;
      }

      // Determine the three lattice points not part of the pentachoron
      // that may contribute.
      // This depends on the closest two pentachoron vertices, including (0,0,0,0).
      T uins = (T)1.0 - inSum;
      if (uins > aScore || uins > bScore) {
        // (0,0,0,0) is one of the closest two pentachoron vertices.

        // The other closest vertex is the closest out of A and B.
        OSN_BYTE c = (bScore > aScore ? bPoint : aPoint);

        if (c != 0x01) {
          xsv_ext0 = xsb - 1;
          xsv_ext1 = xsv_ext2 = xsb;
          dx_ext0 = dx0 + (T)1.0;
          dx_ext1 = dx_ext2 = dx0;
        } else {
          xsv_ext0 = xsv_ext1 = xsv_ext2 = xsb + 1;
          dx_ext0 = dx_ext1 = dx_ext2 = dx0 - (T)1.0;
        }

        if (c != 0x02) {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb;
          dy_ext0 = dy_ext1 = dy_ext2 = dy0;
          if (c != 0x01) {
            ysv_ext1 -= 1;
            dy_ext1 += (T)1.0;
          } else {
            ysv_ext0 -= 1;
            dy_ext0 += (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb + 1;
          dy_ext0 = dy_ext1 = dy_ext2 = dy0 - (T)1.0;
        }

        if (c != 0x04) {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb;
          dz_ext0 = dz_ext1 = dz_ext2 = dz0;
          if (c & 0x03) {
            zsv_ext1 -= 1;
            dz_ext1 += (T)1.0;
          } else {
            zsv_ext2 -= 1;
            dz_ext2 += (T)1.0;
          }
        } else {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb + 1;
          dz_ext0 = dz_ext1 = dz_ext2 = dz0 - (T)1.0;
        }

        if (c != 0x08) {
          wsv_ext0 = wsv_ext1 = wsb;
          wsv_ext2 = wsb - 1;
          dw_ext0 = dw_ext1 = dw0;
          dw_ext2 = dw0 + (T)1.0;
        } else {
          wsv_ext0 = wsv_ext1 = wsv_ext2 = wsb + 1;
          dw_ext0 = dw_ext1 = dw_ext2 = dw0 - (T)1.0;
        }
      } else {
        // (0,0,0,0) is not one of the closest two pentachoron vertices.

        // The three extra vertices are determined by the closest two.
        OSN_BYTE c = (aPoint | bPoint);

        if (!(c & 0x01)) {
          xsv_ext0 = xsv_ext2 = xsb;
          xsv_ext1 = xsb - 1;
          dx_ext0 = dx0 - (SQUISH_CONSTANT * 2);
          dx_ext1 = dx0 + (T)1.0 - SQUISH_CONSTANT;
          dx_ext2 = dx0 - SQUISH_CONSTANT;
        } else {
          xsv_ext0 = xsv_ext1 = xsv_ext2 = xsb + 1;
          dx_ext0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dx_ext1 = dx_ext2 = dx0 - (T)1.0 - SQUISH_CONSTANT;
        }

        if (!(c & 0x02)) {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb;
          dy_ext0 = dy0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy_ext2 = dy0 - SQUISH_CONSTANT;
          if (c & 0x01) {
            ysv_ext1 -= 1;
            dy_ext1 += (T)1.0;
          } else {
            ysv_ext2 -= 1;
            dy_ext2 += (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb + 1;
          dy_ext0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy_ext2 = dy0 - (T)1.0 - SQUISH_CONSTANT;
        }

        if (!(c & 0x04)) {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb;
          dz_ext0 = dz0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz_ext2 = dz0 - SQUISH_CONSTANT;
          if (c & 0x03) {
            zsv_ext1 -= 1;
            dz_ext1 += (T)1.0;
          } else {
            zsv_ext2 -= 1;
            dz_ext2 += (T)1.0;
          }
        } else {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb + 1;
          dz_ext0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz_ext2 = dz0 - (T)1.0 - SQUISH_CONSTANT;
        }

        if (!(c & 0x08)) {
          wsv_ext0 = wsv_ext1 = wsb;
          wsv_ext2 = wsb - 1;
          dw_ext0 = dw0 - (SQUISH_CONSTANT * 2);
          dw_ext1 = dw0 - SQUISH_CONSTANT;
          dw_ext2 = dw0 + (T)1.0 - SQUISH_CONSTANT;
        } else {
          wsv_ext0 = wsv_ext1 = wsv_ext2 = wsb + 1;
          dw_ext0 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dw_ext1 = dw_ext2 = dw0 - (T)1.0 - SQUISH_CONSTANT;
        }
      }

      // Contribution (0,0,0,0).
      T attn0 = (dx0 * dx0) + (dy0 * dy0) + (dz0 * dz0) + (dw0 * dw0);
      value = std::pow(std::max((T)2.0 - attn0, (T)0.0), 4) * extrapolate(xsb, ysb, zsb, wsb, dx0, dy0, dz0, dw0);

      // Contribution (1,0,0,0).
      T dx1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
      T dy1 = dy0 - SQUISH_CONSTANT;
      T dz1 = dz0 - SQUISH_CONSTANT;
      T dw1 = dw0 - SQUISH_CONSTANT;
      T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1) + (dw1 * dw1);
      value += std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb, wsb, dx1, dy1, dz1, dw1);

      // Contribution (0,1,0,0).
      T dx2 = dx0 - SQUISH_CONSTANT;
      T dy2 = dy0 - (T)1.0 - SQUISH_CONSTANT;
      T dz2 = dz1;
      T dw2 = dw1;
      T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2) + (dw2 * dw2);
      value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb, wsb, dx2, dy2, dz2, dw2);

      // Contribution (0,0,1,0).
      T dx3 = dx2;
      T dy3 = dy1;
      T dz3 = dz0 - (T)1.0 - SQUISH_CONSTANT;
      T dw3 = dw1;
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3) + (dw3 * dw3);
      value += std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb, ysb, zsb + 1, wsb, dx3, dy3, dz3, dw3);

      // Contribution (0,0,0,1).
      T dx4 = dx2;
      T dy4 = dy1;
      T dz4 = dz1;
      T dw4 = dw0 - (T)1.0 - SQUISH_CONSTANT;
      T attn4 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4) + (dw4 * dw4);
      value += std::pow(std::max((T)2.0 - attn4, (T)0.0), 4) * extrapolate(xsb, ysb, zsb, wsb + 1, dx4, dy4, dz4, dw4);

    } else if (inSum >= 3.0) {
      // Inside the pentachoron (4-simplex) at (1,1,1,1).

      // Determine which two of (1,1,1,0), (1,1,0,1), (1,0,1,1), (0,1,1,1) are closest.
      OSN_BYTE aPoint = 0x0E;
      T aScore = xins;
      OSN_BYTE bPoint = 0x0D;
      T bScore = yins;
      if (aScore <= bScore && zins < bScore) {
        bPoint = 0x0B;
        bScore = zins;
      } else if (aScore > bScore && zins < aScore) {
        aPoint = 0x0B;
        aScore = zins;
      }
      if (aScore <= bScore && wins < bScore) {
        bPoint = 0x07;
        bScore = wins;
      } else if (aScore > bScore && wins < aScore) {
        aPoint = 0x07;
        aScore = wins;
      }

      // Determine the three lattice points not part of the pentachoron that may contribute.
      // This depends on the closest two pentachoron vertices, including (0,0,0,0).
      T uins = 4.0 - inSum;
      if (uins < aScore || uins < bScore) {
        // (1,1,1,1) is one of the closest two pentachoron vertices.

        // The other closest vertex is the closest out of A and B.
        OSN_BYTE c = (bScore < aScore ? bPoint : aPoint);
        if (c & 0x01) {
          xsv_ext0 = xsb + 2;
          xsv_ext1 = xsv_ext2 = xsb + 1;
          dx_ext0 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 4);
          dx_ext1 = dx_ext2 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 4);
        } else {
          xsv_ext0 = xsv_ext1 = xsv_ext2 = xsb;
          dx_ext0 = dx_ext1 = dx_ext2 = dx0 - (SQUISH_CONSTANT * 4);
        }

        if (c & 0x02) {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb + 1;
          dy_ext0 = dy_ext1 = dy_ext2 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 4);
          if (c & 0x01) {
            ysv_ext1 += 1;
            dy_ext1 -= (T)1.0;
          } else {
            ysv_ext0 += 1;
            dy_ext0 -= (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb;
          dy_ext0 = dy_ext1 = dy_ext2 = dy0 - (SQUISH_CONSTANT * 4);
        }

        if (c & 0x04) {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb + 1;
          dz_ext0 = dz_ext1 = dz_ext2 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 4);
          if ((c & 0x03) != 0x03) {
            if (!(c & 0x03)) {
              zsv_ext0 += 1;
              dz_ext0 -= (T)1.0;
            } else {
              zsv_ext1 += 1;
              dz_ext1 -= (T)1.0;
            }
          } else {
            zsv_ext2 += 1;
            dz_ext2 -= (T)1.0;
          }
        } else {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb;
          dz_ext0 = dz_ext1 = dz_ext2 = dz0 - (SQUISH_CONSTANT * 4);
        }

        if (c & 0x08) {
          wsv_ext0 = wsv_ext1 = wsb + 1;
          wsv_ext2 = wsb + 2;
          dw_ext0 = dw_ext1 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 4);
          dw_ext2 = dw0 - (T)2.0 - (SQUISH_CONSTANT * 4);
        } else {
          wsv_ext0 = wsv_ext1 = wsv_ext2 = wsb;
          dw_ext0 = dw_ext1 = dw_ext2 = dw0 - (SQUISH_CONSTANT * 4);
        }
      } else {
        // (1,1,1,1) is not one of the closest two pentachoron vertices.

        OSN_BYTE c = aPoint & bPoint;
        if (c & 0x01) {
          xsv_ext0 = xsv_ext2 = xsb + 1;
          xsv_ext1 = xsb + 2;
          dx_ext0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dx_ext1 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 3);
          dx_ext2 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
        } else {
          xsv_ext0 = xsv_ext1 = xsv_ext2 = xsb;
          dx_ext0 = dx0 - (SQUISH_CONSTANT * 2);
          dx_ext1 = dx_ext2 = dx0 - (SQUISH_CONSTANT * 3);
        }

        if (c & 0x02) {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb + 1;
          dy_ext0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy_ext2 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          if (c & 0x01) {
            ysv_ext2 += 1;
            dy_ext2 -= (T)1.0;
          } else {
            ysv_ext1 += 1;
            dy_ext1 -= (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysv_ext2 = ysb;
          dy_ext0 = dy0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy_ext2 = dy0 - (SQUISH_CONSTANT * 3);
        }

        if (c & 0x04) {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb + 1;
          dz_ext0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz_ext2 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          if (c & 0x03) {
            zsv_ext2 += 1;
            dz_ext2 -= (T)1.0;
          } else {
            zsv_ext1 += 1;
            dz_ext1 -= (T)1.0;
          }
        } else {
          zsv_ext0 = zsv_ext1 = zsv_ext2 = zsb;
          dz_ext0 = dz0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz_ext2 = dz0 - (SQUISH_CONSTANT * 3);
        }

        if (c & 0x08) {
          wsv_ext0 = wsv_ext1 = wsb + 1;
          wsv_ext2 = wsb + 2;
          dw_ext0 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dw_ext1 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          dw_ext2 = dw0 - (T)2.0 - (SQUISH_CONSTANT * 3);
        } else {
          wsv_ext0 = wsv_ext1 = wsv_ext2 = wsb;
          dw_ext0 = dw0 - (SQUISH_CONSTANT * 2);
          dw_ext1 = dw_ext2 = dw0 - (SQUISH_CONSTANT * 3);
        }
      }

      // Contribution (1,1,1,0).
      T dx4 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T dy4 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T dz4 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T dw4 = dw0 - (SQUISH_CONSTANT * 3);
      T attn4 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4) + (dw4 * dw4);
      value = std::pow(std::max((T)2.0 - attn4, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb + 1, wsb, dx4, dy4, dz4, dw4);

      // Contribution (1,1,0,1).
      T dx3 = dx4;
      T dy3 = dy4;
      T dz3 = dz0 - (SQUISH_CONSTANT * 3);
      T dw3 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3) + (dw3 * dw3);
      value += std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb, wsb + 1, dx3, dy3, dz3, dw3);

      // Contribution (1,0,1,1).
      T dx2 = dx4;
      T dy2 = dy0 - (SQUISH_CONSTANT * 3);
      T dz2 = dz4;
      T dw2 = dw3;
      T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2) + (dw2 * dw2);
      value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb + 1, wsb + 1, dx2, dy2, dz2, dw2);

      // Contribution (0,1,1,1).
      T dx1 = dx0 - (SQUISH_CONSTANT * 3);
      T dy1 = dy4;
      T dz1 = dz4;
      T dw1 = dw3;
      T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1) + (dw1 * dw1);
      value += std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb + 1, wsb + 1, dx1, dy1, dz1, dw1);

      // Contribution (1,1,1,1).
      dx0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 4);
      dy0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 4);
      dz0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 4);
      dw0 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 4);
      T attn0 = (dx0 * dx0) + (dy0 * dy0) + (dz0 * dz0) + (dw0 * dw0);
      value += std::pow(std::max((T)2.0 - attn0, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb + 1, wsb + 1, dx0, dy0, dz0, dw0);

    } else if (inSum <= (T)2.0) {
      // Inside the first dispentachoron (rectified 4-simplex).

      T aScore;
      OSN_BYTE aPoint;
      bool aIsBiggerSide = true;
      T bScore;
      OSN_BYTE bPoint;
      bool bIsBiggerSide = true;

      // Decide between (1,1,0,0) and (0,0,1,1).
      if (xins + yins > zins + wins) {
        aPoint = 0x03;
        aScore = xins + yins;
      } else {
        aPoint = 0x0C;
        aScore = zins + wins;
      }

      // Decide between (1,0,1,0) and (0,1,0,1).
      if (xins + zins > yins + wins) {
        bPoint = 0x05;
        bScore = xins + zins;
      } else {
        bPoint = 0x0A;
        bScore = yins + wins;
      }

      // Closer of (1,0,0,1) and (0,1,1,0) will replace the further of A and B, if closer.
      if (xins + wins > yins + zins) {
        T score = xins + wins;
        if (aScore >= bScore && score > bScore) {
          bPoint = 0x09;
          bScore = score;
        } else if (aScore < bScore && score > aScore) {
          aPoint = 0x09;
          aScore = score;
        }
      } else {
        T score = yins + zins;
        if (aScore >= bScore && score > bScore) {
          bPoint = 0x06;
          bScore = score;
        } else if (aScore < bScore && score > aScore) {
          aPoint = 0x06;
          aScore = score;
        }
      }

      // Decide if (1,0,0,0) is closer.
      T p1 = (T)2.0 - inSum + xins;
      if (aScore >= bScore && p1 > bScore) {
        bPoint = 0x01;
        bScore = p1;
        bIsBiggerSide = false;
      } else if (aScore < bScore && p1 > aScore) {
        aPoint = 0x01;
        aScore = p1;
        aIsBiggerSide = false;
      }

      // Decide if (0,1,0,0) is closer.
      T p2 = (T)2.0 - inSum + yins;
      if (aScore >= bScore && p2 > bScore) {
        bPoint = 0x02;
        bScore = p2;
        bIsBiggerSide = false;
      } else if (aScore < bScore && p2 > aScore) {
        aPoint = 0x02;
        aScore = p2;
        aIsBiggerSide = false;
      }

      // Decide if (0,0,1,0) is closer.
      T p3 = (T)2.0 - inSum + zins;
      if (aScore >= bScore && p3 > bScore) {
        bPoint = 0x04;
        bScore = p3;
        bIsBiggerSide = false;
      } else if (aScore < bScore && p3 > aScore) {
        aPoint = 0x04;
        aScore = p3;
        aIsBiggerSide = false;
      }

      // Decide if (0,0,0,1) is closer.
      T p4 = (T)2.0 - inSum + wins;
      if (aScore >= bScore && p4 > bScore) {
        bPoint = 0x08;
        bScore = p4;
        bIsBiggerSide = false;
      } else if (aScore < bScore && p4 > aScore) {
        aPoint = 0x08;
        aScore = p4;
        aIsBiggerSide = false;
      }

      // Where each of the two closest points are determines how the extra three vertices are calculated.
      if (aIsBiggerSide == bIsBiggerSide) {
        if (aIsBiggerSide) {
          // Both closest points are on the bigger side.

          OSN_BYTE c1 = aPoint | bPoint;
          OSN_BYTE c2 = aPoint & bPoint;
          if (!(c1 & 0x01)) {
            xsv_ext0 = xsb;
            xsv_ext1 = xsb - 1;
            dx_ext0 = dx0 - (SQUISH_CONSTANT * 3);
            dx_ext1 = dx0 + (T)1.0 - (SQUISH_CONSTANT * 2);
          } else {
            xsv_ext0 = xsv_ext1 = xsb + 1;
            dx_ext0 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            dx_ext1 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          }

          if (!(c1 & 0x02)) {
            ysv_ext0 = ysb;
            ysv_ext1 = ysb - 1;
            dy_ext0 = dy0 - (SQUISH_CONSTANT * 3);
            dy_ext1 = dy0 + (T)1.0 - (SQUISH_CONSTANT * 2);
          } else {
            ysv_ext0 = ysv_ext1 = ysb + 1;
            dy_ext0 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            dy_ext1 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          }

          if (!(c1 & 0x04)) {
            zsv_ext0 = zsb;
            zsv_ext1 = zsb - 1;
            dz_ext0 = dz0 - (SQUISH_CONSTANT * 3);
            dz_ext1 = dz0 + (T)1.0 - (SQUISH_CONSTANT * 2);
          } else {
            zsv_ext0 = zsv_ext1 = zsb + 1;
            dz_ext0 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            dz_ext1 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          }

          if (!(c1 & 0x08)) {
            wsv_ext0 = wsb;
            wsv_ext1 = wsb - 1;
            dw_ext0 = dw0 - (SQUISH_CONSTANT * 3);
            dw_ext1 = dw0 + (T)1.0 - (SQUISH_CONSTANT * 2);
          } else {
            wsv_ext0 = wsv_ext1 = wsb + 1;
            dw_ext0 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            dw_ext1 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          }

          // One combination is a permutation of (0,0,0,2) based on c2.
          xsv_ext2 = xsb;
          ysv_ext2 = ysb;
          zsv_ext2 = zsb;
          wsv_ext2 = wsb;
          dx_ext2 = dx0 - (SQUISH_CONSTANT * 2);
          dy_ext2 = dy0 - (SQUISH_CONSTANT * 2);
          dz_ext2 = dz0 - (SQUISH_CONSTANT * 2);
          dw_ext2 = dw0 - (SQUISH_CONSTANT * 2);
          if (c2 & 0x01) {
            xsv_ext2 += 2;
            dx_ext2 -= (T)2.0;
          } else if (c2 & 0x02) {
            ysv_ext2 += 2;
            dy_ext2 -= (T)2.0;
          } else if (c2 & 0x04) {
            zsv_ext2 += 2;
            dz_ext2 -= (T)2.0;
          } else {
            wsv_ext2 += 2;
            dw_ext2 -= (T)2.0;
          }
        } else {
          // Both closest points are on the smaller side.

          // One of the two extra points is (0,0,0,0).
          xsv_ext2 = xsb;
          ysv_ext2 = ysb;
          zsv_ext2 = zsb;
          wsv_ext2 = wsb;
          dx_ext2 = dx0;
          dy_ext2 = dy0;
          dz_ext2 = dz0;
          dw_ext2 = dw0;

          // The other two points are based on the omitted axes.
          OSN_BYTE c = aPoint | bPoint;
          if (!(c & 0x01)) {
            xsv_ext0 = xsb - 1;
            xsv_ext1 = xsb;
            dx_ext0 = dx0 + (T)1.0 - SQUISH_CONSTANT;
            dx_ext1 = dx0 - SQUISH_CONSTANT;
          } else {
            xsv_ext0 = xsv_ext1 = xsb + 1;
            dx_ext0 = dx_ext1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
          }

          if (!(c & 0x02)) {
            ysv_ext0 = ysv_ext1 = ysb;
            dy_ext0 = dy_ext1 = dy0 - SQUISH_CONSTANT;
            if (c & 0x01) {
              ysv_ext0 -= 1;
              dy_ext0 += (T)1.0;
            } else {
              ysv_ext1 -= 1;
              dy_ext1 += (T)1.0;
            }
          } else {
            ysv_ext0 = ysv_ext1 = ysb + 1;
            dy_ext0 = dy_ext1 = dy0 - (T)1.0 - SQUISH_CONSTANT;
          }

          if (!(c & 0x04)) {
            zsv_ext0 = zsv_ext1 = zsb;
            dz_ext0 = dz_ext1 = dz0 - SQUISH_CONSTANT;
            if ((c & 0x03) == 0x03)
            {
              zsv_ext0 -= 1;
              dz_ext0 += (T)1.0;
            } else {
              zsv_ext1 -= 1;
              dz_ext1 += (T)1.0;
            }
          } else {
            zsv_ext0 = zsv_ext1 = zsb + 1;
            dz_ext0 = dz_ext1 = dz0 - (T)1.0 - SQUISH_CONSTANT;
          }

          if (!(c & 0x08)) {
            wsv_ext0 = wsb;
            wsv_ext1 = wsb - 1;
            dw_ext0 = dw0 - SQUISH_CONSTANT;
            dw_ext1 = dw0 + (T)1.0 - SQUISH_CONSTANT;
          } else {
            wsv_ext0 = wsv_ext1 = wsb + 1;
            dw_ext0 = dw_ext1 = dw0 - (T)1.0 - SQUISH_CONSTANT;
          }
        }
      } else {
        // One point on each side.

        OSN_BYTE c1, c2;
        if (aIsBiggerSide) {
          c1 = aPoint;
          c2 = bPoint;
        } else {
          c1 = bPoint;
          c2 = aPoint;
        }

        // Two contributions are the bigger-sided point with each 0 replaced with -1.
        if (!(c1 & 0x01)) {
          xsv_ext0 = xsb - 1;
          xsv_ext1 = xsb;
          dx_ext0 = dx0 + (T)1.0 - SQUISH_CONSTANT;
          dx_ext1 = dx0 - SQUISH_CONSTANT;
        } else {
          xsv_ext0 = xsv_ext1 = xsb + 1;
          dx_ext0 = dx_ext1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
        }

        if (!(c1 & 0x02)) {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy_ext1 = dy0 - SQUISH_CONSTANT;
          if ((c1 & 0x01) == 0x01) {
            ysv_ext0 -= 1;
            dy_ext0 += (T)1.0;
          } else {
            ysv_ext1 -= 1;
            dy_ext1 += (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy_ext1 = dy0 - (T)1.0 - SQUISH_CONSTANT;
        }

        if (!(c1 & 0x04)) {
          zsv_ext0 = zsv_ext1 = zsb;
          dz_ext0 = dz_ext1 = dz0 - SQUISH_CONSTANT;
          if ((c1 & 0x03) == 0x03) {
            zsv_ext0 -= 1;
            dz_ext0 += (T)1.0;
          } else {
            zsv_ext1 -= 1;
            dz_ext1 += (T)1.0;
          }
        } else {
          zsv_ext0 = zsv_ext1 = zsb + 1;
          dz_ext0 = dz_ext1 = dz0 - (T)1.0 - SQUISH_CONSTANT;
        }

        if (!(c1 & 0x08)) {
          wsv_ext0 = wsb;
          wsv_ext1 = wsb - 1;
          dw_ext0 = dw0 - SQUISH_CONSTANT;
          dw_ext1 = dw0 + (T)1.0 - SQUISH_CONSTANT;
        } else {
          wsv_ext0 = wsv_ext1 = wsb + 1;
          dw_ext0 = dw_ext1 = dw0 - (T)1.0 - SQUISH_CONSTANT;
        }

        // One contribution is a permutation of (0,0,0,2) based on the smaller-sided point.
        xsv_ext2 = xsb;
        ysv_ext2 = ysb;
        zsv_ext2 = zsb;
        wsv_ext2 = wsb;
        dx_ext2 = dx0 - (SQUISH_CONSTANT * 2);
        dy_ext2 = dy0 - (SQUISH_CONSTANT * 2);
        dz_ext2 = dz0 - (SQUISH_CONSTANT * 2);
        dw_ext2 = dw0 - (SQUISH_CONSTANT * 2);
        if ((c2 & 0x01) != 0) {
          xsv_ext2 += 2;
          dx_ext2 -= (T)2.0;
        } else if ((c2 & 0x02) != 0) {
          ysv_ext2 += 2;
          dy_ext2 -= (T)2.0;
        } else if ((c2 & 0x04) != 0) {
          zsv_ext2 += 2;
          dz_ext2 -= (T)2.0;
        } else {
          wsv_ext2 += 2;
          dw_ext2 -= (T)2.0;
        }
      }

      //Contribution (1,0,0,0)
      T dx1 = dx0 - (T)1.0 - SQUISH_CONSTANT;
      T dy1 = dy0 - SQUISH_CONSTANT;
      T dz1 = dz0 - SQUISH_CONSTANT;
      T dw1 = dw0 - SQUISH_CONSTANT;
      T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1) + (dw1 * dw1);
      value += std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb, wsb, dx1, dy1, dz1, dw1);

      //Contribution (0,1,0,0)
      T dx2 = dx0 - SQUISH_CONSTANT;
      T dy2 = dy0 - (T)1.0 - SQUISH_CONSTANT;
      T dz2 = dz1;
      T dw2 = dw1;
      T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2) + (dw2 * dw2);
      value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb, wsb, dx2, dy2, dz2, dw2);

      //Contribution (0,0,1,0)
      T dx3 = dx2;
      T dy3 = dy1;
      T dz3 = dz0 - (T)1.0 - SQUISH_CONSTANT;
      T dw3 = dw1;
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3) + (dw3 * dw3);
      value += std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb, ysb, zsb + 1, wsb, dx3, dy3, dz3, dw3);

      //Contribution (0,0,0,1)
      T dx4 = dx2;
      T dy4 = dy1;
      T dz4 = dz1;
      T dw4 = dw0 - (T)1.0 - SQUISH_CONSTANT;
      T attn4 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4) + (dw4 * dw4);
      value += std::pow(std::max((T)2.0 - attn4, (T)0.0), 4) * extrapolate(xsb, ysb, zsb, wsb + 1, dx4, dy4, dz4, dw4);

      //Contribution (1,1,0,0)
      T dx5 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dy5 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dz5 = dz0 - (SQUISH_CONSTANT * 2);
      T dw5 = dw0 - (SQUISH_CONSTANT * 2);
      T attn5 = (dx5 * dx5) + (dy5 * dy5) + (dz5 * dz5) + (dw5 * dw5);
      value += std::pow(std::max((T)2.0 - attn5, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb, wsb, dx5, dy5, dz5, dw5);

      //Contribution (1,0,1,0)
      T dx6 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dy6 = dy0 - (SQUISH_CONSTANT * 2);
      T dz6 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dw6 = dw0 - (SQUISH_CONSTANT * 2);
      T attn6 = (dx6 * dx6) + (dy6 * dy6) + (dz6 * dz6) + (dw6 * dw6);
      value += std::pow(std::max((T)2.0 - attn6, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb + 1, wsb, dx6, dy6, dz6, dw6);

      //Contribution (1,0,0,1)
      T dx7 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dy7 = dy0 - (SQUISH_CONSTANT * 2);
      T dz7 = dz0 - (SQUISH_CONSTANT * 2);
      T dw7 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T attn7 = (dx7 * dx7) + (dy7 * dy7) + (dz7 * dz7) + (dw7 * dw7);
      value += std::pow(std::max((T)2.0 - attn7, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb, wsb + 1, dx7, dy7, dz7, dw7);

      // Contribution (0,1,1,0).
      T dx8 = dx0 - (SQUISH_CONSTANT * 2);
      T dy8 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dz8 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dw8 = dw0 - (SQUISH_CONSTANT * 2);
      T attn8 = (dx8 * dx8) + (dy8 * dy8) + (dz8 * dz8) + (dw8 * dw8);
      value += std::pow(std::max((T)2.0 - attn8, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb + 1, wsb, dx8, dy8, dz8, dw8);

      // Contribution (0,1,0,1).
      T dx9 = dx0 - (SQUISH_CONSTANT * 2);
      T dy9 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dz9 = dz0 - (SQUISH_CONSTANT * 2);
      T dw9 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T attn9 = (dx9 * dx9) + (dy9 * dy9) + (dz9 * dz9) + (dw9 * dw9);
      value += std::pow(std::max((T)2.0 - attn9, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb, wsb + 1, dx9, dy9, dz9, dw9);

      // Contribution (0,0,1,1).
      T dx10 = dx0 - 2 * SQUISH_CONSTANT;
      T dy10 = dy0 - 2 * SQUISH_CONSTANT;
      T dz10 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T dw10 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
      T attn10 = (dx10 * dx10) + (dy10 * dy10) + (dz10 * dz10) + (dw10 * dw10);
      value += std::pow(std::max((T)2.0 - attn10, (T)0.0), 4) * extrapolate(xsb, ysb, zsb + 1, wsb + 1, dx10, dy10, dz10, dw10);

    } else {
      // Inside the second dispentachoron (rectified 4-simplex).

      OSN_BYTE aPoint, bPoint;
      T aScore, bScore;
      bool aIsBiggerSide(true), bIsBiggerSide(true);

      // Decide between (0,0,1,1) and (1,1,0,0).
      if (xins + yins < zins + wins) {
        aPoint = 0x0C;
        aScore = xins + yins;
      } else {
        aPoint = 0x03;
        aScore = zins + wins;
      }

      //Decide between (0,1,0,1) and (1,0,1,0).
      if (xins + zins < yins + wins) {
        bPoint = 0x0A;
        bScore = xins + zins;
      } else {
        bPoint = 0x05;
        bScore = yins + wins;
      }

      // The closer of (0,1,1,0) and (1,0,0,1) will replace the further of a and b, if closer.
      if (xins + wins < yins + zins) {
        T score(xins + wins);
        if (aScore <= bScore && score < bScore) {
          bPoint = 0x06;
          bScore = score;
        } else if (aScore > bScore && score < aScore) {
          aPoint = 0x06;
          aScore = score;
        }
      } else {
        T score(yins + zins);
        if (aScore <= bScore && score < bScore) {
          bPoint = 0x09;
          bScore = score;
        } else if (aScore > bScore && score < aScore) {
          aPoint = 0x09;
          aScore = score;
        }
      }

      // Decide if (0,1,1,1) is closer.
      {
        T p1 = 3.0 - inSum + xins;
        if (aScore <= bScore && p1 < bScore) {
          bPoint = 0x0E;
          bScore = p1;
          bIsBiggerSide = false;
        } else if (aScore > bScore && p1 < aScore) {
          aPoint = 0x0E;
          aScore = p1;
          aIsBiggerSide = false;
        }
      }

      // Decide if (1,0,1,1) is closer.
      {
        T p2 = 3.0 - inSum + yins;
        if (aScore <= bScore && p2 < bScore) {
          bPoint = 0x0D;
          bScore = p2;
          bIsBiggerSide = false;
        } else if (aScore > bScore && p2 < aScore) {
          aPoint = 0x0D;
          aScore = p2;
          aIsBiggerSide = false;
        }
      }

      // Decide if (1,1,0,1) is closer.
      {
        T p3 = 3.0 - inSum + zins;
        if (aScore <= bScore && p3 < bScore) {
          bPoint = 0x0B;
          bScore = p3;
          bIsBiggerSide = false;
        } else if (aScore > bScore && p3 < aScore) {
          aPoint = 0x0B;
          aScore = p3;
          aIsBiggerSide = false;
        }
      }

      // Decide if (1,1,1,0) is closer.
      {
        T p4 = 3.0 - inSum + wins;
        if (aScore <= bScore && p4 < bScore) {
          bPoint = 0x07;
          bScore = p4;
          bIsBiggerSide = false;
        } else if (aScore > bScore && p4 < aScore) {
          aPoint = 0x07;
          aScore = p4;
          aIsBiggerSide = false;
        }
      }

      // Where each of the two closest points are determines how the extra three vertices are calculated.
      if (aIsBiggerSide == bIsBiggerSide) {
        if (aIsBiggerSide) {
          // Both closest points are on the bigger side.

          OSN_BYTE c1 = aPoint & bPoint;
          OSN_BYTE c2 = aPoint | bPoint;

          // Two contributions are permutations of (0,0,0,1) and (0,0,0,2) based on c1.
          xsv_ext0 = xsv_ext1 = xsb;
          ysv_ext0 = ysv_ext1 = ysb;
          zsv_ext0 = zsv_ext1 = zsb;
          wsv_ext0 = wsv_ext1 = wsb;
          dx_ext0 = dx0 - SQUISH_CONSTANT;
          dy_ext0 = dy0 - SQUISH_CONSTANT;
          dz_ext0 = dz0 - SQUISH_CONSTANT;
          dw_ext0 = dw0 - SQUISH_CONSTANT;
          dx_ext1 = dx0 - (SQUISH_CONSTANT * 2);
          dy_ext1 = dy0 - (SQUISH_CONSTANT * 2);
          dz_ext1 = dz0 - (SQUISH_CONSTANT * 2);
          dw_ext1 = dw0 - (SQUISH_CONSTANT * 2);

          if (c1 & 0x01) {
            xsv_ext0 += 1;
            dx_ext0 -= (T)1.0;
            xsv_ext1 += 2;
            dx_ext1 -= (T)2.0;
          } else if (c1 & 0x02) {
            ysv_ext0 += 1;
            dy_ext0 -= (T)1.0;
            ysv_ext1 += 2;
            dy_ext1 -= (T)2.0;
          } else if (c1 & 0x04) {
            zsv_ext0 += 1;
            dz_ext0 -= (T)1.0;
            zsv_ext1 += 2;
            dz_ext1 -= (T)2.0;
          } else {
            wsv_ext0 += 1;
            dw_ext0 -= (T)1.0;
            wsv_ext1 += 2;
            dw_ext1 -= (T)2.0;
          }

          // One contribution is a permutation of (1,1,1,-1) based on c2.
          xsv_ext2 = xsb + 1;
          ysv_ext2 = ysb + 1;
          zsv_ext2 = zsb + 1;
          wsv_ext2 = wsb + 1;
          dx_ext2 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dy_ext2 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dz_ext2 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          dw_ext2 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
          if (!(c2 & 0x01)) {
            xsv_ext2 -= 2;
            dx_ext2 += (T)2.0;
          } else if (!(c2 & 0x02)) {
            ysv_ext2 -= 2;
            dy_ext2 += (T)2.0;
          } else if (!(c2 & 0x04)) {
            zsv_ext2 -= 2;
            dz_ext2 += (T)2.0;
          } else {
            wsv_ext2 -= 2;
            dw_ext2 += (T)2.0;
          }
        } else {
          // Both closest points are on the smaller side.

          // One of the two extra points is (1,1,1,1).
          xsv_ext2 = xsb + 1;
          ysv_ext2 = ysb + 1;
          zsv_ext2 = zsb + 1;
          wsv_ext2 = wsb + 1;
          dx_ext2 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 4);
          dy_ext2 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 4);
          dz_ext2 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 4);
          dw_ext2 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 4);

          // The other two points are based on the shared axes.
          OSN_BYTE c = aPoint & bPoint;
          if (c & 0x01) {
            xsv_ext0 = xsb + 2;
            xsv_ext1 = xsb + 1;
            dx_ext0 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 3);
            dx_ext1 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          } else {
            xsv_ext0 = xsv_ext1 = xsb;
            dx_ext0 = dx_ext1 = dx0 - (SQUISH_CONSTANT * 3);
          }

          if (c & 0x02) {
            ysv_ext0 = ysv_ext1 = ysb + 1;
            dy_ext0 = dy_ext1 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            if (!(c & 0x01)) {
              ysv_ext0 += 1;
              dy_ext0 -= (T)1.0;
            } else {
              ysv_ext1 += 1;
              dy_ext1 -= (T)1.0;
            }
          } else {
            ysv_ext0 = ysv_ext1 = ysb;
            dy_ext0 = dy_ext1 = dy0 - (SQUISH_CONSTANT * 3);
          }

          if (c & 0x04) {
            zsv_ext0 = zsv_ext1 = zsb + 1;
            dz_ext0 = dz_ext1 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            if (!(c & 0x03)) {
              zsv_ext0 += 1;
              dz_ext0 -= (T)1.0;
            } else {
              zsv_ext1 += 1;
              dz_ext1 -= (T)1.0;
            }
          } else {
            zsv_ext0 = zsv_ext1 = zsb;
            dz_ext0 = dz_ext1 = dz0 - (SQUISH_CONSTANT * 3);
          }

          if (c & 0x08) {
            wsv_ext0 = wsb + 1;
            wsv_ext1 = wsb + 2;
            dw_ext0 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 3);
            dw_ext1 = dw0 - (T)2.0 - (SQUISH_CONSTANT * 3);
          } else {
            wsv_ext0 = wsv_ext1 = wsb;
            dw_ext0 = dw_ext1 = dw0 - (SQUISH_CONSTANT * 3);
          }
        }
      } else {
        // One point on each "side".

        OSN_BYTE c1, c2;
        if (aIsBiggerSide) {
          c1 = aPoint;
          c2 = bPoint;
        } else {
          c1 = bPoint;
          c2 = aPoint;
        }

        // Two contributions are the bigger-sided point with each 1 replaced with 2.
        if (c1 & 0x01) {
          xsv_ext0 = xsb + 2;
          xsv_ext1 = xsb + 1;
          dx_ext0 = dx0 - (T)2.0 - (SQUISH_CONSTANT * 3);
          dx_ext1 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
        } else {
          xsv_ext0 = xsv_ext1 = xsb;
          dx_ext0 = dx_ext1 = dx0 - (SQUISH_CONSTANT * 3);
        }

        if (c1 & 0x02) {
          ysv_ext0 = ysv_ext1 = ysb + 1;
          dy_ext0 = dy_ext1 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          if (!(c1 & 0x01)) {
            ysv_ext0 += 1;
            dy_ext0 -= (T)1.0;
          } else {
            ysv_ext1 += 1;
            dy_ext1 -= (T)1.0;
          }
        } else {
          ysv_ext0 = ysv_ext1 = ysb;
          dy_ext0 = dy_ext1 = dy0 - (SQUISH_CONSTANT * 3);
        }

        if (c1 & 0x04) {
          zsv_ext0 = zsv_ext1 = zsb + 1;
          dz_ext0 = dz_ext1 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          if (!(c1 & 0x03)) {
            zsv_ext0 += 1;
            dz_ext0 -= (T)1.0;
          } else {
            zsv_ext1 += 1;
            dz_ext1 -= (T)1.0;
          }
        } else {
          zsv_ext0 = zsv_ext1 = zsb;
          dz_ext0 = dz_ext1 = dz0 - 3 * (SQUISH_CONSTANT * 3);
        }

        if (c1 & 0x08) {
          wsv_ext0 = wsb + 1;
          wsv_ext1 = wsb + 2;
          dw_ext0 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 3);
          dw_ext1 = dw0 - (T)2.0 - (SQUISH_CONSTANT * 3);
        } else {
          wsv_ext0 = wsv_ext1 = wsb;
          dw_ext0 = dw_ext1 = dw0 - (SQUISH_CONSTANT * 3);
        }

        //  One contribution is a permutation of (1,1,1,-1) based on the smaller-sided point.
        xsv_ext2 = xsb + 1;
        ysv_ext2 = ysb + 1;
        zsv_ext2 = zsb + 1;
        wsv_ext2 = wsb + 1;
        dx_ext2 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        dy_ext2 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        dz_ext2 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        dw_ext2 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        if (!(c2 & 0x01)) {
          xsv_ext2 -= 2;
          dx_ext2 += (T)2.0;
        } else if (!(c2 & 0x02)) {
          ysv_ext2 -= 2;
          dy_ext2 += (T)2.0;
        } else if (!(c2 & 0x04)) {
          zsv_ext2 -= 2;
          dz_ext2 += (T)2.0;
        } else {
          wsv_ext2 -= 2;
          dw_ext2 += (T)2.0;
        }
      }

      // Contribution (1,1,1,0).
      T dx4 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T dy4 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T dz4 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T dw4 = dw0 - (SQUISH_CONSTANT * 3);
      T attn4 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4) + (dw4 * dw4);
      value += std::pow(std::max((T)2.0 - attn4, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb + 1, wsb, dx4, dy4, dz4, dw4);

      //Contribution (1,1,0,1).
      T dx3 = dx4;
      T dy3 = dy4;
      T dz3 = dz0 - (SQUISH_CONSTANT * 3);
      T dw3 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 3);
      T attn3 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3) + (dw3 * dw3);
      value += std::pow(std::max((T)2.0 - attn3, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb, wsb + 1, dx3, dy3, dz3, dw3);

      // Contribution (1,0,1,1).
      {
        T dx2 = dx4;
        T dy2 = dy0 - (SQUISH_CONSTANT * 3);
        T dz2 = dz4;
        T dw2 = dw3;
        T attn2 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2) + (dw2 * dw2);
        value += std::pow(std::max((T)2.0 - attn2, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb + 1, wsb + 1, dx2, dy2, dz2, dw2);
      }

      // Contribution (0,1,1,1).
      {
        T dx1 = dx0 - (SQUISH_CONSTANT * 3);
        T dz1 = dz4;
        T dy1 = dy4;
        T dw1 = dw3;
        T attn1 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1) + (dw1 * dw1);
        value += std::pow(std::max((T)2.0 - attn1, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb + 1, wsb + 1, dx1, dy1, dz1, dw1);
      }

      // Contribution (1,1,0,0).
      {
        T dx5 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dy5 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dz5 = dz0 - (SQUISH_CONSTANT * 2);
        T dw5 = dw0 - (SQUISH_CONSTANT * 2);
        T attn5 = (dx5 * dx5) + (dy5 * dy5) + (dz5 * dz5) + (dw5 * dw5);
        value += std::pow(std::max((T)2.0 - attn5, (T)0.0), 4) * extrapolate(xsb + 1, ysb + 1, zsb, wsb, dx5, dy5, dz5, dw5);
      }

      // Contribution (1,0,1,0).
      {
        T dx6 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dy6 = dy0 - (SQUISH_CONSTANT * 2);
        T dz6 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dw6 = dw0 - (SQUISH_CONSTANT * 2);
        T attn6 = (dx6 * dx6) + (dy6 * dy6) + (dz6 * dz6) + (dw6 * dw6);
        value += std::pow(std::max((T)2.0 - attn6, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb + 1, wsb, dx6, dy6, dz6, dw6);
      }

      // Contribution (1,0,0,1).
      {
        T dx7 = dx0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dy7 = dy0 - (SQUISH_CONSTANT * 2);
        T dz7 = dz0 - (SQUISH_CONSTANT * 2);
        T dw7 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T attn7 = (dx7 * dx7) + (dy7 * dy7) + (dz7 * dz7) + (dw7 * dw7);
        value += std::pow(std::max((T)2.0 - attn7, (T)0.0), 4) * extrapolate(xsb + 1, ysb, zsb, wsb + 1, dx7, dy7, dz7, dw7);
      }

      // Contribution (0,1,1,0).
      {
        T dx8 = dx0 - (SQUISH_CONSTANT * 2);
        T dy8 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dz8 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dw8 = dw0 - (SQUISH_CONSTANT * 2);
        T attn8 = (dx8 * dx8) + (dy8 * dy8) + (dz8 * dz8) + (dw8 * dw8);
        value += std::pow(std::max((T)2.0 - attn8, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb + 1, wsb, dx8, dy8, dz8, dw8);
      }

      // Contribution (0,1,0,1).
      {
        T dx9 = dx0 - (SQUISH_CONSTANT * 2);
        T dy9 = dy0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dz9 = dz0 - (SQUISH_CONSTANT * 2);
        T dw9 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T attn9 = (dx9 * dx9) + (dy9 * dy9) + (dz9 * dz9) + (dw9 * dw9);
        value += std::pow(std::max((T)2.0 - attn9, (T)0.0), 4) * extrapolate(xsb, ysb + 1, zsb, wsb + 1, dx9, dy9, dz9, dw9);
      }

      // Contribution (0,0,1,1).
      {
        T dx10 = dx0 - (SQUISH_CONSTANT * 2);
        T dy10 = dy0 - (SQUISH_CONSTANT * 2);
        T dz10 = dz0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T dw10 = dw0 - (T)1.0 - (SQUISH_CONSTANT * 2);
        T attn10 = (dx10 * dx10) + (dy10 * dy10) + (dz10 * dz10) + (dw10 * dw10);
        value += std::pow(std::max((T)2.0 - attn10, (T)0.0), 4) * extrapolate(xsb, ysb, zsb + 1, wsb + 1, dx10, dy10, dz10, dw10);
      }
    }

    // First extra vertex.
    T attn_ext0 = (dx_ext0 * dx_ext0) + (dy_ext0 * dy_ext0) + (dz_ext0 * dz_ext0) + (dw_ext0 * dw_ext0);
    value += std::pow(std::max((T)2.0 - attn_ext0, (T)0.0), 4) * extrapolate(xsv_ext0, ysv_ext0, zsv_ext0, wsv_ext0, dx_ext0, dy_ext0, dz_ext0, dw_ext0);

    // Second extra vertex.
    T attn_ext1 = (dx_ext1 * dx_ext1) + (dy_ext1 * dy_ext1) + (dz_ext1 * dz_ext1) + (dw_ext1 * dw_ext1);
    value += std::pow(std::max((T)2.0 - attn_ext1, (T)0.0), 4) * extrapolate(xsv_ext1, ysv_ext1, zsv_ext1, wsv_ext1, dx_ext1, dy_ext1, dz_ext1, dw_ext1);

    // Third extra vertex.
    T attn_ext2 = (dx_ext2 * dx_ext2) + (dy_ext2 * dy_ext2) + (dz_ext2 * dz_ext2) + (dw_ext2 * dw_ext2);
    value += std::pow(std::max((T)2.0 - attn_ext2, (T)0.0), 4) * extrapolate(xsv_ext2, ysv_ext2, zsv_ext2, wsv_ext2, dx_ext2, dy_ext2, dz_ext2, dw_ext2);

    return (value / NORM_CONSTANT);
  }

};


// Array of gradient values for 4D. They approximate the directions to the
// vertices of a disprismatotesseractihexadecachoron from its center, skewed so that the
// tetrahedral and cubic facets can be inscribed in spheres of the same radius.
// Gradient set 2014-10-06.
template <typename T>
const int Noise<4, T>::gradients [] = {
   3, 1, 1, 1,   1, 3, 1, 1,   1, 1, 3, 1,   1, 1, 1, 3,
  -3, 1, 1, 1,  -1, 3, 1, 1,  -1, 1, 3, 1,  -1, 1, 1, 3,
   3,-1, 1, 1,   1,-3, 1, 1,   1,-1, 3, 1,   1,-1, 1, 3,
  -3,-1, 1, 1,  -1,-3, 1, 1,  -1,-1, 3, 1,  -1,-1, 1, 3,
   3, 1,-1, 1,   1, 3,-1, 1,   1, 1,-3, 1,   1, 1,-1, 3,
  -3, 1,-1, 1,  -1, 3,-1, 1,  -1, 1,-3, 1,  -1, 1,-1, 3,
   3,-1,-1, 1,   1,-3,-1, 1,   1,-1,-3, 1,   1,-1,-1, 3,
  -3,-1,-1, 1,  -1,-3,-1, 1,  -1,-1,-3, 1,  -1,-1,-1, 3,
   3, 1, 1,-1,   1, 3, 1,-1,   1, 1, 3,-1,   1, 1, 1,-3,
  -3, 1, 1,-1,  -1, 3, 1,-1,  -1, 1, 3,-1,  -1, 1, 1,-3,
   3,-1, 1,-1,   1,-3, 1,-1,   1,-1, 3,-1,   1,-1, 1,-3,
  -3,-1, 1,-1,  -1,-3, 1,-1,  -1,-1, 3,-1,  -1,-1, 1,-3,
   3, 1,-1,-1,   1, 3,-1,-1,   1, 1,-3,-1,   1, 1,-1,-3,
  -3, 1,-1,-1,  -1, 3,-1,-1,  -1, 1,-3,-1,  -1, 1,-1,-3,
   3,-1,-1,-1,   1,-3,-1,-1,   1,-1,-3,-1,   1,-1,-1,-3,
  -3,-1,-1,-1,  -1,-3,-1,-1,  -1,-1,-3,-1,  -1,-1,-1,-3
};

}

#else
#pragma message("OpenSimplexNoise.hh included multiple times")
#endif

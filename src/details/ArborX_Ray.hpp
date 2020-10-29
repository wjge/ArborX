/****************************************************************************
 * Copyright (c) 2012-2020 by the ArborX authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the ArborX library. ArborX is                       *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
#ifndef ARBORX_RAY_HPP
#define ARBORX_RAY_HPP

#include <ArborX_Box.hpp>
#include <ArborX_DetailsKokkosExt.hpp>
#include <ArborX_Point.hpp>

#include <Kokkos_Macros.hpp>

#include <cassert>
#include <cmath>

namespace ArborX
{
struct Ray
{
  using Vector = Point; // will regret this later

  using Scalar = std::decay_t<decltype(std::declval<Vector>()[0])>;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr Ray() = default;

  KOKKOS_FUNCTION
  Ray(Point const &origin, Vector const &direction)
      : _origin(origin)
      , _direction(direction)
  {
    normalize(_direction);
    for (int d = 0; d < 3; ++d)
    {
      _invdir[d] = 1.0 / _direction[d];
    }
  }

  KOKKOS_FUNCTION
  static Scalar norm(Vector const &v)
  {
    Scalar sq{};
    for (int d = 0; d < 3; ++d)
      sq += v[d] * v[d];
    return std::sqrt(sq);
  }

  KOKKOS_FUNCTION static void normalize(Vector &v)
  {
    auto const magv = norm(v);
    assert(magv > 0);
    for (int d = 0; d < 3; ++d)
      v[d] /= magv;
  }

  KOKKOS_FUNCTION
  constexpr Point &origin() { return _origin; }

  KOKKOS_FUNCTION
  constexpr Point const &origin() const { return _origin; }

  KOKKOS_FUNCTION
  constexpr Vector &direction() { return _direction; }

  KOKKOS_FUNCTION
  constexpr Vector const &direction() const { return _direction; }

  KOKKOS_FUNCTION
  constexpr Vector &invdir() { return _invdir; }

  KOKKOS_FUNCTION
  constexpr Vector const &invdir() const { return _invdir; }

  Point _origin = {};
  Vector _direction = {{0.f, 0.f, 0.f}};
  Vector _invdir = {{0.f, 0.f, 0.f}};
};

namespace Details
{
// The ray-box intersection algorithm is based on Majercik, A., et al. 2018.
// Their 'efficient slag' algorithm checks the intersections both in front and
// behind the ray. The algorithm here checks the intersections in front of the
// ray. Another reference is from Williams, A., et al. 2005
// and the website (key word: A minimal ray-tracer: rendering simple shapes).
//
// The NaN values will appear when minCorner[d] or maxCorner[d] == origin[d] and
// inv_ray_dir[d]==inf or -inf. The algorithm here adds an extra check for the
// NaN values to make it more robust. The compiler needs to follow the IEEE
// standard.
//
KOKKOS_INLINE_FUNCTION
bool intersects(Ray const &ray, Box const &box)
{
  auto const &minCorner = box.minCorner();
  auto const &maxCorner = box.maxCorner();
  auto const &origin = ray.origin();
  auto const &inv_ray_dir = ray.invdir();

  float max_min;
  float min_max;

  int init = 0;

  for (int d = 0; d < 3; ++d)
  {
    float tmin;
    float tmax;

    if (inv_ray_dir[d] >= 0)
    {
      tmin = (minCorner[d] - origin[d]) * inv_ray_dir[d];
      tmax = (maxCorner[d] - origin[d]) * inv_ray_dir[d];
    }
    else
    {
      tmin = (maxCorner[d] - origin[d]) * inv_ray_dir[d];
      tmax = (minCorner[d] - origin[d]) * inv_ray_dir[d];
    }

    if (d == init)
    {
      if (tmin == tmin && tmax == tmax)
      {
        max_min = tmin;
        min_max = tmax;
      }
      else
      {
        init = d + 1;
        continue;
      }
    }
    else
    {
      max_min = max_min < tmin ? tmin : max_min;
      min_max = min_max > tmax ? tmax : min_max;
    }
  }
  return max_min <= min_max && (min_max >= 0);
}

} // namespace Details
} // namespace ArborX
#endif

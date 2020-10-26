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
      // direction can be zero
      _invdir[d] = 1.0 / (_direction[d]);
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
// behind the ray. The function here checks the intersections in front of the
// ray. The potential issue is the division of zero, when _direction[d] is
// zero. The IEEE standard guarantees the algorithm works for the
// infinities, which is discussed in more detail in Williams, A., et al. 2005
// and the website (key word: A minimal ray-tracer: rendering simple shapes).
//
// NaN values will appear when minCorner[d] or maxCorner[d] == origin[d] and
// inv_ray_dir[d]==inf or -inf. The current implementation assumes any
// comparison "> <" with NaN returns false.
KOKKOS_INLINE_FUNCTION
static bool intersects(Ray const &ray, Box const &box)
{
  auto const &minCorner = box.minCorner();
  auto const &maxCorner = box.maxCorner();
  auto const &origin = ray.origin();
  auto const &inv_ray_dir = ray.invdir();

  float max_min;
  float min_max;

  for (int d = 0; d < 3; ++d)
  {
    float tmin;
    float tmax;
    // Still not sure how robust this is, as it deals with nan and inf. For
    // example, replacing if() with
    //
    //     t0 = (minCorner[d] - origin[d]) * inv_ray_dir[d];
    //     t1 = (maxCorner[d] - origin[d]) * inv_ray_dir[d];
    //     tmin = min(t0, t1);
    //     tmax = max(t0, t1);
    //
    // does not work.
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

    if (d == 0)
    {
      max_min = tmin;
      min_max = tmax;
    }
    else
    {
      // max_min = KokkosExt::max(max_min, tmin);
      // min_max = KokkosExt::min(min_max, tmax);
      if (max_min < tmin)
        max_min = tmin;
      if (min_max > tmax)
        min_max = tmax;
    }
  }

  return max_min <= min_max && (min_max >= 0);
}
} // namespace Details

} // namespace ArborX

#endif

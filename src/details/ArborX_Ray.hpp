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

  class Normalized
  {
    Vector v_;

  public:
    KOKKOS_FUNCTION
    explicit Normalized(Vector const &v)
        : v_{v}
    {
    }

    KOKKOS_FUNCTION
    operator Vector &() { return v_; }

    KOKKOS_FUNCTION
    operator Vector const &() const { return v_; }
  };

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
  Ray(Point const &origin, Normalized const &direction)
      : _origin(origin)
      , _direction(direction)
  {
  }

  KOKKOS_FUNCTION
  static constexpr Scalar norm_sq(Vector const &v)
  {
    Scalar sq{};
    for (int d = 0; d < 3; ++d)
      sq += v[d] * v[d];
    return sq;
  }

  KOKKOS_FUNCTION
  static constexpr Scalar norm(Vector const &v)
  {
    return std::sqrt(norm_sq(v));
  }

  KOKKOS_FUNCTION static void normalize(Vector &v)
  {
    auto const n = norm(v);
    assert(n > 0);
    for (int d = 0; d < 3; ++d)
      v[d] /= n;
  }

  KOKKOS_FUNCTION
  constexpr Point &origin() { return _origin; }

  KOKKOS_FUNCTION
  constexpr Point const &origin() const { return _origin; }

  KOKKOS_FUNCTION
  constexpr Point &direction() { return _direction; }

  KOKKOS_FUNCTION
  constexpr Point const &direction() const { return _direction; }

  KOKKOS_FUNCTION
  constexpr Point &invdir() { return _invdir; }

  KOKKOS_FUNCTION
  constexpr Point const &invdir() const { return _invdir; }

  Point _origin = {};
  Point _direction = {{0.f, 0.f, 0.f}};
  Point _invdir = {{0.f, 0.f, 0.f}};
};

namespace Details
{
namespace RayTracing
{
// ray-box intersection
// The ray-box intersection algorithm is based on Majercik, A., et al. 2018.
// Their 'efficient slag' algorithm checks the intersections both in front and
// behind the ray. The function here checks the intersections in front of the
// ray. The potential issue is the division of zero, when _direction[d] is
// zero. The IEEE standard guarantees the algorithm works for the
// infinities, which is discussed in more detail in Williams, A., et al. 2005
// and the website (key word: A minimal ray-tracer: rendering simple shapes).
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
    float tmin, tmax;
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
      max_min = max(max_min, tmin);
      min_max = min(min_max, tmax);
    }
  }

  return max_min <= min_max && (min_max >= 0);
}
} // namespace RayTracing
} // namespace Details

} // namespace ArborX

#endif

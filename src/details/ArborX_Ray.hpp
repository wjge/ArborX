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

  KOKKOS_FUNCTION static void normalize (Vector &v)
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

  Point _origin = {};
  Point _direction = {{0.f, 0.f, 0.f}};
};
} // namespace ArborX

#endif

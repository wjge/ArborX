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

namespace ArborX
{
struct Ray
{
  KOKKOS_DEFAULTED_FUNCTION
  constexpr Ray() = default;

  KOKKOS_INLINE_FUNCTION
  constexpr Ray(Point const &origin, Point const &direction)
      : _origin(origin)
      , _direction(direction)
  {
  }

  KOKKOS_INLINE_FUNCTION
  constexpr Point &origin() { return _origin; }

  KOKKOS_INLINE_FUNCTION
  constexpr Point const &origin() const { return _origin; }

  KOKKOS_INLINE_FUNCTION
  constexpr Point &direction() { return _direction; }

  KOKKOS_INLINE_FUNCTION
  constexpr Point const &direction() const { return _direction; }

  Point _origin = {};
  Point _direction = {{0.f, 0.f, 0.f}};
};
} // namespace ArborX

#endif

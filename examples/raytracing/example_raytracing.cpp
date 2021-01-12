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

#include <ArborX.hpp>

#include <Kokkos_Core.hpp>

#include <iostream>
#include <random>
#include <vector>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = ExecutionSpace::memory_space;

namespace ArborX
{
using Primitives = Kokkos::View<Box *, MemorySpace>;
using Predicates = Kokkos::View<Ray *, MemorySpace>;

template <>
struct AccessTraits<Predicates, PredicatesTag>
{
  KOKKOS_FUNCTION static std::size_t size(const Predicates &predicates)
  {
    return predicates.extent(0);
  }
  KOKKOS_FUNCTION static auto get(Predicates const &predicates, std::size_t i)
  {
    return intersects(predicates(i));
  }
  using memory_space = MemorySpace;
};
} // namespace ArborX

struct Callback
{
  ArborX::Primitives boxes;

  KOKKOS_FUNCTION float dotproduct(ArborX::Point const &v1,
                                   ArborX::Point const &v2) const
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  KOKKOS_FUNCTION float solveQuadratic(const float &b, const float &c) const
  {
    float delta = b * b - 4.0 * c;
    if (delta <= 0.0)
      return 0.0;
    else
    {
      float q = (b > 0) ? -0.5 * (b + std::sqrt(delta))
                        : -0.5 * (b - std::sqrt(delta));
      return std::abs(q - c / q);
    }
  }

  KOKKOS_FUNCTION float overlap(ArborX::Ray const &ray,
                                ArborX::Box const &box) const
  {
    ArborX::Point mincorner = box.minCorner();
    ArborX::Point maxcorner = box.maxCorner();

    ArborX::Point center{(float)0.5 * (mincorner[0] + maxcorner[0]),
                         (float)0.5 * (mincorner[1] + maxcorner[1]),
                         (float)0.5 * (mincorner[2] + maxcorner[2])};

    float r = 0.5 * (maxcorner[0] - mincorner[0]);

    ArborX::Point orig = ray.origin();

    ArborX::Point L{orig[0] - center[0], orig[1] - center[1],
                    orig[2] - center[2]};

    //  for normalized direction vector a = 1.0
    auto dir = ray.direction();

    float b = 2.0 * dotproduct(dir, L);
    float c = dotproduct(L, L) - r * r;

    return solveQuadratic(b, c);
  }

  template <typename Predicate, typename OutputFunctor>
  KOKKOS_FUNCTION void operator()(Predicate predicate, int primitive,
                                  OutputFunctor const &out) const
  {
    auto const &ray = ArborX::getGeometry(predicate);
    auto const &box = boxes(primitive);
    float length = overlap(ray, box);

#ifndef __SYCL_DEVICE_ONLY__
    printf("ray hit box %d with the overlap %f. \n", primitive, length);
#endif
    out(primitive);
  }
};

int main(int argc, char *argv[])
{
  Kokkos::ScopeGuard guard(argc, argv);

  int const n = 10;
  std::vector<ArborX::Box> boxes;
  // Fill vector with random points in [-1, 1]^3
  std::uniform_real_distribution<float> dis{-1., 1.};
  std::default_random_engine gen;
  auto rd = [&]() { return dis(gen); };

  std::generate_n(std::back_inserter(boxes), n, [&]() {
    float x0 = rd();
    float y0 = rd();
    float z0 = rd();

    float x1 = x0 + 1.0;
    float y1 = y0 + 1.0;
    float z1 = z0 + 1.0;

    return ArborX::Box{{x0, y0, z0}, {x1, y1, z1}};
  });

  ArborX::Primitives boxes_view("boxes_view", n);
  auto boxes_host = Kokkos::create_mirror_view(boxes_view);
  for (int i = 0; i < n; ++i)
  {
    boxes_host[i] = boxes[i];
  }
  Kokkos::deep_copy(boxes_view, boxes_host);
  ArborX::BVH<MemorySpace> bvh{ExecutionSpace{}, boxes_view};

  int const m = 1;
  ArborX::Predicates rays_view("rays_view", m);

  {
    Kokkos::View<int *, MemorySpace> values("values", 0);
    Kokkos::View<int *, MemorySpace> offsets("offsets", 0);
    ArborX::query(bvh, ExecutionSpace{}, rays_view, Callback{boxes_view},
                  values, offsets);
  }
  return 0;
}

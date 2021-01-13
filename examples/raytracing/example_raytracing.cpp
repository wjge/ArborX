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
#include <Kokkos_ScatterView.hpp>

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

  Kokkos::View<float *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Atomic>>
      accumulator;

  KOKKOS_FUNCTION float dotproduct(ArborX::Point const &v1,
                                   ArborX::Point const &v2) const
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  KOKKOS_FUNCTION float solveQuadratic(const float &b, const float &c) const
  {
    float delta = b * b - 4.0 * c;
    if (delta <= 0.0)
      // NO intersection
      return 0.0;
    else
    {
      float q = (b > 0) ? -0.5 * (b + std::sqrt(delta))
                        : -0.5 * (b - std::sqrt(delta));

      // printf("q=%f, c/q = %f \n", q, c / q);

      if (q < 0)
      {
        // both intersections are BEHIND
        if (c > 0)
          return 0.0;
        // one point is behind
        else
          return c / q;
      }
      else
      {
        // both intersections are AHEAD
        if (c > 0)
          return std::abs(q - c / q);
        // one point is behind
        else
          return q;
      }
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

  template <typename Predicate>
  KOKKOS_FUNCTION void operator()(Predicate const &predicate,
                                  int primitive) const
  {
    auto const &ray = ArborX::getGeometry(predicate);
    auto const &box = boxes(primitive);
    float length = overlap(ray, box);
    // auto i = getData(predicate);
#ifndef __SYCL_DEVICE_ONLY__
    printf("ray hit box %d with the overlap %f. \n", primitive, length);
#endif
  }
};

int main(int argc, char *argv[])
{
  Kokkos::ScopeGuard guard(argc, argv);

  // size of the domain (LxLxL)
  const float L = 100.0;

  // number of primitives
  int const n = 100;
  std::vector<ArborX::Box> boxes;
  std::uniform_real_distribution<float> uniform{0.0, 1.0};
  std::default_random_engine gen;
  auto rd_uniform = [&]() { return uniform(gen); };

  // radius
  const float mu_R = 1.0;
  const float sigma_R = mu_R / 3.0;

  std::normal_distribution<> normal{mu_R, sigma_R};
  auto rd_normal = [&]() { return std::max(normal(gen), 0.0); };

  std::generate_n(std::back_inserter(boxes), n, [&]() {
    float r = rd_normal();

    float x0 = rd_uniform() * L - r;
    float y0 = rd_uniform() * L - r;
    float z0 = rd_uniform() * L - r;

    float x1 = x0 + 2.0 * r;
    float y1 = y0 + 2.0 * r;
    float z1 = z0 + 2.0 * r;

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

  // initialize predicates
  int const m = 10000;
  float const PI = 4.0 * std::atan(1.0);

  std::vector<ArborX::Ray> rays;
  std::generate_n(std::back_inserter(rays), m, [&]() {
    float x = rd_uniform() * L;
    float y = rd_uniform() * L;
    float z = 0.0;

    float sinazimuth = std::sin(rd_uniform() * PI * 2.0);
    float cosazimuth = std::sqrt(1 - sinazimuth * sinazimuth);
    float cospolar = 1.0 - rd_uniform();
    float sinpolar = std::sqrt(1.0 - cospolar * cospolar);

    float dirx = sinpolar * cosazimuth;
    float diry = sinpolar * sinazimuth;
    float dirz = cospolar;

    return ArborX::Ray{{x, y, z}, {dirx, diry, dirz}};
  });

  ArborX::Predicates rays_view("rays_view", m);
  auto rays_host = Kokkos::create_mirror_view(rays_view);
  for (int i = 0; i < m; ++i)
  {
    rays_host[i] = rays[i];
  }

  Kokkos::deep_copy(rays_view, rays_host);

  {
    // accumulator
    Kokkos::View<float[m], MemorySpace> overlap_view("overlap_view");

    bvh.query(ExecutionSpace{}, rays_view, Callback{boxes_view, overlap_view});
  }
  return 0;
}

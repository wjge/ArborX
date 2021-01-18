/****************************************************************************
 * Copyright (c) 2017-2021 by the ArborX authors                            *
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

template <typename MemorySpace>
struct ArborX::AccessTraits<Kokkos::View<ArborX::Sphere *, MemorySpace>,
                            ArborX::PrimitivesTag>
{
  using Primitives = Kokkos::View<ArborX::Sphere *, MemorySpace>;
  KOKKOS_FUNCTION static std::size_t size(const Primitives &primitives)
  {
    return primitives.extent(0);
  }
  KOKKOS_FUNCTION static ArborX::Box get(Primitives const &primitives,
                                         std::size_t i)
  {
    auto const &sphere = primitives(i);
    auto const &c = sphere.centroid();
    auto const &r = sphere.radious();
    return {{c[0] - r, c[1] - r, c[2] - r}, {c[0] + r, c[1] + r, c[2] + r}};
  }
  using memory_space = MemorySpace;
};

template <typename MemorySpace>
struct ArborX::AccessTraits<Kokkos::View<ArborX::Ray *, MemorySpace>,
                            ArborX::PredicatesTag>
{
  using Predicates = Kokkos::View<ArborX::Ray *, MemorySpace>;
  KOKKOS_FUNCTION static std::size_t size(const Predicates &predicates)
  {
    return predicates.extent(0);
  }
  KOKKOS_FUNCTION static auto get(Predicates const &predicates, std::size_t i)
  {
    return attach(intersects(predicates(i)), (int)i);
  }
  using memory_space = MemorySpace;
};

template <typename MemorySpace>
struct Callback
{
  using Primitives = Kokkos::View<ArborX::Sphere *, MemorySpace>;
  using Vector = ArborX::Point;

  Primitives spheres;

  Kokkos::View<float *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Atomic>>
      accumulator;

  KOKKOS_FUNCTION float dotproduct(Vector const &v1, Vector const &v2) const
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

  KOKKOS_FUNCTION float overlap_distance(ArborX::Ray const &ray,
                                         ArborX::Sphere const &sphere) const
  {
    auto center = sphere.centroid();
    auto r = sphere.radius();

    ArborX::Point orig = ray.origin();

    Vector L{orig[0] - center[0], orig[1] - center[1], orig[2] - center[2]};

    //  for normalized direction vector a = 1.0
    auto dir = ray.direction();

    float b = 2.0 * dotproduct(dir, L);
    float c = dotproduct(L, L) - r * r;

    return solveQuadratic(b, c);
  }

  template <typename Predicate>
  KOKKOS_FUNCTION void operator()(Predicate const &predicate,
                                  int primitive_index) const
  {
    auto const &ray = ArborX::getGeometry(predicate);
    auto const &sphere = spheres(primitive_index);
    float length = overlap_distance(ray, sphere);
    int i = getData(predicate);
    accumulator(i) += length;
  }
};

int main(int argc, char *argv[])
{
  using ExecutionSpace = Kokkos::DefaultExecutionSpace;
  using MemorySpace = ExecutionSpace::memory_space;

  using Primitives = Kokkos::View<ArborX::Sphere *, MemorySpace>;
  using Predicates = Kokkos::View<ArborX::Ray *, MemorySpace>;

  Kokkos::ScopeGuard guard(argc, argv);

  // size of the domain (LxLxL)
  const float L = 100.0;

  // number of primitives
  int const n = 100;
  std::vector<ArborX::Sphere> spheres;
  std::uniform_real_distribution<float> uniform{0.0, 1.0};
  std::default_random_engine gen;
  auto rd_uniform = [&]() { return uniform(gen); };

  // radius
  const float mu_R = 1.0;
  const float sigma_R = mu_R / 3.0;

  std::normal_distribution<> normal{mu_R, sigma_R};
  auto rd_normal = [&]() { return std::max(normal(gen), 0.0); };

  Primitives spheres_view("spheres_view", n);
  auto spheres_host = Kokkos::create_mirror_view(spheres_view);
  for (int i = 0; i < n; ++i)
  {
    spheres_host[i].centroid() = {rd_uniform(), rd_uniform(), rd_uniform()};
    spheres_host[i].radius() = rd_normal();
  }
  Kokkos::deep_copy(spheres_view, spheres_host);
  ArborX::BVH<MemorySpace> bvh{ExecutionSpace{}, spheres_view};

  // initialize predicates
  int const m = 10000;

  std::vector<ArborX::Ray> rays;
  std::generate_n(std::back_inserter(rays), m, [&]() {
    float x = rd_uniform() * L;
    float y = rd_uniform() * L;
    float z = 0.0;

    float sinazimuth = std::sin(rd_uniform() * 2.0 * M_PI);
    float cosazimuth = std::sqrt(1 - sinazimuth * sinazimuth);
    float cospolar = 1.0 - rd_uniform();
    float sinpolar = std::sqrt(1.0 - cospolar * cospolar);

    float dirx = sinpolar * cosazimuth;
    float diry = sinpolar * sinazimuth;
    float dirz = cospolar;

    return ArborX::Ray{{x, y, z}, {dirx, diry, dirz}};
  });

  Predicates rays_view("rays_view", m);
  auto rays_host = Kokkos::create_mirror_view(rays_view);
  for (int i = 0; i < m; ++i)
  {
    rays_host[i] = rays[i];
  }

  Kokkos::deep_copy(rays_view, rays_host);

  {
    // accumulator
    Kokkos::View<float *, MemorySpace> overlap_view("overlap_view", m);

    bvh.query(ExecutionSpace{}, rays_view,
              Callback<MemorySpace>{spheres_view, overlap_view});
  }
  return 0;
}

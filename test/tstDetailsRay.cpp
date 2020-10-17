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
#include <ArborX_Box.hpp>
#include <ArborX_DetailsAlgorithms.hpp>
#include <ArborX_DetailsKokkosExt.hpp> // sgn
#include <ArborX_Ray.hpp>

#include <boost/test/unit_test.hpp>

#define BOOST_TEST_MODULE DetailsRay

namespace ArborX
{
namespace Details
{
Point &operator+=(Point &rhs, Point const &lhs)
{
  for (int d = 0; d < 3; ++d)
    rhs[d] += lhs[d];
  return rhs;
}
Point operator*(Ray::Scalar a, Point p)
{
  for (int d = 0; d < 3; ++d)
    p[d] *= a;
  return p;
}

/*old algorithm by dalg
//KOKKOS_INLINE_FUNCTION
bool intersects(Ray const &r, Box const &b)
{
  using KokkosExt::min;
  //using KokkosExt::sgn;
  for (int d = 0; d < 3; ++d)
  {
    auto const sign = sgn(r.direction()[d]);
    if (sign == 0)
      continue;
    for (auto proj : {b.minCorner()[d], b.maxCorner()[d]})
    {
      auto dot = (proj - r.origin()[d]);
      dot *= r.direction()[d];
      if (dot < 0)
        continue;
      auto test{r.origin()};
      test += dot * r.direction();
      test[d] = proj;
      if (intersects(test, b))
        return true;
    }
  }
  return false;
}
*/

KOKKOS_FUNCTION
bool intersects(Ray const &ray, Box const &box)
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

} // namespace Details
} // namespace ArborX

BOOST_AUTO_TEST_CASE(intersects_box)
{
  using ArborX::Box;
  using ArborX::Point;
  using ArborX::Ray;
  using ArborX::Details::intersects;

  Box unit_box1{{0, 0, 0}, {1, 1, 1}};
  Box unit_box2{{-1, -1, -1}, {0, 0, 0}};

  // box1
  // origin is within the box
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {0, 1, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {0, 0, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {1, 1, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {1, 1, 1}}, unit_box1));

  // origin is outside the box
  // hit the center of the face
  BOOST_TEST(intersects(Ray{{-1, .5, .5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, .5, .5}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, .5, .5}, {0, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, .5, .5}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, .5, .5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, .5, .5}, {0, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{-1, 1.5, .5}, {1, -1, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, .5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{-1, 1.5, 1.5}, {1, -1, -1}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, 1.5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, 1.5}, {1, -1, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, 1.5}, {1, 0, -1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, -.5, .5}, {-1, 1, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, -.5, .5}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, -.5, 1.5}, {-1, 1, -1}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, -.5, 1.5}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, -.5, 1.5}, {-1, 1, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, -.5, 1.5}, {-1, 0, -1}}, unit_box1));
  // hit the center of an edge
  BOOST_TEST(intersects(Ray{{-1, 0, .5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 0, .5}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 0, .5}, {0, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, 0, .5}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, 0, .5}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, 0, .5}, {0, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{-1, -1, .5}, {1, 1, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{-1, 1, .5}, {1, -1, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{-1, -2, .5}, {1, 2, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, 2, .5}, {-1, -1, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, -1, .5}, {-1, 1, 0}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, 1, .5}, {-1, -1, 0}}, unit_box1));
  // hit a corner
  BOOST_TEST(intersects(Ray{{-1, 1, 1}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 1, 1}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{-1, 1, 1}, {0, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, 1, 1}, {-1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, 1, 1}, {1, 0, 0}}, unit_box1));
  BOOST_TEST(!intersects(Ray{{2, 1, 1}, {0, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{-1, -1, -1}, {1, 1, 1}}, unit_box1));
  BOOST_TEST(intersects(Ray{{2, 2, 2}, {-1, -1, -1}}, unit_box1));

  BOOST_TEST(!intersects(Ray{{1, 2, 3}, {4, 5, 6}}, unit_box1));
  BOOST_TEST(intersects(Ray{{1, 2, 3}, {-1, -2, -3}}, unit_box1));

  // box2
  // origin is within the box
  BOOST_TEST(intersects(Ray{{-.5, -.5, -.5}, {1, 0, 0}}, unit_box2));
  BOOST_TEST(intersects(Ray{{-.5, -.5, -.5}, {0, 1, 0}}, unit_box2));
  BOOST_TEST(intersects(Ray{{-.5, -.5, -.5}, {0, 0, 1}}, unit_box2));

  // origin is outside the box
  // hit the center of the face
  BOOST_TEST(intersects(Ray{{1, -.5, -.5}, {-1, 0, 0}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{1, -.5, -.5}, {1, 0, 0}}, unit_box2));
  BOOST_TEST(intersects(Ray{{-.5, 1, -.5}, {0, -1, 0}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{-.5, 1, -.5}, {0, 1, 0}}, unit_box2));
  BOOST_TEST(intersects(Ray{{-.5, -.5, 1}, {0, 0, -1}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{-.5, -.5, 1}, {0, 0, 1}}, unit_box2));

  // hit the center of an edge
  BOOST_TEST(intersects(Ray{{.5, .5, -1}, {-1, -1, 0}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{.5, .5, -1}, {1, 1, 0}}, unit_box2));
  BOOST_TEST(intersects(Ray{{.5, -1, .5}, {-1, 0, -1}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{.5, -1, .5}, {1, 0, 1}}, unit_box2));
  BOOST_TEST(intersects(Ray{{-1, .5, .5}, {0, -1, -1}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{-1, .5, .5}, {0, 1, 1}}, unit_box2));

  // hit a corner
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {-1, -1, -1}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{.5, .5, .5}, {1, 1, 1}}, unit_box2));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {-1, -1, -1}}, unit_box2));
  BOOST_TEST(!intersects(Ray{{.5, .5, .5}, {1, 1, 1}}, unit_box2));
}

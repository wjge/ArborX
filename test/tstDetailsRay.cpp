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
#include <ArborX_Ray.hpp>

#include <boost/test/unit_test.hpp>

#define BOOST_TEST_MODULE DetailsRay

BOOST_AUTO_TEST_CASE(intersects_box)
{
  using ArborX::Box;
  using ArborX::Point;
  using ArborX::Ray;
  constexpr auto &intersects = ArborX::Details::RayTracing::intersects;

  float pinf = 1.0 / 0.0;
  float ninf = -1.0 / 0.0;
  float zero_pinf = 0 * pinf;
  float zero_ninf = 0 * ninf;
  BOOST_TEST(!(zero_pinf == zero_pinf));
  BOOST_TEST(!(zero_ninf == zero_ninf));
  BOOST_TEST(!(ninf < zero_pinf));
  BOOST_TEST(!(pinf < zero_pinf));
  BOOST_TEST(!(ninf > zero_ninf));
  BOOST_TEST(!(pinf > zero_ninf));
  BOOST_TEST(!(1.0 < zero_pinf));
  BOOST_TEST(!(1.0 > zero_pinf));
  BOOST_TEST(!(-1.0 > zero_ninf));
  BOOST_TEST(!(-1.0 > zero_ninf));

  Box unit_box{{0, 0, 0}, {1, 1, 1}};

  // origin is within the box
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {0, 1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {0, 0, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {1, 1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{.5, .5, .5}, {1, 1, 1}}, unit_box));

  // origin is outside the box
  // hit the center of the face
  BOOST_TEST(intersects(Ray{{-1, .5, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, .5, .5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, .5, .5}, {0, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, .5, .5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, .5, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, .5, .5}, {0, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{1, .5, .5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{1, .5, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{.5, 1, .5}, {0, -1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{.5, 1, .5}, {0, 1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{.5, .5, 1}, {0, 0, -1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{.5, .5, 1}, {0, 0, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{-1, 1.5, .5}, {1, -1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{-1, 1.5, 1.5}, {1, -1, -1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, 1.5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, 1.5}, {1, -1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 1.5, 1.5}, {1, 0, -1}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, -.5, .5}, {-1, 1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, -.5, .5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, -.5, 1.5}, {-1, 1, -1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, -.5, 1.5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, -.5, 1.5}, {-1, 1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, -.5, 1.5}, {-1, 0, -1}}, unit_box));
  // hit the center of an edge
  BOOST_TEST(intersects(Ray{{-1, 0, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 0, .5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 0, .5}, {0, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, 0, .5}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, 0, .5}, {1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, 0, .5}, {0, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{-1, -1, .5}, {1, 1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{-1, 1, .5}, {1, -1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{-1, -2, .5}, {1, 2, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, 2, .5}, {-1, -1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, -1, .5}, {-1, 1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, 1, .5}, {-1, -1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{0, 0, .5}, {1, 1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{0, 0, .5}, {-1, -1, 0}}, unit_box));
  BOOST_TEST(intersects(Ray{{1, 1, .5}, {-1, -1, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{1, 1, .5}, {1, 1, 0}}, unit_box));

  // hit a corner
  BOOST_TEST(intersects(Ray{{-1, 1.5, 1.5}, {1, -1, -1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 1, 1}, {-1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{-1, 1, 1}, {0, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{1, 1, 1}, {-1, -1, -1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{1, 1, 1}, {1, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{0, 0, 0}, {1, 1, 1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{0, 0, 0}, {-1, -1, -1}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, 1, 1}, {1, 0, 0}}, unit_box));
  BOOST_TEST(!intersects(Ray{{2, 1, 1}, {0, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{-1, -1, -1}, {1, 1, 1}}, unit_box));
  BOOST_TEST(intersects(Ray{{2, 2, 2}, {-1, -1, -1}}, unit_box));

  BOOST_TEST(!intersects(Ray{{1, 2, 3}, {4, 5, 6}}, unit_box));
  BOOST_TEST(intersects(Ray{{1, 2, 3}, {-1, -2, -3}}, unit_box));
}

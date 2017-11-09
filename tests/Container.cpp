/** @file
    @brief Tests

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>

    SPDX-License-Identifier:BSD-3-Clause
*/

// Copyright 2017 Sensics, Inc.

#include <subdiv2d/SubdivContainer.h>

#include "catch.hpp"

using namespace sensics::subdiv2d;
using SubdivDoubleContainer = SubdivContainer<double>;
TEST_CASE("Container constructor behavior", "[SubdivContainer]") {
    REQUIRE_NOTHROW(SubdivDoubleContainer(Rect(0, 0, 1, 1)));
    SubdivDoubleContainer subdiv(Rect(0, 0, 1, 1));
    THEN("Should be able to insert a point.") { REQUIRE_NOTHROW(subdiv.insert(Point2f(0, 0), 1.0)); }
}

TEST_CASE("Simple retrieval", "[SubdivContainer]") {
    SubdivDoubleContainer subdiv(Rect(0, 0, 1, 1));

    const Point2f ZeroZero(0, 0);

    REQUIRE(!subdiv.lookup(ZeroZero));

    subdiv.insert(ZeroZero, 1.0);

    REQUIRE(subdiv.lookup(ZeroZero));
    REQUIRE(1.0 == subdiv.get(ZeroZero));
}

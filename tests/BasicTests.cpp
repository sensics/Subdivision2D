/** @file
    @brief Tests

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>

    SPDX-License-Identifier:BSD-3-Clause
*/

// Copyright 2017 Sensics, Inc.

#include <subdiv2d/Subdivision2D.h>

#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

using namespace sensics::subdiv2d;
TEST_CASE("Default constructor behavior", "[Subdivision2d]") {
    REQUIRE_NOTHROW(Subdiv2D());
    Subdiv2D subdiv;
    THEN("Should not be able to insert a point with un-initialized subdivision structure") {
        REQUIRE_THROWS(subdiv.insert(Point2f(0, 0)));
    }
    {
        auto edges = subdiv.getEdgeList();
        std::cerr << "Now have " << edges.size() << " edges" << std::endl;
    }
    GIVEN("uninitialized object with initDelauney called on it") {
        REQUIRE_NOTHROW(subdiv.initDelaunay(Rect(0, 0, 1, 1)));
        {
            auto edges = subdiv.getEdgeList();
            std::cerr << "Now have " << edges.size() << " edges" << std::endl;
        }
        THEN("Should be able to insert a point.") { REQUIRE_NOTHROW(subdiv.insert(Point2f(0, 0))); }
    }
}
TEST_CASE("Constructor from Rect behavior", "[Subdivision2d]") {

    REQUIRE_NOTHROW(Subdiv2D(Rect(0, 0, 1, 1)));
    Subdiv2D subdiv(Rect(0, 0, 1, 1));
    {
        auto edges = subdiv.getEdgeList();
        std::cerr << "Now have " << edges.size() << " edges" << std::endl;
    }
    THEN("Should be able to insert a point immediately") { REQUIRE_NOTHROW(subdiv.insert(Point2f(0, 0))); }
}

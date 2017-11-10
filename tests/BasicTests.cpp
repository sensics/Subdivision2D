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

TEST_CASE("Locate vertices", "[Subdivision2d]") {

    using Catch::Matchers::Equals;
    using Catch::Matchers::VectorContains;
    Subdiv2D subdiv(Rect(0, 0, 5, 5));
    const auto OneOne = Point2f(1, 1);
    const auto FourOne = Point2f(4, 1);
    const auto OneFour = Point2f(1, 4);
    subdiv.insert(OneOne);
    subdiv.insert(FourOne);
    subdiv.insert(OneFour);
    WHEN("Looking up any of the input points") {
        THEN("we should get exactly one vertex") {
            REQUIRE(1 == subdiv.locateVertexIds(OneOne).size());
            REQUIRE(1 == subdiv.locateVertexIds(FourOne).size());
            REQUIRE(1 == subdiv.locateVertexIds(OneFour).size());
        }
        THEN("we should get the vertex we look up") {
            REQUIRE_THAT(subdiv.locateVertices(OneOne), VectorContains(OneOne));
            REQUIRE_THAT(subdiv.locateVertices(FourOne), VectorContains(FourOne));
            REQUIRE_THAT(subdiv.locateVertices(OneFour), VectorContains(OneFour));
        }
    }

    WHEN("Looking up a point along an edge") {
        THEN("we should get two vertices") {
            REQUIRE(2 == subdiv.locateVertexIds(Point2f(1, 2)).size());
            REQUIRE(2 == subdiv.locateVertexIds(Point2f(2, 1)).size());
        }
        THEN("the vertices returned should be the expected given ones") {
            // REQUIRE_THAT(subdiv.locateVertices(Point2f(1, 2)), Equals(std::vector<Point2f>({OneOne, OneFour})));
            REQUIRE_THAT(subdiv.locateVertices(Point2f(1, 2)), VectorContains(OneOne));
            REQUIRE_THAT(subdiv.locateVertices(Point2f(1, 2)), VectorContains(OneFour));

            // REQUIRE_THAT(subdiv.locateVertices(Point2f(2, 1)), Equals(std::vector<Point2f>({OneOne, FourOne})));
            REQUIRE_THAT(subdiv.locateVertices(Point2f(2, 1)), VectorContains(OneOne));
            REQUIRE_THAT(subdiv.locateVertices(Point2f(2, 1)), VectorContains(FourOne));
        }
    }

    WHEN("Looking up an interior point") {
        THEN("we should get three points") { REQUIRE(3 == subdiv.locateVertexIds(Point2f(2, 2)).size()); }
        THEN("we should get all three given points") {
            REQUIRE_THAT(subdiv.locateVertices(Point2f(2, 2)), VectorContains(OneOne));
            REQUIRE_THAT(subdiv.locateVertices(Point2f(2, 2)), VectorContains(FourOne));
            REQUIRE_THAT(subdiv.locateVertices(Point2f(2, 2)), VectorContains(OneFour));
        }
    }
}

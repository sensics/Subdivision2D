/** @file
    @brief Implementation

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
    GIVEN("uninitialized object with initDelauney called on it") {
        REQUIRE_NOTHROW(subdiv.initDelaunay(Rect(0, 0, 1, 1)));
        THEN("Should be able to insert a point.") { REQUIRE_NOTHROW(subdiv.insert(Point2f(0, 0))); }
    }
}
TEST_CASE("Constructor from Rect behavior", "[Subdivision2d]") {

    REQUIRE_NOTHROW(Subdiv2D(Rect(0, 0, 1, 1)));
    Subdiv2D subdiv(Rect(0, 0, 1, 1));
    THEN("Should be able to insert a point immediately") { REQUIRE_NOTHROW(subdiv.insert(Point2f(0, 0))); }
}

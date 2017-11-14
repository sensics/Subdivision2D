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

// Internal Includes
#include "AngleMeshUtils.h"
#include "GenericExtremaFinder.h"
#include <subdiv2d/SubdivContainer.h>

// Library/third-party includes
#include <Eigen/Core>

#include "EigenStdArrayInterop.h"

// Standard includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

static const float STEPS = 5;

using namespace sensics::subdiv2d;

inline bool orderAfromBwrtCenter(Point2f const& a, Point2f const& b, Point2f const& center) {
    // based on https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
    const auto axdiff = a.x - center.x;
    const auto bxdiff = b.x - center.x;
    if (axdiff >= 0 && bxdiff < 0) {
        return true;
    }
    if (axdiff < 0 && bxdiff >= 0) {
        return false;
    }
    const auto aydiff = a.y - center.y;
    const auto bydiff = b.y - center.y;
    if (axdiff == 0 && bxdiff == 0) {
        if (aydiff >= 0 || bydiff >= 0) {
            return a.y > b.y;
        }
        return b.y > a.y;
    }
    // compute the cross product of vectors (center -> a) x (center -> b)
    const auto det = axdiff * bydiff - bxdiff * aydiff;
    if (det < 0)
        return true;
    if (det > 0)
        return false;

    // points a and b are on the same line from the center
    // check which point is closer to the center
    const auto d1 = axdiff * axdiff + aydiff * aydiff;
    const auto d2 = bxdiff * bxdiff + bydiff * bydiff;
    return d1 > d2;
}

class TrianglePointOrderComparison {
  public:
    /// Construction from a provided "center" point
    explicit TrianglePointOrderComparison(Point2f const& center) : center_(center) {}
    /// Construction from three points, which are averaged to find the center.
    TrianglePointOrderComparison(Point2f const& a, Point2f const& b, Point2f const& c) : center_((a + b + c) / 3.f) {}

    bool operator()(Point2f const& a, Point2f const& b) const { return orderAfromBwrtCenter(a, b, center_); }

  private:
    Point2f center_;
};

struct ScreenData {
    bool populated = false;
    Point2d screen;
    DataOrigin origin;
};
using Subdiv = SubdivContainer<ScreenData>;
using Vertex = Subdiv::Vertex;
using Vertices = Subdiv::Vertices;

void sort(Point2f const& center, Vertices& vertices) {
    if (vertices.size() < 3) {
        /// only bother sorting containers of 3 or more.
        return;
    }
    TrianglePointOrderComparison comparison(center);
    std::sort(vertices.begin(), vertices.end(),
              [&](Vertex const& a, Vertex const& b) { return comparison(a.location, b.location); });
}

bool hasUnusableData(Vertices const& vertices) {
    for (auto& v : vertices) {
        if (v.status != VertexStatus::AdditionalVertex || !v.hasValue) {
            return true;
        }
        if (!v.value.populated) {
            std::cerr << "Hey, we got an unpopulated value that wasn't caught by earlier checks!" << std::endl;
            return true;
        }
    }
    return false;
}
bool computeWeights(Vertices& vertices, Point2f const& pt) {
    if (hasUnusableData(vertices)) {
        return false;
    }
    if (vertices.size() == 1) {
        vertices.front().weight = 1;
        return true;
    }
    if (vertices.size() == 2) {
        // on the edge.
        auto totalLength = (vertices.back().location - vertices.front().location).norm();
        auto fractionalLength = (pt - vertices.front().location).norm() / totalLength;
        /// fractionalLength is now the weight of vertices.back()
        vertices.front().weight = fractionalLength;
        vertices.back().weight = 1. - fractionalLength;
        return true;
    }
    if (vertices.size() == 3) {
        // in a triangle
        // sort(pt, vertices); /// @todo remove once we have a guarantee of order in the returned vertices.
        auto totalDoubleArea = doubleTriangleArea(vertices[0].location, vertices[1].location, vertices[2].location);
        vertices[0].weight = doubleTriangleArea(vertices[1].location, vertices[2].location, pt) / totalDoubleArea;
        vertices[1].weight = doubleTriangleArea(vertices[0].location, pt, vertices[2].location) / totalDoubleArea;
        vertices[2].weight = doubleTriangleArea(vertices[0].location, vertices[1].location, pt) / totalDoubleArea;
        std::ostringstream os;
        os << "Barycentric coordinates: ";
        static const auto Precision = 4;
        static const auto ColWidth = 10;
        os << std::setprecision(Precision) << std::setw(ColWidth) << vertices[0].weight;
        os << std::setprecision(Precision) << std::setw(ColWidth) << vertices[1].weight;
        os << std::setprecision(Precision) << std::setw(ColWidth) << vertices[2].weight;
        std::cout << __func__ << ": " << os.str() << std::endl;
        return true;
    }
    assert(false && "This shouldn't be reached...");
    return false;
}

Point2d interpolate(Vertices const& vertices) {
    Eigen::Vector2d accum = Eigen::Vector2d::Zero();
    for (auto& v : vertices) {
        accum += (ei::map(v.value.screen) * v.weight);
    }
    Point2d ret;
    ei::map(ret) = accum;
    return ret;
}
int main(int argc, char* argv[]) {
    std::string fn = argv[1];
    auto myFile = std::ifstream(fn);
    auto measurements = readInputMeasurements(fn, myFile);

    Subdiv triangulationData(Rect(-90, -90, 180, 180));
    GenericExtremaFinder<float> longitudeExtrema;
    GenericExtremaFinder<float> latitudeExtrema;
    for (auto& meas : measurements.measurements) {
        ScreenData val;
        val.populated = true;
        val.screen = meas.screen;
        val.origin = meas.getOrigin(measurements);
        auto longitude = static_cast<float>(meas.viewAnglesDegrees.longitude());
        auto latitude = static_cast<float>(meas.viewAnglesDegrees.latitude());
        longitudeExtrema.process(longitude);
        latitudeExtrema.process(latitude);
        triangulationData.insert(Point2f(longitude, latitude), val);
    }

    std::cout << "Longitude range: " << longitudeExtrema << std::endl;
    std::cout << "Latitude range: " << latitudeExtrema << std::endl;

    auto longitudeRange = longitudeExtrema.getMax() - longitudeExtrema.getMin();
    auto latitudeRange = latitudeExtrema.getMax() - latitudeExtrema.getMin();
    auto step = std::min(longitudeRange / STEPS, latitudeRange / STEPS);
    for (std::size_t xStep = 0; xStep * step + longitudeExtrema.getMin() <= longitudeExtrema.getMax(); ++xStep) {
        auto xLong = xStep * step + longitudeExtrema.getMin();
        for (std::size_t yStep = 0; yStep * step + latitudeExtrema.getMin() <= latitudeExtrema.getMax(); ++yStep) {
            auto yLat = yStep * step + latitudeExtrema.getMin();
            const auto pt = Point2f(xLong, yLat);
            auto neighborhood = triangulationData.findNeighborhood(pt);
            auto canInterpolate = computeWeights(neighborhood, pt);
            if (canInterpolate) {
                auto interpolated = interpolate(neighborhood);
            }
        }
    }

    return 0;
}

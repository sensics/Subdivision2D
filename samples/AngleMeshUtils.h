/** @file
    @brief Header

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

#ifndef INCLUDED_AngleMeshUtils_h_GUID_53BAD7DD_E37D_4F89_C221_7A71F60162CC
#define INCLUDED_AngleMeshUtils_h_GUID_53BAD7DD_E37D_4F89_C221_7A71F60162CC

// Internal Includes
// - none

// Library/third-party includes
// - none

// Standard includes
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using Point2d = std::array<double, 2>;

/// Gentle wrapper around Point2d assigning longitude and latitude meaning (respectively) to the elements.
struct LongLat {
    Point2d longLat;
    /// angle in x
    double& longitude() { return longLat[0]; }
    /// angle in x (read-only)
    double longitude() const { return longLat[0]; }

    /// angle in y
    double& latitude() { return longLat[1]; }
    /// angle in y (read-only)
    double latitude() const { return longLat[1]; }
};

/// Convenient storage for the input source (typically file name) and line number associated with some measurement or
/// (more typically) a derived quantity.
struct DataOrigin {
    std::string inputSource;
    std::size_t lineNumber;

    /// Is the origin of this data known? (or default-constructed aka unknown?)
    bool known() const { return !inputSource.empty(); }
};

/// Stream insertion operator for DataOrigin
inline std::ostream& operator<<(std::ostream& os, DataOrigin const& orig) {
    if (orig.known()) {
        std::ostringstream oss;
        oss << orig.inputSource << ":" << orig.lineNumber;
        os << oss.str();
    } else {
        os << "(unknown)";
    }
    return os;
}

struct InputMeasurements;
struct InputMeasurement {
    /// In arbitrary units
    Point2d screen;

    /// in degrees (either field angles or longitude/latitude, depending on config option)
    LongLat viewAnglesDegrees;

    /// Line number in loaded file
    std::size_t lineNumber = 0;

    /// Provide the parent container, get a single object usable to refer to the source of this measurement.
    DataOrigin getOrigin(InputMeasurements const& parent) const;
};

struct InputMeasurements {
    /// Filename that these measurements were loaded from.
    std::string inputSource;
    /// Collection of measurements as they were loaded.
    std::vector<InputMeasurement> measurements;

    /// Is the measurement collection empty?
    bool empty() const { return measurements.empty(); }
    /// How many measurements do we have?
    std::size_t size() const { return measurements.size(); }
};

inline DataOrigin InputMeasurement::getOrigin(InputMeasurements const& parent) const {
    // Out of line to allow declaration of InputMeasurements.
    return DataOrigin{parent.inputSource, lineNumber};
}

/// Reads the four-column whitespace-delimited mapping file.
/// Returns empty vector if it fails to read anything.
InputMeasurements readInputMeasurements(std::string const& inputSource, std::istream& in);

#endif // INCLUDED_AngleMeshUtils_h_GUID_53BAD7DD_E37D_4F89_C221_7A71F60162CC

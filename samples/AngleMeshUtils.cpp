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

// Library/third-party includes
// - none

// Standard includes
#include <cassert>
#include <iostream>

namespace {

enum class LineParseResult { Success, LongitudeError, LatitudeError, ScreenXError, ScreenYError };
static inline const char* to_string(LineParseResult result) {
    switch (result) {
    case LineParseResult::Success:
        return "success (no error)";
    case LineParseResult::LongitudeError:
        return "longitude (x angle, first column) error";
    case LineParseResult::LatitudeError:
        return "latitude (y angle, second column) error";
    case LineParseResult::ScreenXError:
        return "screen x position (third column) error";
    case LineParseResult::ScreenYError:
        return "screen y position (fourth column) error";
    default:
        assert(0 && "Should never happen!");
        return "Unknown result! Should be impossible!";
    }
}

static inline LineParseResult parseInputMeasurementLine(std::string const& line, InputMeasurement& meas) {
    std::istringstream lineStream(line);
    // Read the mapping info from the input file.
    if (!(lineStream >> meas.viewAnglesDegrees.longitude())) {
        return LineParseResult::LongitudeError;
    }
    if (!(lineStream >> meas.viewAnglesDegrees.latitude())) {
        return LineParseResult::LatitudeError;
    }
    if (!(lineStream >> meas.screen[0])) {
        return LineParseResult::ScreenXError;
    }
    if (!(lineStream >> meas.screen[1])) {
        return LineParseResult::ScreenYError;
    }
    return LineParseResult::Success;
}
} // namespace

InputMeasurements readInputMeasurements(std::string const& inputSource, std::istream& in) {
    InputMeasurements ret;
    ret.inputSource = inputSource;
    std::string line;
    std::getline(in, line);
    // Read first line as initial condition,
    // termination condition is reading a line failed the stream,
    // and iteration is reading another line and incrementing line number.
    // We read each line individually so that additional fields can be added at the end of the line and not mess up this
    // parsing.
    for (std::size_t lineNumber = 1; !in.eof() && in.good(); std::getline(in, line), ++lineNumber) {
        // Read the mapping info from the input file.
        InputMeasurement meas;
        meas.lineNumber = lineNumber;
        auto parseResult = parseInputMeasurementLine(line, meas);
        if (LineParseResult::Success == parseResult) {
            // if we didn't fail sometime here...
            ret.measurements.push_back(meas);
        } else {
            std::cerr << "Failed to parse line  " << lineNumber << " with a \"" << to_string(parseResult) << "\": '"
                      << line << "'" << std::endl;
            std::cerr << "longitude: " << meas.viewAnglesDegrees.longitude() << "\n"
                      << "latitude: " << meas.viewAnglesDegrees.latitude() << "\n"
                      << "screen: " << meas.screen[0] << "\t " << meas.screen[1] << std::endl;
        }
    }
    std::cerr << "Read " << ret.size() << " lines..." << std::endl;

    return ret;
}

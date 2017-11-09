/** @file
    @brief Header containing suppor types. Point2 is written from scratch. Rect_ is based on the OpenCV equivalent.

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

#ifndef INCLUDED_Types_h_GUID_EE6A183B_3E3E_4B8B_28FB_A6114F219F51
#define INCLUDED_Types_h_GUID_EE6A183B_3E3E_4B8B_28FB_A6114F219F51

// Internal Includes
// - none

// Library/third-party includes
// - none

// Standard includes
#include <cfloat>
#include <limits>

namespace sensics {
namespace subdiv2d {
    template <typename T> struct Point2 {
        using value_type = T;
        Point2(value_type xVal, value_type yVal) : x(xVal), y(yVal) {}
        Point2() : x(std::numeric_limits<T>::(max)()), y(std::numeric_limits<T>::(max)()) {}

        value_type x;
        value_type y;

        Point2& operator-=(Point2 const& rhs) {
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }
    };
    template <typename T> inline Point2<T> operator-(Point2<T> const& lhs, Point2<T> const& rhs) {
        auto ret = Point2<T>(lhs);
        ret -= rhs;
        return ret;
    }

    using Point2f = Point2<float>;

} // namespace subdiv2d
} // namespace sensics
#endif // INCLUDED_Types_h_GUID_EE6A183B_3E3E_4B8B_28FB_A6114F219F51

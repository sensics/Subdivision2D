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

#ifndef INCLUDED_GenericExtremaFinder_h_GUID_F8106110_8A9A_4AA0_0313_A42B04FFD1D0
#define INCLUDED_GenericExtremaFinder_h_GUID_F8106110_8A9A_4AA0_0313_A42B04FFD1D0

// Internal Includes
// - none

// Library/third-party includes
// - none

// Standard includes
#include <cassert>
#include <utility>

template <typename ValueType, typename Comparator = std::less<ValueType>> class GenericExtremaFinder {
  public:
    using comparator_type = Comparator;
    using value_type = ValueType;
    GenericExtremaFinder() : compare_() {}
    GenericExtremaFinder(comparator_type&& compare) : compare_(std::move(compare)) {}
    GenericExtremaFinder(comparator_type const& compare) : compare_(compare) {}
    bool valid() const { return valid_; }
    value_type const& getMin() const {
        assert(valid_ && "Only makes sense to get min value if you've actually processed any elements!");
        return minVal_;
    }
    value_type const& getMax() const {
        assert(valid_ && "Only makes sense to get max value if you've actually processed any elements!");
        return maxVal_;
    }
    void process(value_type const& val) {
        if (!valid_) {
            minVal_ = val;
            maxVal_ = val;
            valid_ = true;
            return;
        }

        if (compare_(val, minVal_)) {
            /// new is less than min
            minVal_ = val;
        }

        if (compare_(maxVal_, val)) {
            /// max is less than new
            maxVal_ = val;
        }
    }

  private:
    comparator_type compare_;
    bool valid_ = false;

    value_type minVal_;
    value_type maxVal_;
};

template <typename Stream, typename ValueType, typename Comparator>
Stream& operator<<(Stream& os, GenericExtremaFinder<ValueType, Comparator> const& extrema) {
    os << "[";
    if (extrema.valid()) {
        os << extrema.getMin() << ", " << extrema.getMax();
    } else {
        os << "NO DATA";
    }
    os << "]";
    return os;
}

#endif // INCLUDED_GenericExtremaFinder_h_GUID_F8106110_8A9A_4AA0_0313_A42B04FFD1D0

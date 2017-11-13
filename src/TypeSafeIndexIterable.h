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

#ifndef INCLUDED_TypeSafeIndexIterable_h_GUID_3184D13C_A6AD_4889_FED8_D6D44FF89FEE
#define INCLUDED_TypeSafeIndexIterable_h_GUID_3184D13C_A6AD_4889_FED8_D6D44FF89FEE

// Internal Includes
#include "TypeSafeIndex.h"

// Library/third-party includes
// - none

// Standard includes
#include <iterator>

namespace sensics {
namespace detail {
    template <typename Tag> class TypeSafeIndexRange {
      public:
        using value_type = TypeSafeIndex<Tag>;
        using value_wrapped_type = typename value_type::value_type;
        /// Input iterator because dereference returns value_type itself, rather than a reference type.
        using iterator_tag = std::input_iterator_tag;
        class iterator : public std::iterator<iterator_tag, value_type, value_wrapped_type, const value_wrapped_type*,
                                              value_wrapped_type> {
          public:
            /// produces a universal end iterator
            iterator() = default;
            iterator(value_type val, value_type pastTheEnd) : val_(val), end_(pastTheEnd) {
                if (!val) {
                    reset();
                    return;
                }
                if (!pastTheEnd) {
                    reset();
                    return;
                }
            }

            value_type operator*() const {
                assert(!isEndIterator() && "Can't dereference end iterator!");
                return val_;
            }

            /// utility function, not required for iterators in general, just this one to make it simpler.
            bool isEndIterator() const { return (val_ == end_) || (!val_.valid() && !end_.valid()); }

            /// Pre-increment.
            iterator& operator++() {
                assert(!isEndIterator() && "Can't increment end iterator!");
                val_ = value_type(val_.get() + 1);
                return *this;
            }
            /// Post-increment, icky.
            iterator operator++(int) {
                iterator ret(*this);
                assert(!isEndIterator() && "Can't increment end iterator!");
                val_ = value_type(val_.get() + 1);
                return ret;
            }

            /// test for equality
            bool operator==(iterator const& rhs) const {
                return (isEndIterator() && rhs.isEndIterator()) || (val_ == rhs.val_);
            }

            /// test for inequality
            bool operator!=(iterator const& rhs) const {
                if (isEndIterator() && rhs.isEndIterator()) {
                    /// only case where val_ might differ but they're actually equal.
                    return false;
                }
                return (val_ != rhs.val_);
            }

          private:
            void reset() {
                val_ = value_type();
                end_ = value_type();
            }
            /// @name State
            /// @brief current value and "beyond the end" value.
            /// Either both are full or both are empty.
            /// If both are empty, or they are equal
            /// @{
            value_type val_;
            value_type end_;
            /// @}
        };

        TypeSafeIndexRange(value_type first, value_type pastEnd) : begin_(first), end_(pastEnd) {
            assert(!first.empty() && "Can't iterate from an empty start");
            assert(!pastEnd.empty() && "Can't iterate to an empty end");
            assert((first < pastEnd || first == pastEnd) && "First must be less than end!");
        }
        iterator begin() const { return iterator(begin_, end_); }

        iterator end() const {
            // could also do iterator(end_, end_)
            return iterator();
        }

      private:
        value_type begin_;
        value_type end_;
    };
    /// Returns an iterable range containing integer indices [first.value(), last.value()]   (that is, all values
    /// inclusive of the endpoints)
    template <typename Tag>
    inline TypeSafeIndexRange<Tag> inclusiveIndexRange(TypeSafeIndex<Tag> first, TypeSafeIndexRange<Tag> last) {
        auto pastTheEnd = TypeSafeIndex<Tag>(last.get() + 1);
        return TypeSafeIdRange<Tag>(first, pastTheEnd);
    }

    /// Returns an iterable range containing integer indices [first.value(), pastTheEnd.value())   (that is, a half-open
    /// range)
    template <typename Tag>
    inline TypeSafeIndexRange<Tag> indexRange(TypeSafeIndex<Tag> first, TypeSafeIndex<Tag> pastTheEnd) {
        return TypeSafeIndexRange<Tag>(first, pastTheEnd);
    }

    template <typename TypeSafeIdType>
    using RangeOfTypeSafeIndex = TypeSafeIndexRange<typename TypeSafeIdType::tag_type>;
} // namespace detail
using detail::inclusiveIndexRange;
using detail::indexRange;

} // namespace sensics
#endif // INCLUDED_TypeSafeIndexIterable_h_GUID_3184D13C_A6AD_4889_FED8_D6D44FF89FEE

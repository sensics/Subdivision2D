/** @file
    @brief Header

    Related to osvr::util::TypeSafeId in spirit, but totally separate implementation.


    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
// SPDX-License-Identifier:BSD-3-Clause

#ifndef INCLUDED_TypeSafeIndex_h_GUID_DC93F7B0_1715_41A3_9B68_93531C657526
#define INCLUDED_TypeSafeIndex_h_GUID_DC93F7B0_1715_41A3_9B68_93531C657526

// Internal Includes
// - none

// Library/third-party includes
// - none

// Standard includes
#include <cstddef>

namespace sensics {
namespace detail {
    /// specialize to change the name used in IO operations to indicate this ID.
    /// Default, not defined in this file, uses typeid and name() on your tag.
    template <typename Tag> struct TypeSafeIndexNameTrait;

    /// function for easier use of TypeSafeIndexNameTrait
    template <typename Tag> constexpr static inline const char* getTypeSafeIndexName() {
        return TypeSafeIndexNameTrait<Tag>::get();
    }

    /// specialize to change the contained type for a class using the given tag.
    template <typename Tag> struct TypeSafeIndexValueTypeTrait { using type = std::size_t; };

    /// alias for easier use of TypeSafeIndexValueTypeTrait
    template <typename Tag> using TypeSafeIndexValueType = typename TypeSafeIndexValueTypeTrait<Tag>::type;

    /// specialize to change which (ideally "invalid") value to use as the initial value for a class using the given
    /// tag.
    template <typename Tag> struct TypeSafeIndexInitValueTrait {
        static constexpr TypeSafeIndexValueType<Tag> get() { return static_cast<TypeSafeIndexValueType<Tag>>(0); }
    };

    /// function for easier use of TypeSafeIndexInitValueTrait
    template <typename Tag> constexpr static inline TypeSafeIndexValueType<Tag> getTypeSafeIndexInitValue() {
        return TypeSafeIndexInitValueTrait<Tag>::get();
    }

    /// specialize to change the predicate to determine if a value is "valid" for a class using the given tag.
    /// default is just equality comparison with the invalid value, which allows TypeSafeIndexFastEquality.
    template <typename Tag> struct TypeSafeIndexIsValidTrait {
        /// don't put this in your specialization - since you wouldn't be specializing if you were comparing equality
        /// with just one value...
        using can_use_fast_equality = void;
        static constexpr bool get(TypeSafeIndexValueType<Tag> val) { return val != getTypeSafeIndexInitValue<Tag>(); }
    };

    /// Generalized equality: first checks validity (invalid == invalid, invalid != valid, always) before comparing
    /// values.
    template <typename Tag> struct TypeSafeIndexGeneralEquality;
    /// Marginally faster alternative, but only valid where there's exactly 1 invalid value. (Otherwise two invalid
    /// values might compare to be unequal.)
    template <typename Tag> struct TypeSafeIndexFastEquality;

    /// Select general equality by default.
    template <typename Tag, typename Dummy = void> struct TypeSafeIndexEquality : TypeSafeIndexGeneralEquality<Tag> {};

    /// Only when there's a "using can_use_fast_equality = void", enable fast equality.
    template <typename Tag>
    struct TypeSafeIndexEquality<Tag, typename TypeSafeIndexIsValidTrait<Tag>::can_use_fast_equality>
        : TypeSafeIndexFastEquality<Tag> {};

    /// Main template class - pass a tag type (can be incomplete) to create mutually-incompatible types.
    template <typename Tag> class TypeSafeIndex {
      public:
        using value_type = TypeSafeIndexValueType<Tag>;
        using tag_type = Tag;
        using type = TypeSafeIndex<Tag>;
        /// default construct
        TypeSafeIndex() = default;
        /// explicit construction with a value of the right type.
        explicit TypeSafeIndex(value_type val) : val_(val) {}
        /// copy construct
        TypeSafeIndex(type const& other) : val_(other.val_) {}
        /// assign
        TypeSafeIndex& operator=(type const& other) {
            val_ = other.val_;
            return *this;
        }
        /// swap
        void swap(type& other) { std::swap(val_, other.val_); }

        /// Is the contained value valid?
        bool valid() const { return TypeSafeIndexIsValidTrait<Tag>::get(val_); }
        /// Is the contained value valid? Explicit bool conversion for testing in conditionals, etc.
        explicit operator bool() const { return valid(); }

        /// Get contained value, whether or not it is valid
        value_type get() const { return val_; }

        /// Get contained value, whether or not it is valid
        value_type value() const { return val_; }

      private:
        value_type val_ = getTypeSafeIndexInitValue<Tag>();
    };

    /// free function for swap.
    template <typename Tag> static inline void swap(TypeSafeIndex<Tag>& lhs, TypeSafeIndex<Tag>& rhs) {
        lhs.swap(rhs) :
    }

    // Implementation of general equality.
    template <typename Tag> struct TypeSafeIndexGeneralEquality {
        using type = TypeSafeIndex<Tag>;
        static bool get(type const& lhs, type const& rhs) {
            if (!lhs.valid() && !rhs.valid()) {
                // Both are invalid: then they are equal.
                return true;
            }
            if (!lhs.valid() || !rhs.valid()) {
                // Both invalid already pruned out - so here, only one would be invalid,
                // so unequal.
                return false;
            }
            // Now compare the values: they're both valid by now.
            return lhs.get() == rhs.get();
        }
    };
    // Implementation of fast equality.
    template <typename Tag> struct TypeSafeIndexFastEquality {
        using type = TypeSafeIndex<Tag>;
        static bool get(type const& lhs, type const& rhs) {
            /// if there's only one invalid value, this handles all the cases in one comparison.
            return lhs.get() == rhs.get();
        }
    };
    template <typename Tag>
    static inline bool operator==(TypeSafeIndex<Tag> const& lhs, TypeSafeIndex<Tag> const& rhs) {
        return TypeSafeIndexEquality<Tag>::get(lhs, rhs);
    }
    template <typename Tag>
    static inline bool operator!=(TypeSafeIndex<Tag> const& lhs, TypeSafeIndex<Tag> const& rhs) {
        return !TypeSafeIndexEquality<Tag>::get(lhs, rhs);
    }

} // namespace detail

} // namespace sensics

#endif // INCLUDED_TypeSafeIndex_h_GUID_DC93F7B0_1715_41A3_9B68_93531C657526

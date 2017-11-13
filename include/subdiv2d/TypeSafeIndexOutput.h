/** @file
    @brief Header

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
// SPDX-License-Identifier:BSD-3-Clause

#ifndef INCLUDED_TypeSafeIndexOutput_h_GUID_89C636F4_7DFC_418C_11E5_703E440926EB
#define INCLUDED_TypeSafeIndexOutput_h_GUID_89C636F4_7DFC_418C_11E5_703E440926EB

// Internal Includes
#include "TypeSafeIndex.h"

// Library/third-party includes
// - none

// Standard includes
#include <sstream>
#include <typeinfo>

namespace sensics {
namespace detail {
    /// specialize to change the name used in IO operations to indicate this ID.
    template <typename Tag> struct TypeSafeIndexNameTrait {
        static constexpr const char* get() { return typeid(Tag).name(); }
    };

    template <typename Tag> static inline std::ostream& operator<<(std::ostream& os, TypeSafeIndex<Tag> const& val) {
        std::ostringstream oss;
        oss << getTypeSafeIndexName<Tag>() << "(";
        if (val.valid()) {
            oss << val.value();
        } else {
            oss << "NULL";
        }
        oss << ")";
        os << oss.str();
        return os;
    }

} // namespace detail

} // namespace sensics
#endif // INCLUDED_TypeSafeIndexOutput_h_GUID_89C636F4_7DFC_418C_11E5_703E440926EB

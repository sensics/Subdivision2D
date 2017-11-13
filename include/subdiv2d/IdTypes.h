/** @file
    @brief Header for the various ID types that are not interchangeable but derive from indices.

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
// SPDX-License-Identifier:BSD-3-Clause

#ifndef INCLUDED_IdTypes_h_GUID_331FBAA0_4F87_4690_ED9C_83D68CF54785
#define INCLUDED_IdTypes_h_GUID_331FBAA0_4F87_4690_ED9C_83D68CF54785

// Internal Includes
#include "TypeSafeIndex.h"
#include "TypeSafeIndexOutput.h"

// Library/third-party includes
// - none

// Standard includes
// - none

namespace sensics {

namespace subdiv2d {
    namespace detail {
        struct VertexTag;
        struct EdgeTag;
        struct QuadEdgeTag;
    } // namespace detail
} // namespace subdiv2d

namespace detail {
    template <> struct TypeSafeIndexNameTrait<subdiv2d::detail::VertexTag> {
        static constexpr const char* get() { return "VertexId"; }
    };

    template <> struct TypeSafeIndexNameTrait<subdiv2d::detail::EdgeTag> {
        static constexpr const char* get() { return "EdgeId"; }
    };

    template <> struct TypeSafeIndexNameTrait<subdiv2d::detail::QuadEdgeTag> {
        static constexpr const char* get() { return "QuadEdgeId"; }
    };

    /// Right now, using int for Vertex ids, and <= 0 as invalid.
    template <> struct TypeSafeIndexValueTypeTrait<subdiv2d::detail::VertexTag> { using type = int; };
    template <> struct TypeSafeIndexIsValidTrait<subdiv2d::detail::VertexTag> {
        static constexpr bool get(TypeSafeIndexValueType<subdiv2d::detail::VertexTag> val) { return val <= 0; }
    };

    /// Right now, using int for Edge ids, and <= 0 as invalid.
    template <> struct TypeSafeIndexValueTypeTrait<subdiv2d::detail::EdgeTag> { using type = int; };
    template <> struct TypeSafeIndexIsValidTrait<subdiv2d::detail::EdgeTag> {
        static constexpr bool get(TypeSafeIndexValueType<subdiv2d::detail::EdgeTag> val) { return val <= 0; }
    };
    /// Right now, using int for QuadEdge ids, and <= 0 as invalid.
    template <> struct TypeSafeIndexValueTypeTrait<subdiv2d::detail::QuadEdgeTag> { using type = int; };

} // namespace detail

namespace subdiv2d {
    using VertexId = ::sensics::detail::TypeSafeIndex<detail::VertexTag>;
    using EdgeId = ::sensics::detail::TypeSafeIndex<detail::EdgeTag>;
    using QuadEdgeId = ::sensics::detail::TypeSafeIndex<detail::QuadEdgeTag>;

    static const VertexId InvalidVertex = VertexId();
    static const EdgeId InvalidEdge = EdgeId();
} // namespace subdiv2d

} // namespace sensics
#endif // INCLUDED_IdTypes_h_GUID_331FBAA0_4F87_4690_ED9C_83D68CF54785

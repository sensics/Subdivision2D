/** @file
    @brief Implementation

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
// SPDX-License-Identifier:BSD-3-Clause

// Internal Includes
#include <subdiv2d/SubdivContainer.h>

// Library/third-party includes
// - none

// Standard includes
// - none

namespace sensics {
namespace subdiv2d {
    SubdivContainerBase::SubdivContainerBase(Rect bounds) {}

    VertexStatus SubdivContainerBase::categorizeVertex(VertexId id) {
        if (!id.valid()) {
            return VertexStatus::Unpopulated;
        }
        if (id.get() < NumDummyVertices) {
            return VertexStatus::OuterBoundingVertex;
        }
        return VertexStatus::AdditionalVertex;
    }

} // namespace subdiv2d
} // namespace sensics

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
    SubdivContainerBase::SubdivContainerBase(Rect bounds) : subdiv_(bounds) {}

    VertexStatus SubdivContainerBase::categorizeVertex(VertexId id) {
        if (!id.valid()) {
            return VertexStatus::Unpopulated;
        }
        if (id.get() < NumDummyVertices) {
            return VertexStatus::OuterBoundingVertex;
        }
        return VertexStatus::AdditionalVertex;
    }

    VertexValueId SubdivContainerBase::getValueId(VertexId id) {
        if (!id.valid() || id.get() < NumDummyVertices) {
            return InvalidVertexValueId;
        }
        return VertexValueId(static_cast<size_t>(id.get() - NumDummyVertices));
    }

    bool SubdivContainerBase::hasValue(VertexValueId valueId) const {
        if (!valueId) {
            return false;
        }
        return (valueId.get() < valuesSet_.size()) && (valuesSet_[valueId.get()]);
    }

    std::size_t SubdivContainerBase::set(VertexValueId valueId) {
        std::size_t ret = NoResizingNeededSentinel;
        if (valueId.get() >= valuesSet_.size()) {
            ret = valueId.get() + 1;
            valuesSet_.resize(ret, false);
        }
        valuesSet_[valueId.get()] = true;
        return ret;
    }

    SubdivContainerBase::BaseVector SubdivContainerBase::locateNeighborhood(Point2f const& pt) {
        BaseVector ret;
#if 0
        auto vertexIds = subdiv_.locateVertexIdsArray(pt);
#else
        auto vertexIds = subdiv_.locateVertexIdsForInterpolationArray(pt);
#endif
        std::cout << "Vertices ";
        for (auto v : vertexIds) {
            std::cout << " " << v;
            if (!v) {
                /// invalid vertex id
                continue;
            }
            ret.push_back(ContainerVertexBase());
            populate(v, ret.back());
        }
        std::cout << std::endl;
        return ret;
    }

    void SubdivContainerBase::populate(VertexId id, ContainerVertexBase& data) {
        data.status = categorizeVertex(id);
        data.id = id;
        if (data.status != VertexStatus::Unpopulated) {
            data.location = subdiv_.getVertex(id);
        }
        if (data.status >= VertexStatus::AdditionalVertex) {
            data.valueId = getValueId(id);
            data.hasValue = hasValue(data.valueId);
        }
    }

    VertexValueId SubdivContainerBase::lookup(Point2f const& pt) {
        if (subdiv_.empty()) {
            return InvalidVertexValueId;
        }
        auto edge = InvalidEdge;
        auto vertex = InvalidVertex;
        const auto status = subdiv_.locate(pt, edge, vertex);
        if (status == PtLoc::PTLOC_VERTEX) {
            return getValueId(vertex);
        }
        return InvalidVertexValueId;
    }

    std::pair<VertexValueId, std::size_t> SubdivContainerBase::insert(Point2f const& pt) {
        std::size_t outSize = NoResizingNeededSentinel;
        auto ptId = subdiv_.insert(pt);
        if (!ptId) {
            return std::make_pair(InvalidVertexValueId, outSize);
        }

        auto valueId = getValueId(ptId);
        if (!valueId) {
            return std::make_pair(valueId, outSize);
        }
        outSize = set(valueId);
        return std::make_pair(valueId, outSize);
    }

} // namespace subdiv2d
} // namespace sensics

/** @file
    @brief Header defining a template class container associating values with points on a 2D subdivision (triangulation)

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr

    SPDX-License-Identifier:BSD-3-Clause
*/

// Copyright 2017 Sensics, Inc.

#ifndef INCLUDED_SubdivContainer_h_GUID_B7F5EEE6_E316_4EA6_5A20_EDADEE3D99A9
#define INCLUDED_SubdivContainer_h_GUID_B7F5EEE6_E316_4EA6_5A20_EDADEE3D99A9

// Internal Includes
#include "FixedMaxSizeArray.h"
#include "IdTypes.h"
#include "Subdivision2D.h"

// Library/third-party includes
// - none

// Standard includes
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sensics {
namespace subdiv2d {
    using VertexValueId = ::sensics::detail::TypeSafeIndex<detail::VertexValueTag>;

    enum class VertexStatus {
        /// Invalid/null vertex
        /// Neither value nor location.
        /// Sometimes used for padding so we can use a fixed-size container instead of a std::vector when we know a max
        /// return size.
        Unpopulated,
        /// One of the initial dummy vertices added implicitly to form an outer triangle.
        /// These vertices have locations (out of bounds, but known), but no values associated.
        OuterBoundingVertex,
        /// A non-dummy vertex - may still be a "virtual" vertex however (?)
        /// These are the only vertices that may have associated values.
        AdditionalVertex
    };
    struct ContainerVertexBase {
        VertexStatus status = VertexStatus::Unpopulated;
        VertexId id = InvalidVertex;
        VertexValueId valueId = InvalidVertexValueId;
        bool hasValue = false;
        Point2f location;
    };
    template <typename T> struct ContainerVertexValue : ContainerVertexBase {
        ContainerVertexValue() = default;
        ContainerVertexValue(ContainerVertexBase const& base) : ContainerVertexBase(base) {}

        using value_type = T;
        value_type value;
    };

    static const std::size_t MaxNeighborhoodSize = 3;
    template <typename T> using ContainerVertices = MaxSizeVector<ContainerVertexValue<T>, MaxNeighborhoodSize>;

    class SubdivContainerBase {
      public:
        explicit SubdivContainerBase(Rect bounds);

      protected:
        VertexStatus categorizeVertex(VertexId id);
        VertexValueId getValueId(VertexId id);
        bool hasValue(VertexValueId valueId) const;

        using BaseVector = MaxSizeVector<ContainerVertexBase, MaxNeighborhoodSize>;
        BaseVector locateNeighborhood(Point2f const& pt);

        // Lookup a vertex by location. If it exists, a valid VertexValueId will be returned. (It may be that no value
        // is set for that id - that's a separate question/call)
        VertexValueId lookup(Point2f const& pt);

        static const std::size_t NoResizingNeededSentinel = 0;

        // Insert a new point into the subdivision. If no errors occur, a valid VertexValueId will be returned. outSize
        // contains the size to resize your vector to, or 0 (NoResizingNeededSentinel) if no resize is needed.
        std::pair<VertexValueId, std::size_t> insert(Point2f const& pt);

      private:
        void populate(VertexId id, ContainerVertexBase& data);
        /// returns the size the value vector should be resized to, or 0 if no resizing needed.
        std::size_t set(VertexValueId valueId);
        static const int NumDummyVertices = 4;
        Subdiv2D subdiv_;
        std::vector<bool> valuesSet_;
    };

    template <typename T> class SubdivContainer : public SubdivContainerBase {
      public:
        using value_type = T;
        using pointer_type = T*;
        using Vertex = ContainerVertexValue<value_type>;
        using Vertices = ContainerVertices<value_type>;

        explicit SubdivContainer(Rect bounds);

        /// Insert a new point into the subdivision, along with its associated value. If it is outside the bounds, a
        /// runtime error is raised. If it is an already-existing point, the value will be replaced.
        void insert(Point2f const& pt, value_type const& val);

        /// Returns true if the given point is a vertex in the subdivision and outVal has been set. Returns false in all
        /// other cases.
        bool lookup(Point2f const& pt, value_type& outVal);

        /// Returns true if the given point is a vertex in the subdivision. Iff outPtr is not nullptr and it will return
        /// true, outPtr is assigned the associated value.  Returns false in all other cases (point not found/out of
        /// bounds/etc).
        bool lookup(Point2f const& pt, pointer_type outPtr = nullptr);

        /// Gets the value associated with a point. If no value is associated with the given point, throws a runtime
        /// error.
        value_type get(Point2f const& pt);

        /// If the point is a vertex in the subdivision, return just the point and its value in outVertices.
        /// If the point is on an edge, return the vertices and values at either end of the edge.
        /// If the point is in some facet, return the three vertices and values of that facet.
        /// Otherwise (out of bounds, etc), return an empty container.
        ContainerVertices<value_type> findNeighborhood(Point2f const& pt);
#if 0
        /// If the point is a vertex in the subdivision, return just the point and its value in outVertices.
        /// If the point is on an edge, return the vertices and values at either end of the edge.
        /// If the point is in some facet, return the three vertices and values of that facet.
        /// Otherwise (out of bounds, etc), return false without touching outVertices.
        bool findNearest(Point2f const& pt, Vertices& outVertices);
#endif
      private:
        using Base = SubdivContainerBase;
        bool get_(VertexValueId valueId, value_type& outVal);
        bool get_(VertexValueId valueId, pointer_type outPtr = nullptr);
        bool get_(VertexId vertexId, value_type& outVal);
        bool get_(VertexId vertexId, pointer_type outPtr = nullptr);
#if 0
        bool setFromVertex_(Point2f const& pt, VertexId vertexId, Vertices& outVertices);
        bool setFromEdge_(EdgeId edgeId, Vertices& outVertices);
        bool setFromFacet_(Point2f const& pt, EdgeId edgeId, Vertices& outVertices);
#endif
        std::vector<value_type> associatedValues_;
    };

    template <typename T> inline SubdivContainer<T>::SubdivContainer(Rect bounds) : SubdivContainerBase(bounds) {}
    template <typename T> inline void SubdivContainer<T>::insert(Point2f const& pt, value_type const& val) {
        std::size_t newSize;
        VertexValueId valueId;
        std::tie(valueId, newSize) = Base::insert(pt);
        if (!valueId) {
            return;
        }
        if (newSize != NoResizingNeededSentinel) {
            associatedValues_.resize(newSize);
        }
        associatedValues_[valueId.get()] = val;
    }
    template <typename T> inline bool SubdivContainer<T>::lookup(Point2f const& pt, value_type& outVal) {
        return lookup(pt, &outVal);
    }
    template <typename T> inline bool SubdivContainer<T>::lookup(Point2f const& pt, pointer_type outPtr) {
        auto valueId = Base::lookup(pt);
        if (valueId && hasValue(valueId)) {
            return get_(valueId, outPtr);
        }
        return false;
    }
    template <typename T> inline typename SubdivContainer<T>::value_type SubdivContainer<T>::get(Point2f const& pt) {
        value_type ret;
        if (!lookup(pt, ret)) {
            throw std::runtime_error("Vertex not found with a value in subdivision");
        }
        return ret;
    }
    template <typename T> inline ContainerVertices<T> SubdivContainer<T>::findNeighborhood(Point2f const& pt) {
        auto baseNeighborhood = Base::locateNeighborhood(pt);
        Vertices ret;
        for (auto& baseData : baseNeighborhood) {
            ret.push_back(baseData);
            if (baseData.hasValue) {
                get_(baseData.valueId, ret.back().value);
            }
        }
        return ret;
    }

    template <typename T> inline bool SubdivContainer<T>::get_(VertexValueId valueId, value_type& outVal) {
        return get_(valueId, &outVal);
    }
    template <typename T> inline bool SubdivContainer<T>::get_(VertexValueId valueId, pointer_type outPtr) {
        if (!hasValue(valueId)) {
            return false;
        }
        if (outPtr) {
            *outPtr = associatedValues_[valueId.get()];
        }
        return true;
    }
    template <typename T> inline bool SubdivContainer<T>::get_(VertexId vertexId, value_type& outVal) {
        return get_(getValueId(vertexId), &outVal);
    }
    template <typename T> inline bool SubdivContainer<T>::get_(VertexId vertexId, pointer_type outPtr) {
        return get_(getValueId(vertexId), outPtr);
    }
#if 0
    template <typename T> inline bool SubdivContainer<T>::findNearest(Point2f const& pt, Vertices& outVertices) {
        auto vertices = subdiv_.locateVertexIdsArray(pt);
        std::cout << "Vertices ";
        for (auto& v : vertices) {
            std::cout << " " << v;
        }
        std::cout << std::endl;
        auto locateResult = subdiv_.locate(pt);
        switch (std::get<PtLoc>(locateResult)) {
        case PtLoc::PTLOC_INSIDE:
            return setFromFacet_(pt, std::get<EdgeId>(locateResult), outVertices);
        case PtLoc::PTLOC_ON_EDGE:
            return setFromEdge_(std::get<EdgeId>(locateResult), outVertices);
        case PtLoc::PTLOC_VERTEX:
            return setFromVertex_(pt, std::get<VertexId>(locateResult), outVertices);
        default:
            return false;
        }
    }
    template <typename T>
    inline bool SubdivContainer<T>::setFromVertex_(Point2f const& pt, VertexId vertexId, Vertices& outVertices) {

        Vertex v;
        v.location = pt;
        if (!get_(vertexId, v.value)) {
            return false;
        }
        outVertices = {std::move(v)};
        return true;
    }
    template <typename T> inline bool SubdivContainer<T>::setFromEdge_(EdgeId edgeId, Vertices& outVertices) {
        Vertex org;
        int orgId = subdiv_.edgeOrg(edgeId, &org.location);
        if (!get_(orgId, org.value)) {
            return false;
        }
        Vertex dst;
        int dstId = subdiv_.edgeDst(edgeId, &dst.location);
        if (!get_(dstId, dst.value)) {
            return false;
        }
        outVertices = {std::move(org), std::move(dst)};
        return true;
    }
    template <typename T>
    inline bool SubdivContainer<T>::setFromFacet_(Point2f const& pt, EdgeId edgeId, Vertices& outVertices) {
        /// @todo
        return false;
    }
#endif
} // namespace subdiv2d
} // namespace sensics
#endif // INCLUDED_SubdivContainer_h_GUID_B7F5EEE6_E316_4EA6_5A20_EDADEE3D99A9

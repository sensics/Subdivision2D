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
#include "Subdivision2D.h"

// Library/third-party includes
// - none

// Standard includes
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sensics {
namespace subdiv2d {
    template <typename T> class SubdivContainer {
      public:
        using value_type = T;
        using pointer_type = T*;
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

        struct Vertex {
            Point2f location;
            value_type value;
        };
        using Vertices = std::vector<Vertex>;

        /// If the point is a vertex in the subdivision, return just the point and its value in outVertices.
        /// If the point is on an edge, return the vertices and values at either end of the edge.
        /// If the point is in some facet, return the three vertices and values of that facet.
        /// Otherwise (out of bounds, etc), return false without touching outVertices.
        bool findNearest(Point2f const& pt, Vertices& outVertices);

      private:
        static const int NumDummyVertices = 4;
        bool get_(int vertexId, value_type& outVal);
        bool get_(int vertexId, pointer_type outPtr = nullptr);
        void set_(int vertexId, value_type const& val);
        bool setFromVertex_(Point2f const& pt, int vertexId, Vertices& outVertices);
        bool setFromEdge_(int edgeId, Vertices& outVertices);
        bool setFromFacet_(Point2f const& pt, int edgeId, Vertices& outVertices);
        Subdiv2D subdiv_;
        std::vector<value_type> associatedValues_;
    };

    template <typename T> inline SubdivContainer<T>::SubdivContainer(Rect bounds) : subdiv_(bounds) {}
    template <typename T> inline void SubdivContainer<T>::insert(Point2f const& pt, value_type const& val) {
        auto ptId = subdiv_.insert(pt);
        set_(ptId, val);
    }
    template <typename T> inline bool SubdivContainer<T>::lookup(Point2f const& pt, value_type& outVal) {
        return lookup(pt, &outVal);
    }
    template <typename T> inline bool SubdivContainer<T>::lookup(Point2f const& pt, pointer_type outPtr) {
        int edge = -1;
        int vertex = -1;
        const auto ret = subdiv_.locate(pt, edge, vertex);
        if (ret == PtLoc::PTLOC_VERTEX) {
            return get_(vertex, outPtr);
        }
        return false;
    }
    template <typename T> inline typename SubdivContainer<T>::value_type SubdivContainer<T>::get(Point2f const& pt) {
        value_type ret;
        if (!lookup(pt, ret)) {
            throw std::runtime_error("Vertex not found in subdivision");
        }
        return ret;
    }
    template <typename T> inline bool SubdivContainer<T>::findNearest(Point2f const& pt, Vertices& outVertices) {
        auto vertices = subdiv_.locateVertexIdsArray(pt);
        std::cout << "Vertices ";
        for (auto& v : vertices) {
            std::cout << " " << v;
        }
        std::cout << std::endl;
        int edge;
        int vertex;
        auto locateResult = subdiv_.locate(pt, edge, vertex);
        switch (locateResult) {
        case PtLoc::PTLOC_INSIDE:
            return setFromFacet_(pt, edge, outVertices);
        case PtLoc::PTLOC_ON_EDGE:
            return setFromEdge_(edge, outVertices);
        case PtLoc::PTLOC_VERTEX:
            return setFromVertex_(pt, vertex, outVertices);
        default:
            return false;
        }
    }
    template <typename T> inline bool SubdivContainer<T>::get_(int vertexId, value_type& outVal) {
        return get_(vertexId, &outVal);
    }
    template <typename T> inline bool SubdivContainer<T>::get_(int vertexId, pointer_type outPtr) {
        if (vertexId < NumDummyVertices) {
            return false;
        }
        auto idx = static_cast<std::size_t>(vertexId) - NumDummyVertices;
        if (idx >= associatedValues_.size()) {
            // a vertex we don't possibly have a value for.
            return false;
        }
        if (outPtr) {
            *outPtr = associatedValues_[idx];
        }
        return true;
    }
    template <typename T> inline void SubdivContainer<T>::set_(int vertexId, value_type const& val) {
        std::size_t idx = static_cast<std::size_t>(vertexId) - NumDummyVertices;

        if (associatedValues_.size() <= idx) {
            associatedValues_.resize(idx + 1);
        }
        associatedValues_[idx] = val;
    }
    template <typename T>
    inline bool SubdivContainer<T>::setFromVertex_(Point2f const& pt, int vertexId, Vertices& outVertices) {

        Vertex v;
        v.location = pt;
        if (!get_(vertexId, v.value)) {
            return false;
        }
        outVertices = {std::move(v)};
        return true;
    }
    template <typename T> inline bool SubdivContainer<T>::setFromEdge_(int edgeId, Vertices& outVertices) {
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
    inline bool SubdivContainer<T>::setFromFacet_(Point2f const& pt, int edgeId, Vertices& outVertices) {
        /// @todo
        return false;
    }
} // namespace subdiv2d
} // namespace sensics
#endif // INCLUDED_SubdivContainer_h_GUID_B7F5EEE6_E316_4EA6_5A20_EDADEE3D99A9

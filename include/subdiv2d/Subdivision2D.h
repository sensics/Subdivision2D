/** @file
    @brief Header (based on portions of imgproc.hpp from OpenCV) defining a 2D subdivision class.

    @date 2017

    OpenCV web site clarifies that OpenCV is provided under the 3-clause BSD license.

    SPDX-License-Identifier:BSD-3-Clause
*/

// Copyright 2017 Sensics, Inc.
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.

/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#ifndef INCLUDED_Subdivision2D_h_GUID_469ED0D3_51F3_4C94_4A97_52C530C16C2A
#define INCLUDED_Subdivision2D_h_GUID_469ED0D3_51F3_4C94_4A97_52C530C16C2A

// Internal Includes
#include "Types.h"

// Library/third-party includes
// - none

// Standard includes
#include <array>
#include <vector>

namespace sensics {
namespace subdiv2d {
    /** Subdiv2D point location cases */
    enum class PtLoc {
        PTLOC_ERROR = -2,        //!< Point location error
        PTLOC_OUTSIDE_RECT = -1, //!< Point outside the subdivision bounding rect
        PTLOC_INSIDE = 0,        //!< Point inside some facet
        PTLOC_VERTEX = 1,        //!< Point coincides with one of the subdivision vertices
        PTLOC_ON_EDGE = 2        //!< Point on some edge
    };

    using VertexId = int;
    static const VertexId InvalidVertex = 0;

    using EdgeId = int;
    static const EdgeId InvalidEdge = 0;

    namespace detail {

        struct LocateSubResults {
            using VertexArray = std::array<VertexId, 3>;

            PtLoc locateStatus = PtLoc::PTLOC_ERROR;
            bool refined = false;

            EdgeId getEdge() const;
            bool hasOtherEdge() const;
            EdgeId getOtherEdge() const;
            void setEdges(EdgeId e1 = InvalidEdge, EdgeId e2 = InvalidEdge);

            std::size_t numVertices() const;
            VertexArray const& getVertices() const { return vertices; }
            void addVertex(VertexId vertex);
            void setVertices(std::initializer_list<VertexId> const& newVertices);
            bool isVertexInVertices(VertexId vertex) const;

          private:
            EdgeId edge = InvalidEdge;
            EdgeId otherEdge = InvalidEdge;
            VertexArray vertices = {{InvalidVertex, InvalidVertex, InvalidVertex}};
        };

        class EdgeIterationHelper;
    } // namespace detail

    /**
    The Subdiv2D class described in this section is used to perform various planar subdivision on
    a set of 2D points (represented as vector of Point2f). OpenCV subdivides a plane into triangles
    using the Delaunay's algorithm, which corresponds to the dual graph of the Voronoi diagram.
    In the figure below, the Delaunay's triangulation is marked with black lines and the Voronoi
    diagram with red lines.

    ![Delaunay triangulation (black) and Voronoi (red)](pics/delaunay_voronoi.png)

    The subdivisions can be used for the 3D piece-wise transformation of a plane, morphing, fast
    location of points on the plane, building special graphs (such as NNG,RNG), and so forth.
    */
    class Subdiv2D {
      public:
        using Point = Point2f;
        using value_type = Point::value_type;

        /** Subdiv2D edge type navigation (see: getEdge()) */
        enum {
            NEXT_AROUND_ORG = 0x00,
            NEXT_AROUND_DST = 0x22,
            PREV_AROUND_ORG = 0x11,
            PREV_AROUND_DST = 0x33,
            NEXT_AROUND_LEFT = 0x13,
            NEXT_AROUND_RIGHT = 0x31,
            PREV_AROUND_LEFT = 0x20,
            PREV_AROUND_RIGHT = 0x02
        };

        /** creates an empty Subdiv2D object.
        To create a new empty Delaunay subdivision you need to use the initDelaunay() function.
         */
        Subdiv2D();

        /** @overload

        @param rect Rectangle that includes all of the 2D points that are to be added to the subdivision.

        The function creates an empty Delaunay subdivision where 2D points can be added using the function
        insert() . All of the points to be added must be within the specified rectangle, otherwise a runtime
        error is raised.
         */
        Subdiv2D(Rect rect);

        /** @brief Creates a new empty Delaunay subdivision

        @param rect Rectangle that includes all of the 2D points that are to be added to the subdivision.

         */
        void initDelaunay(Rect rect);

        /** @brief Insert a single point into a Delaunay triangulation.

        @param pt Point to insert.

        The function inserts a single point into a subdivision and modifies the subdivision topology
        appropriately. If a point with the same coordinates exists already, no new point is added.
        @returns the ID of the point.

        @note If the point is outside of the triangulation specified rect a runtime error is raised.
         */
        VertexId insert(Point2f pt);

        /** @brief Insert multiple points into a Delaunay triangulation.

        @param ptvec Points to insert.

        The function inserts a vector of points into a subdivision and modifies the subdivision topology
        appropriately.
         */
        void insert(const std::vector<Point2f>& ptvec);

        /** @brief Returns the location of a point within a Delaunay triangulation.

        @param pt Point to locate.
        @param edge Output edge that the point belongs to or is located to the right of it.
        @param vertex Optional output vertex the input point coincides with.

        The function locates the input point within the subdivision and gives one of the triangle edges
        or vertices.

        @returns an integer which specify one of the following five cases for point location:
        -  The point falls into some facet. The function returns PTLOC_INSIDE and edge will contain one of
           edges of the facet.
        -  The point falls onto the edge. The function returns PTLOC_ON_EDGE and edge will contain this edge.
        -  The point coincides with one of the subdivision vertices. The function returns PTLOC_VERTEX and
           vertex will contain a pointer to the vertex.
        -  The point is outside the subdivision reference rectangle. The function returns PTLOC_OUTSIDE_RECT
           and no pointers are filled.
        -  One of input arguments is invalid. A runtime error is raised or, if silent or "parent" error
           processing mode is selected, CV_PTLOC_ERROR is returned.
         */
        PtLoc locate(Point2f pt, EdgeId& edge, VertexId& vertex);

        /** @brief Finds the subdivision vertex closest to the given point.

        @param pt Input point.
        @param nearestPt Output subdivision vertex point.

        The function is another function that locates the input point within the subdivision. It finds the
        subdivision vertex that is the closest to the input point. It is not necessarily one of vertices
        of the facet containing the input point, though the facet (located using locate() ) is used as a
        starting point.

        @returns vertex ID.
         */
        VertexId findNearest(Point2f pt, Point2f* nearestPt = nullptr);

        /** @brief Gets the number of vertices, including virtual ones, dummy ones, and the placeholder. */
        std::size_t getNumVertices() const { return vtx.size(); }

        struct Edge {
            Point2f origin;
            Point2f destination;
        };

        /** @brief Returns a list of all edges.

        @param edgeList Output vector.
         */
        void getEdgeList(std::vector<Edge>& edgeList) const;

        /** @brief Returns a list of all edges. */
        std::vector<Edge> getEdgeList() const {
            std::vector<Edge> ret;
            getEdgeList(ret);
            return ret;
        }

        /** @brief Returns a list of the leading edge ID connected to each triangle.

        @param leadingEdgeList Output vector.

        The function gives one edge ID for each triangle.
         */
        void getLeadingEdgeList(std::vector<EdgeId>& leadingEdgeList) const;

        using Triangle = std::array<Point2f, 3>;

        /** @brief Returns a list of all triangles.

        @param triangleList Output vector.
         */
        void getTriangleList(std::vector<Triangle>& triangleList) const;

        /** @brief Returns a list of all Voroni facets.

        @param idx Vector of vertices IDs to consider. For all vertices you can pass empty vector.
        @param facetList Output vector of the Voroni facets.
        @param facetCenters Output vector of the Voroni facets center points.

         */
        void getVoronoiFacetList(const std::vector<VertexId>& idx, std::vector<std::vector<Point2f> >& facetList,
                                 std::vector<Point2f>& facetCenters);

        /** @brief Returns vertex location from vertex ID.

        @param vertex vertex ID.
        @param firstEdge Optional. The first edge ID which is connected to the vertex.
        @returns vertex (x,y)

         */
        Point2f getVertex(VertexId vertex, EdgeId* firstEdge = 0) const;

        /** @brief Returns one of the edges related to the given edge.

        @param edge Subdivision edge ID.
        @param nextEdgeType Parameter specifying which of the related edges to return.
        The following values are possible:
        -   NEXT_AROUND_ORG next around the edge origin ( eOnext on the picture below if e is the input edge)
        -   NEXT_AROUND_DST next around the edge vertex ( eDnext )
        -   PREV_AROUND_ORG previous around the edge origin (reversed eRnext )
        -   PREV_AROUND_DST previous around the edge destination (reversed eLnext )
        -   NEXT_AROUND_LEFT next around the left facet ( eLnext )
        -   NEXT_AROUND_RIGHT next around the right facet ( eRnext )
        -   PREV_AROUND_LEFT previous around the left facet (reversed eOnext )
        -   PREV_AROUND_RIGHT previous around the right facet (reversed eDnext )

        ![sample output](pics/quadedge.png)

        @returns edge ID related to the input edge.
         */
        EdgeId getEdge(EdgeId edge, int nextEdgeType) const;

        /** @brief Returns next edge around the edge origin.

        @param edge Subdivision edge ID.

        @returns an integer which is next edge ID around the edge origin: eOnext on the
        picture above if e is the input edge).
         */
        EdgeId nextEdge(EdgeId edge) const;

        /** @brief Returns another edge of the same quad-edge.

        @param edge Subdivision edge ID.
        @param rotate Parameter specifying which of the edges of the same quad-edge as the input
        one to return. The following values are possible:
        -   0 - the input edge ( e on the picture below if e is the input edge)
        -   1 - the rotated edge ( eRot )
        -   2 - the reversed edge (reversed e (in green))
        -   3 - the reversed rotated edge (reversed eRot (in green))

        @returns one of the edges ID of the same quad-edge as the input edge.
         */
        EdgeId rotateEdge(EdgeId edge, int rotate) const;
        EdgeId symEdge(EdgeId edge) const;

        /** @brief Returns the edge origin.

        @param edge Subdivision edge ID.
        @param orgpt Output vertex location.

        @returns vertex ID.
         */
        VertexId edgeOrg(EdgeId edge, Point2f* orgpt = 0) const;

        /** @brief Returns the edge destination.

        @param edge Subdivision edge ID.
        @param dstpt Output vertex location.

        @returns vertex ID.
         */
        VertexId edgeDst(EdgeId edge, Point2f* dstpt = 0) const;

        /** @brief Returns the applicable vertex or vertices (1 if on a vertex, 2 if on an edge, 3 if in a facet) for a
         * given point */
        std::vector<VertexId> locateVertices(Point2f const& pt);

        static constexpr value_type MAX_VAL() { return std::numeric_limits<value_type>::max(); }
        static constexpr value_type EPSILON() { return std::numeric_limits<value_type>::epsilon(); }

      private:
        EdgeId newEdge();
        void deleteEdge(EdgeId edge);
        VertexId newPoint(Point2f pt, bool isvirtual, EdgeId firstEdge = 0);
        void deletePoint(VertexId vtx);
        void setEdgePoints(EdgeId edge, VertexId orgPt, VertexId dstPt);
        void splice(EdgeId edgeA, EdgeId edgeB);
        int connectEdges(EdgeId edgeA, EdgeId edgeB);
        void swapEdges(EdgeId edge);
        int isRightOf(Point2f pt, EdgeId edge) const;
        void calcVoronoi();
        void clearVoronoi();
        void checkSubdiv() const;

        detail::LocateSubResults locateSub(Point2f const& pt);
        void locateRefine(Point2f const& pt, detail::LocateSubResults& result);

        struct Vertex {
            Vertex();
            Vertex(Point2f pt, bool _isvirtual, EdgeId _firstEdge = 0);
            bool isvirtual() const;
            bool isfree() const;

            EdgeId firstEdge;
            int type;
            Point2f pt;
        };

        struct QuadEdge {
            QuadEdge();
            QuadEdge(EdgeId edgeidx);
            bool isfree() const;

            int next[4];
            VertexId pt[4];
        };

        //! All of the vertices
        std::vector<Vertex> vtx;
        //! All of the edges
        std::vector<QuadEdge> qedges;
        int freeQEdge = 0;
        VertexId freePoint = 0;
        bool validGeometry = false;

        EdgeId recentEdge = 0;
        //! Top left corner of the bounding rect
        Point2f topLeft;
        //! Bottom right corner of the bounding rect
        Point2f bottomRight;

        friend class detail::EdgeIterationHelper;
    };

} // namespace subdiv2d
} // namespace sensics

#endif // INCLUDED_Subdivision2D_h_GUID_469ED0D3_51F3_4C94_4A97_52C530C16C2A

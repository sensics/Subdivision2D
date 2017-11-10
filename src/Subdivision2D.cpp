/** @file
    @brief Implementation (based on subdivision2d.cpp in OpenCV imgproc) of a 2D subdivision class.

    @date 2017

    OpenCV web site clarifies that OpenCV is provided under the 3-clause BSD license.

    SPDX-License-Identifier:BSD-3-Clause
*/

// Copyright 2017 Sensics, Inc.
// Copyright (C) 2000, Intel Corporation, all rights reserved.
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
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
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
//   * The name of Intel Corporation may not be used to endorse or promote products
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

// Internal Includes
#include "subdiv2d/Subdivision2D.h"
#include "subdiv2d/AssertAndError.h"

// Library/third-party includes
// - none

// Standard includes
#include <cassert>
#include <iostream>

namespace sensics {
namespace subdiv2d {
    namespace detail {
        EdgeId LocateSubResults::getEdge() const { return edge; }
        bool LocateSubResults::hasOtherEdge() const { return otherEdge != InvalidEdge; }
        EdgeId LocateSubResults::getOtherEdge() const { return otherEdge; }
        void LocateSubResults::setEdges(EdgeId e1, EdgeId e2) {
            Subdiv2D_Assert(!(e1 == InvalidEdge && e2 != InvalidEdge));
            edge = e1;
            otherEdge = e2;
        }
        std::size_t LocateSubResults::numVertices() const {
            return vertices.size() - std::count(vertices.begin(), vertices.end(), InvalidVertex);
        }
        void LocateSubResults::addVertex(VertexId vertex) {
            if (vertex == InvalidVertex) {
                return;
            }
            for (auto& v : vertices) {
                if (v == vertex) {
                    // won't add a duplicate.
                    return;
                }
                if (v == InvalidVertex) {
                    // found a blank space.
                    v = vertex;
                    return;
                }
            }
        }
        void LocateSubResults::setVertices(std::initializer_list<VertexId> const& newVertices) {
            vertices = {{InvalidVertex, InvalidVertex, InvalidVertex}};
            for (auto v : newVertices) {
                addVertex(v);
            }
        }
        bool LocateSubResults::isVertexInVertices(VertexId vertex) const {
            return std::find(vertices.begin(), vertices.end(), vertex) != vertices.end();
        }

        class EdgeIterationHelper {
          public:
            EdgeIterationHelper(std::vector<Subdiv2D::QuadEdge> const& qedges, std::size_t start = 4,
                                std::size_t inc = 2);
            /// Returns true while get() is valid
            explicit operator bool() const;
            /// Increments until a non-visted edge is reached or get would no longer be valid.
            EdgeIterationHelper& advance();

            /// Mark additional edge as visited.
            void markVisited(EdgeId edge);

            /// Get the current edge.
            EdgeId get() const;

          private:
            std::size_t i_ = 4;
            const std::size_t inc_;
            const std::size_t n_;
            std::vector<bool> edgemask_;
        };

        EdgeIterationHelper::EdgeIterationHelper(std::vector<Subdiv2D::QuadEdge> const& qedges, std::size_t start,
                                                 std::size_t inc)
            : i_(start), inc_(inc), n_(qedges.size() * 4), edgemask_(n_, false) {}

        EdgeIterationHelper::operator bool() const { return i_ < n_; }

        EdgeIterationHelper& EdgeIterationHelper::advance() {
            i_ += inc_;
            while (*this) {
                if (!edgemask_[i_]) {
                    break;
                }
                i_ += inc_;
            }
            if (*this) {
                markVisited(get());
            }
            return *this;
        }

        void EdgeIterationHelper::markVisited(EdgeId edge) { edgemask_[edge] = true; }

        EdgeId EdgeIterationHelper::get() const { return EdgeId(i_); }
    } // namespace detail

    EdgeId Subdiv2D::nextEdge(EdgeId edge) const {
        Subdiv2D_DbgAssert((size_t)(edge >> 2) < qedges.size());
        return qedges[edge >> 2].next[edge & 3];
    }

    EdgeId Subdiv2D::rotateEdge(EdgeId edge, int rotate) const { return (edge & ~3) + ((edge + rotate) & 3); }

    EdgeId Subdiv2D::symEdge(EdgeId edge) const { return edge ^ 2; }

    EdgeId Subdiv2D::getEdge(EdgeId edge, int nextEdgeType) const {
        Subdiv2D_DbgAssert((size_t)(edge >> 2) < qedges.size());
        edge = qedges[edge >> 2].next[(edge + nextEdgeType) & 3];
        return (edge & ~3) + ((edge + (nextEdgeType >> 4)) & 3);
    }

    VertexId Subdiv2D::edgeOrg(EdgeId edge, Point2f* orgpt) const {
        Subdiv2D_DbgAssert((size_t)(edge >> 2) < qedges.size());
        VertexId vidx = qedges[edge >> 2].pt[edge & 3];
        if (orgpt) {
            *orgpt = getVertex(vidx);
        }
        return vidx;
    }

    VertexId Subdiv2D::edgeDst(EdgeId edge, Point2f* dstpt) const {
        Subdiv2D_DbgAssert((size_t)(edge >> 2) < qedges.size());
        VertexId vidx = qedges[edge >> 2].pt[(edge + 2) & 3];
        if (dstpt) {
            *dstpt = getVertex(vidx);
        }
        return vidx;
    }

    Point2f Subdiv2D::getVertex(VertexId vertex, EdgeId* firstEdge) const {
        Subdiv2D_DbgAssert((size_t)vertex < vtx.size());
        if (firstEdge) {
            *firstEdge = vtx[vertex].firstEdge;
        }
        return vtx[vertex].pt;
    }

    Subdiv2D::Subdiv2D() {}

    Subdiv2D::Subdiv2D(Rect rect) { initDelaunay(rect); }

    Subdiv2D::QuadEdge::QuadEdge() {
        next[0] = next[1] = next[2] = next[3] = InvalidEdge;
        pt[0] = pt[1] = pt[2] = pt[3] = InvalidVertex;
    }

    Subdiv2D::QuadEdge::QuadEdge(EdgeId edgeidx) {
        Subdiv2D_DbgAssert((edgeidx & 3) == 0);
        next[0] = edgeidx;
        next[1] = edgeidx + 3;
        next[2] = edgeidx + 2;
        next[3] = edgeidx + 1;

        pt[0] = pt[1] = pt[2] = pt[3] = InvalidVertex;
    }

    bool Subdiv2D::QuadEdge::isfree() const { return next[0] <= InvalidEdge; }

    Subdiv2D::Vertex::Vertex() {
        firstEdge = InvalidEdge;
        type = -1;
    }

    Subdiv2D::Vertex::Vertex(Point2f _pt, bool _isvirtual, EdgeId _firstEdge) {
        firstEdge = _firstEdge;
        type = (int)_isvirtual;
        pt = _pt;
    }

    bool Subdiv2D::Vertex::isvirtual() const { return type > 0; }

    bool Subdiv2D::Vertex::isfree() const { return type < 0; }

    void Subdiv2D::splice(EdgeId edgeA, EdgeId edgeB) {
        auto& a_next = qedges[edgeA >> 2].next[edgeA & 3];
        auto& b_next = qedges[edgeB >> 2].next[edgeB & 3];
        auto a_rot = rotateEdge(a_next, 1);
        auto b_rot = rotateEdge(b_next, 1);
        auto& a_rot_next = qedges[a_rot >> 2].next[a_rot & 3];
        auto& b_rot_next = qedges[b_rot >> 2].next[b_rot & 3];
        std::swap(a_next, b_next);
        std::swap(a_rot_next, b_rot_next);
    }

    void Subdiv2D::setEdgePoints(EdgeId edge, VertexId orgPt, VertexId dstPt) {
        qedges[edge >> 2].pt[edge & 3] = orgPt;
        qedges[edge >> 2].pt[(edge + 2) & 3] = dstPt;
        vtx[orgPt].firstEdge = edge;
        vtx[dstPt].firstEdge = edge ^ 2;
    }

    EdgeId Subdiv2D::connectEdges(EdgeId edgeA, EdgeId edgeB) {
        auto edge = newEdge();

        splice(edge, getEdge(edgeA, NEXT_AROUND_LEFT));
        splice(symEdge(edge), edgeB);

        setEdgePoints(edge, edgeDst(edgeA), edgeOrg(edgeB));
        return edge;
    }

    void Subdiv2D::swapEdges(EdgeId edge) {
        auto sedge = symEdge(edge);
        auto a = getEdge(edge, PREV_AROUND_ORG);
        auto b = getEdge(sedge, PREV_AROUND_ORG);

        splice(edge, a);
        splice(sedge, b);

        setEdgePoints(edge, edgeDst(a), edgeDst(b));

        splice(edge, getEdge(a, NEXT_AROUND_LEFT));
        splice(sedge, getEdge(b, NEXT_AROUND_LEFT));
    }

    static double triangleArea(Point2f a, Point2f b, Point2f c) {
        return ((double)b.x - a.x) * ((double)c.y - a.y) - ((double)b.y - a.y) * ((double)c.x - a.x);
    }

    int Subdiv2D::isRightOf(Point2f pt, EdgeId edge) const {
        Point2f org, dst;
        edgeOrg(edge, &org);
        edgeDst(edge, &dst);
        double cw_area = triangleArea(pt, dst, org);

        return (cw_area > 0) - (cw_area < 0);
    }

    std::vector<VertexId> Subdiv2D::locateVertices(Point2f const& pt) {
        auto result = locateSub(pt);
        locateRefine(pt, result);
        std::vector<VertexId> ret;
        for (auto v : result.getVertices()) {
            if (v != InvalidVertex) {
                ret.push_back(v);
            }
        }
        return ret;
    }

    EdgeId Subdiv2D::newEdge() {
        if (freeQEdge <= 0) {
            qedges.push_back(QuadEdge());
            freeQEdge = (int)(qedges.size() - 1);
        }
        EdgeId edge = freeQEdge * 4;
        freeQEdge = qedges[edge >> 2].next[1];
        qedges[edge >> 2] = QuadEdge(edge);
        return edge;
    }

    void Subdiv2D::deleteEdge(EdgeId edge) {
        Subdiv2D_DbgAssert((size_t)(edge >> 2) < (size_t)qedges.size());
        splice(edge, getEdge(edge, PREV_AROUND_ORG));
        auto sedge = symEdge(edge);
        splice(sedge, getEdge(sedge, PREV_AROUND_ORG));

        edge >>= 2;
        qedges[edge].next[0] = InvalidEdge;
        qedges[edge].next[1] = freeQEdge;
        freeQEdge = edge;
    }

    VertexId Subdiv2D::newPoint(Point2f pt, bool isvirtual, EdgeId firstEdge) {
        if (freePoint == InvalidVertex) {
            vtx.push_back(Vertex());
            freePoint = static_cast<VertexId>(vtx.size() - 1);
        }
        VertexId vidx = freePoint;
        freePoint = vtx[vidx].firstEdge;
        vtx[vidx] = Vertex(pt, isvirtual, firstEdge);

        return vidx;
    }

    void Subdiv2D::deletePoint(VertexId vidx) {
        Subdiv2D_DbgAssert(static_cast<size_t>(vidx) < vtx.size());
        vtx[vidx].firstEdge = freePoint;
        vtx[vidx].type = -1;
        freePoint = vidx;
    }

    PtLoc Subdiv2D::locate(Point2f pt, EdgeId& _edge, VertexId& _vertex) {
        auto result = locateSub(pt);
        locateRefine(pt, result);

        _edge = result.getEdge();
        if (result.numVertices() == 1) {
            _vertex = result.getVertices().front();
        } else {
            _vertex = InvalidVertex;
        }

        return result.locateStatus;
    }

    inline int isPtInCircle3(Point2f pt, Point2f a, Point2f b, Point2f c) {
        const double eps = Subdiv2D::EPSILON() * 0.125;
        double val = ((double)a.x * a.x + (double)a.y * a.y) * triangleArea(b, c, pt);
        val -= ((double)b.x * b.x + (double)b.y * b.y) * triangleArea(a, c, pt);
        val += ((double)c.x * c.x + (double)c.y * c.y) * triangleArea(a, b, pt);
        val -= ((double)pt.x * pt.x + (double)pt.y * pt.y) * triangleArea(a, b, c);

        return val > eps ? 1 : val < -eps ? -1 : 0;
    }

    VertexId Subdiv2D::insert(Point2f pt) {

        VertexId curr_point = InvalidVertex;
        EdgeId curr_edge = InvalidEdge;
        auto location = locate(pt, curr_edge, curr_point);

        if (location == PtLoc::PTLOC_ERROR) {
            Subdiv2D_Error(Error::StsBadSize, "");
        }

        if (location == PtLoc::PTLOC_OUTSIDE_RECT) {
            Subdiv2D_Error(Error::StsOutOfRange, "");
        }

        if (location == PtLoc::PTLOC_VERTEX) {
            return curr_point;
        }

        if (location == PtLoc::PTLOC_ON_EDGE) {
            auto deleted_edge = curr_edge;
            recentEdge = curr_edge = getEdge(curr_edge, PREV_AROUND_ORG);
            deleteEdge(deleted_edge);
        } else if (location == PtLoc::PTLOC_INSIDE) {
            // pass - this is fine
        } else {
            Subdiv2D_Error_(Error::StsError,
                            ("Subdiv2D::locate returned invalid location = %d", static_cast<int>(location)));
        }

        assert(curr_edge != InvalidEdge);
        validGeometry = false;

        curr_point = newPoint(pt, false);
        auto base_edge = newEdge();
        auto first_point = edgeOrg(curr_edge);
        setEdgePoints(base_edge, first_point, curr_point);
        splice(base_edge, curr_edge);

        do {
            base_edge = connectEdges(curr_edge, symEdge(base_edge));
            curr_edge = getEdge(base_edge, PREV_AROUND_ORG);
        } while (edgeDst(curr_edge) != first_point);

        curr_edge = getEdge(base_edge, PREV_AROUND_ORG);

        const auto max_edges = qedges.size() * 4;
        for (std::size_t i = 0; i < max_edges; ++i) {
            auto temp_edge = getEdge(curr_edge, PREV_AROUND_ORG);

            auto temp_dst = edgeDst(temp_edge);
            auto curr_org = edgeOrg(curr_edge);
            auto curr_dst = edgeDst(curr_edge);

            if (isRightOf(getVertex(temp_dst), curr_edge) > 0 &&
                isPtInCircle3(getVertex(curr_org), getVertex(temp_dst), getVertex(curr_dst), getVertex(curr_point)) <
                    0) {
                swapEdges(curr_edge);
                curr_edge = getEdge(curr_edge, PREV_AROUND_ORG);
            } else if (curr_org == first_point) {
                break;
            } else {
                curr_edge = getEdge(nextEdge(curr_edge), PREV_AROUND_LEFT);
            }
        }

        return curr_point;
    }

    void Subdiv2D::insert(const std::vector<Point2f>& ptvec) {
        for (auto& pt : ptvec) {
            insert(pt);
        }
    }

    void Subdiv2D::initDelaunay(Rect rect) {

        float big_coord = 3.f * std::max(rect.width, rect.height);
        float rx = (float)rect.x;
        float ry = (float)rect.y;

        vtx.clear();
        qedges.clear();

        recentEdge = InvalidEdge;
        validGeometry = false;

        topLeft = Point2f(rx, ry);
        bottomRight = Point2f(rx + rect.width, ry + rect.height);

        Point2f ppA(rx + big_coord, ry);
        Point2f ppB(rx, ry + big_coord);
        Point2f ppC(rx - big_coord, ry - big_coord);

        // Vertex 0: null/dummy - 0 is an invalid vertex ID
        vtx.push_back(Vertex());
        qedges.push_back(QuadEdge());

        freeQEdge = InvalidEdge;
        freePoint = InvalidVertex;

        // Vertex 1: top, way past the right edge.
        int pA = newPoint(ppA, false);
        // Vertex 2: left, way past the bottom edge.
        int pB = newPoint(ppB, false);
        // Vertex 3: Way past the top-left corner
        int pC = newPoint(ppC, false);

        int edge_AB = newEdge();
        int edge_BC = newEdge();
        int edge_CA = newEdge();

        setEdgePoints(edge_AB, pA, pB);
        setEdgePoints(edge_BC, pB, pC);
        setEdgePoints(edge_CA, pC, pA);

        splice(edge_AB, symEdge(edge_CA));
        splice(edge_BC, symEdge(edge_AB));
        splice(edge_CA, symEdge(edge_BC));

        recentEdge = edge_AB;
    }

    void Subdiv2D::clearVoronoi() {

        for (auto& qedge : qedges) {
            qedge.pt[1] = InvalidVertex;
            qedge.pt[3] = InvalidVertex;
        }

        const auto total = vtx.size();
        for (std::size_t i = 0; i < total; ++i) {
            if (vtx[i].isvirtual())
                deletePoint(static_cast<VertexId>(i));
        }

        validGeometry = false;
    }

    static Point2f computeVoronoiPoint(Point2f org0, Point2f dst0, Point2f org1, Point2f dst1) {
        double a0 = dst0.x - org0.x;
        double b0 = dst0.y - org0.y;
        double c0 = -0.5 * (a0 * (dst0.x + org0.x) + b0 * (dst0.y + org0.y));

        double a1 = dst1.x - org1.x;
        double b1 = dst1.y - org1.y;
        double c1 = -0.5 * (a1 * (dst1.x + org1.x) + b1 * (dst1.y + org1.y));

        double det = a0 * b1 - a1 * b0;

        if (det != 0) {
            det = 1. / det;
            return Point2f((float)((b0 * c1 - b1 * c0) * det), (float)((a1 * c0 - a0 * c1) * det));
        }

        return Point2f(Subdiv2D::MAX_VAL(), Subdiv2D::MAX_VAL());
    }

    void Subdiv2D::calcVoronoi() {
        // check if it is already calculated
        if (validGeometry)
            return;

        clearVoronoi();
        // loop through all quad-edges, except for the first 3 (#1, #2, #3 - 0 is reserved for "NULL" pointer)
        const auto total = qedges.size();
        for (std::size_t i = 4; i < total; ++i) {
            QuadEdge& quadedge = qedges[i];

            if (quadedge.isfree()) {
                continue;
            }

            auto edge0 = static_cast<EdgeId>(i * 4);
            Point2f org0, dst0, org1, dst1;

            if (quadedge.pt[3] == InvalidVertex) {
                auto edge1 = getEdge(edge0, NEXT_AROUND_LEFT);
                auto edge2 = getEdge(edge1, NEXT_AROUND_LEFT);

                edgeOrg(edge0, &org0);
                edgeDst(edge0, &dst0);
                edgeOrg(edge1, &org1);
                edgeDst(edge1, &dst1);

                Point2f virt_point = computeVoronoiPoint(org0, dst0, org1, dst1);

                if (std::abs(virt_point.x) < MAX_VAL() * 0.5 && std::abs(virt_point.y) < MAX_VAL() * 0.5) {
                    quadedge.pt[3] = qedges[edge1 >> 2].pt[3 - (edge1 & 2)] = qedges[edge2 >> 2].pt[3 - (edge2 & 2)] =
                        newPoint(virt_point, true);
                }
            }

            if (!quadedge.pt[1] == InvalidVertex) {
                auto edge1 = getEdge(edge0, NEXT_AROUND_RIGHT);
                auto edge2 = getEdge(edge1, NEXT_AROUND_RIGHT);

                edgeOrg(edge0, &org0);
                edgeDst(edge0, &dst0);
                edgeOrg(edge1, &org1);
                edgeDst(edge1, &dst1);

                Point2f virt_point = computeVoronoiPoint(org0, dst0, org1, dst1);

                if (std::abs(virt_point.x) < MAX_VAL() * 0.5 && std::abs(virt_point.y) < MAX_VAL() * 0.5) {
                    quadedge.pt[1] = qedges[edge1 >> 2].pt[1 + (edge1 & 2)] = qedges[edge2 >> 2].pt[1 + (edge2 & 2)] =
                        newPoint(virt_point, true);
                }
            }
        }

        validGeometry = true;
    }

    static int isRightOf2(const Point2f& pt, const Point2f& org, const Point2f& diff) {
        double cw_area = ((double)org.x - pt.x) * diff.y - ((double)org.y - pt.y) * diff.x;
        return (cw_area > 0) - (cw_area < 0);
    }

    VertexId Subdiv2D::findNearest(Point2f pt, Point2f* nearestPt) {

        if (!validGeometry) {
            calcVoronoi();
        }

        VertexId vertex = InvalidVertex;
        EdgeId edge = InvalidEdge;
        auto loc = locate(pt, edge, vertex);

        if (loc != PtLoc::PTLOC_ON_EDGE && loc != PtLoc::PTLOC_INSIDE) {
            return vertex;
        }

        vertex = InvalidVertex;

        Point2f start;
        edgeOrg(edge, &start);
        Point2f diff = pt - start;

        edge = rotateEdge(edge, 1);

        const auto total = vtx.size();
        for (std::size_t i = 0; i < total; ++i) {
            Point2f t;

            for (;;) {
                Subdiv2D_Assert(edgeDst(edge, &t) > 0);
                if (isRightOf2(t, start, diff) >= 0)
                    break;

                edge = getEdge(edge, NEXT_AROUND_LEFT);
            }

            for (;;) {
                Subdiv2D_Assert(edgeOrg(edge, &t) > 0);

                if (isRightOf2(t, start, diff) < 0)
                    break;

                edge = getEdge(edge, PREV_AROUND_LEFT);
            }

            Point2f tempDiff;
            edgeDst(edge, &tempDiff);
            edgeOrg(edge, &t);
            tempDiff -= t;

            if (isRightOf2(pt, t, tempDiff) >= 0) {
                vertex = edgeOrg(rotateEdge(edge, 3));
                break;
            }

            edge = symEdge(edge);
        }

        if (nearestPt && vertex > InvalidVertex) {
            *nearestPt = getVertex(vertex);
        }

        return vertex;
    }

    void Subdiv2D::getEdgeList(std::vector<Edge>& edgeList) const {
        edgeList.clear();
        const auto n = qedges.size();
        for (size_t i = 4; i < n; ++i) {
            if (qedges[i].isfree()) {
                continue;
            }
            const auto& qedge = qedges[i];
            if (qedge.pt[0] > InvalidVertex && qedge.pt[2] > InvalidVertex) {
                Point2f org = getVertex(qedge.pt[0]);
                Point2f dst = getVertex(qedge.pt[2]);
                edgeList.push_back(Edge{org, dst});
            }
        }
    }

    void Subdiv2D::getLeadingEdgeList(std::vector<EdgeId>& leadingEdgeList) const {
        leadingEdgeList.clear();
        for (detail::EdgeIterationHelper helper(qedges, 4, 2); helper; helper.advance()) {
            auto edge = helper.get();
            leadingEdgeList.push_back(edge);
            edge = getEdge(edge, NEXT_AROUND_LEFT);
            helper.markVisited(edge);
            edge = getEdge(edge, NEXT_AROUND_LEFT);
            helper.markVisited(edge);
        }
    }

    void Subdiv2D::getTriangleList(std::vector<Triangle>& triangleList) const {
        triangleList.clear();
        for (detail::EdgeIterationHelper helper(qedges, 4, 2); helper; helper.advance()) {
            auto edge = helper.get();
            Point2f a, b, c;
            edgeOrg(edge, &a);

            edge = getEdge(edge, NEXT_AROUND_LEFT);
            helper.markVisited(edge);
            edgeOrg(edge, &b);

            edge = getEdge(edge, NEXT_AROUND_LEFT);
            helper.markVisited(edge);
            edgeOrg(edge, &c);

            triangleList.push_back(Triangle{{a, b, c}});
        }
    }

    void Subdiv2D::getVoronoiFacetList(const std::vector<VertexId>& idx, std::vector<std::vector<Point2f> >& facetList,
                                       std::vector<Point2f>& facetCenters) {
        calcVoronoi();
        facetList.clear();
        facetCenters.clear();

        std::vector<Point2f> buf;

        size_t i = 0;
        size_t total = idx.size();
        if (idx.empty()) {
            i = 4;
            total = vtx.size();
        }

        for (; i < total; ++i) {
            VertexId k = idx.empty() ? static_cast<VertexId>(i) : idx[i];

            if (vtx[k].isfree() || vtx[k].isvirtual()) {
                continue;
            }
            auto edge = rotateEdge(vtx[k].firstEdge, 1);
            auto t = edge;

            // gather points
            buf.clear();
            do {
                buf.push_back(getVertex(edgeOrg(t)));
                t = getEdge(t, NEXT_AROUND_LEFT);
            } while (t != edge);

            facetList.push_back(buf);
            facetCenters.push_back(getVertex(k));
        }
    }

    void Subdiv2D::checkSubdiv() const {
        int i, j, total = (int)qedges.size();

        for (i = 0; i < total; i++) {
            const QuadEdge& qe = qedges[i];

            if (qe.isfree())
                continue;

            for (j = 0; j < 4; j++) {
                int e = (int)(i * 4 + j);
                int o_next = nextEdge(e);
                int o_prev = getEdge(e, PREV_AROUND_ORG);
                int d_prev = getEdge(e, PREV_AROUND_DST);
                int d_next = getEdge(e, NEXT_AROUND_DST);

                // check points
                Subdiv2D_Assert(edgeOrg(e) == edgeOrg(o_next));
                Subdiv2D_Assert(edgeOrg(e) == edgeOrg(o_prev));
                Subdiv2D_Assert(edgeDst(e) == edgeDst(d_next));
                Subdiv2D_Assert(edgeDst(e) == edgeDst(d_prev));

                if (j % 2 == 0) {
                    Subdiv2D_Assert(edgeDst(o_next) == edgeOrg(d_prev));
                    Subdiv2D_Assert(edgeDst(o_prev) == edgeOrg(d_next));
                    Subdiv2D_Assert(
                        getEdge(getEdge(getEdge(e, NEXT_AROUND_LEFT), NEXT_AROUND_LEFT), NEXT_AROUND_LEFT) == e);
                    Subdiv2D_Assert(
                        getEdge(getEdge(getEdge(e, NEXT_AROUND_RIGHT), NEXT_AROUND_RIGHT), NEXT_AROUND_RIGHT) == e);
                }
            }
        }
    }

    detail::LocateSubResults Subdiv2D::locateSub(Point2f const& pt) {
        if (qedges.size() < 4) {
            Subdiv2D_Error(Error::StsError, "Subdivision is empty");
        }
        if (pt.x < topLeft.x || pt.y < topLeft.y || pt.x >= bottomRight.x || pt.y >= bottomRight.y) {
            Subdiv2D_Error(Error::StsOutOfRange, "");
        }
        const std::size_t maxEdges = qedges.size() * 4;

        detail::LocateSubResults ret;
        int edge = recentEdge;
        Subdiv2D_Assert(edge > 0);

        int right_of_curr = isRightOf(pt, edge);
        if (right_of_curr > 0) {
            edge = symEdge(edge);
            right_of_curr = -right_of_curr;
        }

        for (std::size_t i = 0; i < maxEdges; ++i) {
            int onext_edge = nextEdge(edge);
            int dprev_edge = getEdge(edge, PREV_AROUND_DST);

            std::cout << "edge: " << edge << "\t onext_edge: " << onext_edge << "\t dprev_edge: " << dprev_edge
                      << std::endl;

            int right_of_onext = isRightOf(pt, onext_edge);
            int right_of_dprev = isRightOf(pt, dprev_edge);

            if (right_of_dprev > 0) {
                if (right_of_onext > 0 || (right_of_onext == 0 && right_of_curr == 0)) {
                    ret.locateStatus = PtLoc::PTLOC_INSIDE;
                    ret.setEdges(edge, dprev_edge);
                    break;
                } else {
                    right_of_curr = right_of_onext;
                    edge = onext_edge;
                }
            } else {
                if (right_of_onext > 0) {
                    if (right_of_dprev == 0 && right_of_curr == 0) {
                        ret.locateStatus = PtLoc::PTLOC_INSIDE;
                        ret.setEdges(edge, onext_edge);
                        break;
                    } else {
                        right_of_curr = right_of_dprev;
                        edge = dprev_edge;
                    }
                } else if (right_of_curr == 0 && isRightOf(getVertex(edgeDst(onext_edge)), edge) >= 0) {
                    edge = symEdge(edge);
                } else {
                    right_of_curr = right_of_onext;
                    edge = onext_edge;
                }
            }
        }

        recentEdge = edge;
        return ret;
    }

    static inline float simpleAbsPointDistance(Point2f const& a, Point2f const& b) {
        // think this is the manhattan distance...
        auto diff = a - b;
        return std::abs(diff.x) + std::abs(diff.y);
    }
    void Subdiv2D::locateRefine(Point2f const& pt, detail::LocateSubResults& result) {
        if (result.locateStatus != PtLoc::PTLOC_INSIDE || result.refined) {
            // no further refinement.
            return;
        }
        result.refined = true;

        Point2f org_pt;
        auto orgId = edgeOrg(result.getEdge(), &org_pt);
        auto orgDist = simpleAbsPointDistance(pt, org_pt);
        Point2f dst_pt;
        auto dstId = edgeDst(result.getEdge(), &dst_pt);
        auto dstDist = simpleAbsPointDistance(pt, dst_pt);

        auto edgeDist = simpleAbsPointDistance(org_pt, dst_pt);

        if (orgDist < EPSILON()) {
            result.locateStatus = PtLoc::PTLOC_VERTEX;
            result.setVertices({orgId});
            result.setEdges();
        } else if (dstDist < EPSILON()) {
            result.locateStatus = PtLoc::PTLOC_VERTEX;
            result.setVertices({dstId});
            result.setEdges();
        } else if ((orgDist < edgeDist || dstDist < edgeDist) &&
                   std::abs(triangleArea(pt, org_pt, dst_pt)) < EPSILON()) {
            result.locateStatus = PtLoc::PTLOC_ON_EDGE;
            result.setEdges(result.getEdge());
            result.setVertices({orgId, dstId});
        } else {
            // this case means, really inside.
            result.setVertices({orgId, dstId});
            // Grab the second edge in the vector.
            Subdiv2D_Assert(result.hasOtherEdge());
            auto otherEdge = result.getOtherEdge();
            auto otherOrg = edgeOrg(otherEdge);
            auto otherDst = edgeDst(otherEdge);
            if (!result.isVertexInVertices(otherOrg)) {
                result.addVertex(otherOrg);
            } else if (!result.isVertexInVertices(otherDst)) {
                result.addVertex(otherDst);
            } else {
                Subdiv2D_Error(Error::StsAssert,
                               "Should never happen - our other edge didn't have a useful other vertex.");
            }
        }
    }

} // namespace subdiv2d
} // namespace sensics

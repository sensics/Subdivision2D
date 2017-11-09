/** @file
    @brief Header containing suppor types. Point2 is written from scratch. Rect_ is based on the OpenCV equivalent.

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Copyright (C) 2013, OpenCV Foundation, all rights reserved.
// Copyright 2017 Sensics, Inc.
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
//                          License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Copyright (C) 2013, OpenCV Foundation, all rights reserved.
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

#ifndef INCLUDED_Types_h_GUID_EE6A183B_3E3E_4B8B_28FB_A6114F219F51
#define INCLUDED_Types_h_GUID_EE6A183B_3E3E_4B8B_28FB_A6114F219F51

// Internal Includes
// - none

// Library/third-party includes
// - none

// Standard includes
#include <cfloat>
#include <limits>

namespace sensics {
namespace subdiv2d {
    template <typename T> struct Point2 {
        using value_type = T;
        Point2(value_type xVal, value_type yVal) : x(xVal), y(yVal) {}
        Point2() : x(std::numeric_limits<T>::(max)()), y(std::numeric_limits<T>::(max)()) {}

        value_type x;
        value_type y;

        Point2& operator-=(Point2 const& rhs) {
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }
    };
    template <typename T> inline Point2<T> operator-(Point2<T> const& lhs, Point2<T> const& rhs) {
        auto ret = Point2<T>(lhs);
        ret -= rhs;
        return ret;
    }

    using Point2f = Point2<float>;



    /** @brief Template class for 2D rectangles

    described by the following parameters:
    -   Coordinates of the top-left corner. This is a default interpretation of Rect_::x and Rect_::y
    in OpenCV. Though, in your algorithms you may count x and y from the bottom-left corner.
    -   Rectangle width and height.

    OpenCV typically assumes that the top and left boundary of the rectangle are inclusive, while the
    right and bottom boundaries are not. For example, the method Rect_::contains returns true if

    \f[x  \leq pt.x < x+width,
    y  \leq pt.y < y+height\f]

    Virtually every loop over an image ROI in OpenCV (where ROI is specified by Rect_\<int\> ) is
    implemented as:
    @code
    for(int y = roi.y; y < roi.y + roi.height; y++)
    for(int x = roi.x; x < roi.x + roi.width; x++)
    {
    // ...
    }
    @endcode
    In addition to the class members, the following operations on rectangles are implemented:
    -   \f$\texttt{rect} = \texttt{rect} \pm \texttt{point}\f$ (shifting a rectangle by a certain offset)
    -   \f$\texttt{rect} = \texttt{rect} \pm \texttt{size}\f$ (expanding or shrinking a rectangle by a
    certain amount)
    -   rect += point, rect -= point, rect += size, rect -= size (augmenting operations)
    -   rect = rect1 & rect2 (rectangle intersection)
    -   rect = rect1 | rect2 (minimum area rectangle containing rect1 and rect2 )
    -   rect &= rect1, rect |= rect1 (and the corresponding augmenting operations)
    -   rect == rect1, rect != rect1 (rectangle comparison)

    This is an example how the partial ordering on rectangles can be established (rect1 \f$\subseteq\f$
    rect2):
    @code
    template<typename _Tp> inline bool
    operator <= (const Rect_<_Tp>& r1, const Rect_<_Tp>& r2)
    {
    return (r1 & r2) == r1;
    }
    @endcode
    For your convenience, the Rect_\<\> alias is available: cv::Rect
    */
    template<typename _Tp> class Rect_
    {
    public:
        typedef _Tp value_type;

        //! various constructors
        Rect_();
        Rect_(_Tp _x, _Tp _y, _Tp _width, _Tp _height);
        Rect_(const Rect_& r);
        Rect_(const Point_<_Tp>& org, const Size_<_Tp>& sz);
        Rect_(const Point_<_Tp>& pt1, const Point_<_Tp>& pt2);

        Rect_& operator = (const Rect_& r);
        //! the top-left corner
        Point_<_Tp> tl() const;
        //! the bottom-right corner
        Point_<_Tp> br() const;

        //! size (width, height) of the rectangle
        Size_<_Tp> size() const;
        //! area (width*height) of the rectangle
        _Tp area() const;

        //! conversion to another data type
        template<typename _Tp2> operator Rect_<_Tp2>() const;

        //! checks whether the rectangle contains the point
        bool contains(const Point_<_Tp>& pt) const;

        _Tp x, y, width, height; //< the top-left corner, as well as width and height of the rectangle
    };

    typedef Rect_<int> Rect2i;
    typedef Rect_<float> Rect2f;
    typedef Rect_<double> Rect2d;
    typedef Rect2i Rect;

    template<typename _Tp> inline
        Rect_<_Tp>::Rect_()
        : x(0), y(0), width(0), height(0) {}

    template<typename _Tp> inline
        Rect_<_Tp>::Rect_(_Tp _x, _Tp _y, _Tp _width, _Tp _height)
        : x(_x), y(_y), width(_width), height(_height) {}

    template<typename _Tp> inline
        Rect_<_Tp>::Rect_(const Rect_<_Tp>& r)
        : x(r.x), y(r.y), width(r.width), height(r.height) {}

    template<typename _Tp> inline
        Rect_<_Tp>::Rect_(const Point_<_Tp>& org, const Size_<_Tp>& sz)
        : x(org.x), y(org.y), width(sz.width), height(sz.height) {}

    template<typename _Tp> inline
        Rect_<_Tp>::Rect_(const Point_<_Tp>& pt1, const Point_<_Tp>& pt2)
    {
        x = std::min(pt1.x, pt2.x);
        y = std::min(pt1.y, pt2.y);
        width = std::max(pt1.x, pt2.x) - x;
        height = std::max(pt1.y, pt2.y) - y;
    }

    template<typename _Tp> inline
        Rect_<_Tp>& Rect_<_Tp>::operator = (const Rect_<_Tp>& r)
    {
        x = r.x;
        y = r.y;
        width = r.width;
        height = r.height;
        return *this;
    }

    template<typename _Tp> inline
        Point_<_Tp> Rect_<_Tp>::tl() const
    {
        return Point_<_Tp>(x, y);
    }

    template<typename _Tp> inline
        Point_<_Tp> Rect_<_Tp>::br() const
    {
        return Point_<_Tp>(x + width, y + height);
    }

    template<typename _Tp> inline
        Size_<_Tp> Rect_<_Tp>::size() const
    {
        return Size_<_Tp>(width, height);
    }

    template<typename _Tp> inline
        _Tp Rect_<_Tp>::area() const
    {
        const _Tp result = width * height;
        CV_DbgAssert(!std::numeric_limits<_Tp>::is_integer
            || width == 0 || result / width == height); // make sure the result fits in the return value
        return result;
    }

    template<typename _Tp> template<typename _Tp2> inline
        Rect_<_Tp>::operator Rect_<_Tp2>() const
    {
        return Rect_<_Tp2>(saturate_cast<_Tp2>(x), saturate_cast<_Tp2>(y), saturate_cast<_Tp2>(width), saturate_cast<_Tp2>(height));
    }

    template<typename _Tp> inline
        bool Rect_<_Tp>::contains(const Point_<_Tp>& pt) const
    {
        return x <= pt.x && pt.x < x + width && y <= pt.y && pt.y < y + height;
    }


    template<typename _Tp> static inline
        Rect_<_Tp>& operator += (Rect_<_Tp>& a, const Point_<_Tp>& b)
    {
        a.x += b.x;
        a.y += b.y;
        return a;
    }

    template<typename _Tp> static inline
        Rect_<_Tp>& operator -= (Rect_<_Tp>& a, const Point_<_Tp>& b)
    {
        a.x -= b.x;
        a.y -= b.y;
        return a;
    }

    template<typename _Tp> static inline
        Rect_<_Tp>& operator += (Rect_<_Tp>& a, const Size_<_Tp>& b)
    {
        a.width += b.width;
        a.height += b.height;
        return a;
    }

    template<typename _Tp> static inline
        Rect_<_Tp>& operator -= (Rect_<_Tp>& a, const Size_<_Tp>& b)
    {
        a.width -= b.width;
        a.height -= b.height;
        return a;
    }

    template<typename _Tp> static inline
        Rect_<_Tp>& operator &= (Rect_<_Tp>& a, const Rect_<_Tp>& b)
    {
        _Tp x1 = std::max(a.x, b.x);
        _Tp y1 = std::max(a.y, b.y);
        a.width = std::min(a.x + a.width, b.x + b.width) - x1;
        a.height = std::min(a.y + a.height, b.y + b.height) - y1;
        a.x = x1;
        a.y = y1;
        if (a.width <= 0 || a.height <= 0)
            a = Rect();
        return a;
    }

    template<typename _Tp> static inline
        Rect_<_Tp>& operator |= (Rect_<_Tp>& a, const Rect_<_Tp>& b)
    {
        if (!a.area()) {
            a = b;
        }
        else if (b.area()) {
            _Tp x1 = std::min(a.x, b.x);
            _Tp y1 = std::min(a.y, b.y);
            a.width = std::max(a.x + a.width, b.x + b.width) - x1;
            a.height = std::max(a.y + a.height, b.y + b.height) - y1;
            a.x = x1;
            a.y = y1;
        }
        return a;
    }

    template<typename _Tp> static inline
        bool operator == (const Rect_<_Tp>& a, const Rect_<_Tp>& b)
    {
        return a.x == b.x && a.y == b.y && a.width == b.width && a.height == b.height;
    }

    template<typename _Tp> static inline
        bool operator != (const Rect_<_Tp>& a, const Rect_<_Tp>& b)
    {
        return a.x != b.x || a.y != b.y || a.width != b.width || a.height != b.height;
    }

    template<typename _Tp> static inline
        Rect_<_Tp> operator + (const Rect_<_Tp>& a, const Point_<_Tp>& b)
    {
        return Rect_<_Tp>(a.x + b.x, a.y + b.y, a.width, a.height);
    }

    template<typename _Tp> static inline
        Rect_<_Tp> operator - (const Rect_<_Tp>& a, const Point_<_Tp>& b)
    {
        return Rect_<_Tp>(a.x - b.x, a.y - b.y, a.width, a.height);
    }

    template<typename _Tp> static inline
        Rect_<_Tp> operator + (const Rect_<_Tp>& a, const Size_<_Tp>& b)
    {
        return Rect_<_Tp>(a.x, a.y, a.width + b.width, a.height + b.height);
    }

    template<typename _Tp> static inline
        Rect_<_Tp> operator & (const Rect_<_Tp>& a, const Rect_<_Tp>& b)
    {
        Rect_<_Tp> c = a;
        return c &= b;
    }

    template<typename _Tp> static inline
        Rect_<_Tp> operator | (const Rect_<_Tp>& a, const Rect_<_Tp>& b)
    {
        Rect_<_Tp> c = a;
        return c |= b;
    }
} // namespace subdiv2d
} // namespace sensics
#endif // INCLUDED_Types_h_GUID_EE6A183B_3E3E_4B8B_28FB_A6114F219F51

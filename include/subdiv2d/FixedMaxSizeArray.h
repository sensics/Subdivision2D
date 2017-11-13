/** @file
    @brief Header

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
// SPDX-License-Identifier:BSD-3-Clause

#ifndef INCLUDED_FixedMaxSizeArray_h_GUID_3B5F9209_0B80_4A91_B0D4_5CCAF36EF73A
#define INCLUDED_FixedMaxSizeArray_h_GUID_3B5F9209_0B80_4A91_B0D4_5CCAF36EF73A

// Internal Includes
#include "Subdiv2DConfig.h"

#ifdef SUBDIV2D_USE_BOOST_STATIC_VECTOR

#include <boost/container/static_vector.hpp>
namespace sensics {
namespace subdiv2d {
    template <typename T, std::size_t MaxSize> using MaxSizeVector = boost::container::static_vector<T, MaxSize>;
} // namespace subdiv2d
} // namespace sensics

#else
#error "Haven't written antyhing on this path yet sorry!"
#endif

#endif // INCLUDED_FixedMaxSizeArray_h_GUID_3B5F9209_0B80_4A91_B0D4_5CCAF36EF73A

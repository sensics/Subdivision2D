/** @file
    @brief Implementation file for AssertAndError.h (based on OpenCV core/base.hpp) providing assert and error handling.

    @date 2017
*/

// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Copyright (C) 2013, OpenCV Foundation, all rights reserved.
// Copyright (C) 2014, Itseez Inc., all rights reserved.
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
// Copyright (C) 2014, Itseez Inc., all rights reserved.
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

// Internal Includes
#include "subdiv2d/AssertAndError.h"

// Library/third-party includes
// - none

// Standard includes
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace sensics {
namespace subdiv2d {
    void error(int _code, const std::string& _err, const char* _func, const char* _file, int _line) {
        std::ostringstream os;
        os << "Error: " << _err << " in " << ((_func[0] != '\0') ? _func : "unknown function") << ", file " << _file
           << ", line " << _line;

        std::cerr << os.str() << std::endl;
#ifdef __ANDROID__
        std::string msg = os.str();
        __android_log_print(ANDROID_LOG_ERROR, "subdiv2d::error()", "%s", msg.c_str());
#endif

#ifdef SUBDIV2D_BREAK_ON_ERROR
        static volatile int* p = 0;
        *p = 0;
#endif
        throw std::runtime_error(os.str());
    }
} // namespace subdiv2d
} // namespace sensics

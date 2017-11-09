/** @file
    @brief Header (based on OpenCV core/base.hpp) providing assert and error handling.

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

#ifndef INCLUDED_AssertAndError_h_GUID_DCA9D718_0346_4088_A1AC_7EE0A243A538
#define INCLUDED_AssertAndError_h_GUID_DCA9D718_0346_4088_A1AC_7EE0A243A538

// Internal Includes
// - none

// Library/third-party includes
// - none

// Standard includes
#include <algorithm>
#include <climits>
#include <sstream>
#include <string>

namespace sensics {
namespace subdiv2d {

    namespace Error {
        //! error codes
        enum Code {
            StsOk = 0,             //!< everything is ok
            StsBackTrace = -1,     //!< pseudo error for back trace
            StsError = -2,         //!< unknown /unspecified error
            StsInternal = -3,      //!< internal error (bad state)
            StsNoMem = -4,         //!< insufficient memory
            StsBadArg = -5,        //!< function arg/param is bad
            StsBadFunc = -6,       //!< unsupported function
            StsNoConv = -7,        //!< iteration didn't converge
            StsAutoTrace = -8,     //!< tracing
            HeaderIsNull = -9,     //!< image header is NULL
            BadImageSize = -10,    //!< image size is invalid
            BadOffset = -11,       //!< offset is invalid
            BadDataPtr = -12,      //!<
            BadStep = -13,         //!< image step is wrong, this may happen for a non-continuous matrix.
            BadModelOrChSeq = -14, //!<
            BadNumChannels =
                -15, //!< bad number of channels, for example, some functions accept only single channel matrices.
            BadNumChannel1U = -16,           //!<
            BadDepth = -17,                  //!< input image depth is not supported by the function
            BadAlphaChannel = -18,           //!<
            BadOrder = -19,                  //!< number of dimensions is out of range
            BadOrigin = -20,                 //!< incorrect input origin
            BadAlign = -21,                  //!< incorrect input align
            BadCallBack = -22,               //!<
            BadTileSize = -23,               //!<
            BadCOI = -24,                    //!< input COI is not supported
            BadROISize = -25,                //!< incorrect input roi
            MaskIsTiled = -26,               //!<
            StsNullPtr = -27,                //!< null pointer
            StsVecLengthErr = -28,           //!< incorrect vector length
            StsFilterStructContentErr = -29, //!< incorrect filter structure content
            StsKernelStructContentErr = -30, //!< incorrect transform kernel content
            StsFilterOffsetErr = -31,        //!< incorrect filter offset value
            StsBadSize = -201,               //!< the input/output structure size is incorrect
            StsDivByZero = -202,             //!< division by zero
            StsInplaceNotSupported = -203,   //!< in-place operation is not supported
            StsObjectNotFound = -204,        //!< request can't be completed
            StsUnmatchedFormats = -205,      //!< formats of input/output arrays differ
            StsBadFlag = -206,               //!< flag is wrong or not supported
            StsBadPoint = -207,              //!< bad CvPoint
            StsBadMask = -208,               //!< bad format of mask (neither 8uC1 nor 8sC1)
            StsUnmatchedSizes = -209,        //!< sizes of input/output structures do not match
            StsUnsupportedFormat = -210,     //!< the data format/type is not supported by the function
            StsOutOfRange = -211,            //!< some of parameters are out of range
            StsParseError = -212,            //!< invalid syntax/structure of the parsed file
            StsNotImplemented = -213,        //!< the requested function/feature is not implemented
            StsBadMemBlock = -214,           //!< an allocated block has been corrupted
            StsAssert = -215,                //!< assertion failed
            GpuNotSupported = -216,          //!< no CUDA support
            GpuApiCallError = -217,          //!< GPU API call error
            OpenGlNotSupported = -218,       //!< no OpenGL support
            OpenGlApiCallError = -219,       //!< OpenGL API call error
            OpenCLApiCallError = -220,       //!< OpenCL API call error
            OpenCLDoubleNotSupported = -221,
            OpenCLInitError = -222, //!< OpenCL initialization error
            OpenCLNoAMDBlasFft = -223
        };
    } // namespace Error

//////////////// static assert /////////////////

#define Subdiv2D_StaticAssert(condition, reason) static_assert((condition), reason " " #condition)

#if defined(__GNUC__)
#define SUBDIV2D_NORETURN __attribute__((__noreturn__))
#elif defined(_MSC_VER) && (_MSC_VER >= 1300)
#define SUBDIV2D_NORETURN __declspec(noreturn)
#else
#define SUBDIV2D_NORETURN /* nothing by default */
#endif

#ifndef SUBDIV2D_STATIC_ANALYSIS
#if defined(__KLOCWORK__) || defined(__clang_analyzer__) || defined(__COVERITY__)
#define SUBDIV2D_STATIC_ANALYSIS
#endif
#endif

    /*! @brief Signals an error and raises the exception.

    @param _code - error code (Error::Code)
    @param _err - error description
    @param _func - function name. Available only when the compiler supports getting it
    @param _file - source file name where the error has occurred
    @param _line - line number in the source file where the error has occurred
    @see Subdiv2D_Error, Subdiv2D_Error_, Subdiv2D_ErrorNoReturn, Subdiv2D_ErrorNoReturn_, Subdiv2D_Assert,
    Subdiv2D_DbgAssert
     */
    void error(int _code, const std::string& _err, const char* _func, const char* _file, int _line);

#ifdef __GNUC__
#if defined __clang__ || defined __APPLE__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winvalid-noreturn"
#endif
#endif

    /** same as cv::error, but does not return */
    static inline SUBDIV2D_NORETURN void errorNoReturn(int _code, const std::string& _err, const char* _func,
                                                       const char* _file, int _line) {
        error(_code, _err, _func, _file, _line);
#ifdef __GNUC__
#if !defined __clang__ && !defined __APPLE__
        // this suppresses this warning: "noreturn" function does return [enabled by default]
        __builtin_trap();
// or use infinite loop: for (;;) {}
#endif
#endif
    }
#ifdef __GNUC__
#if defined __clang__ || defined __APPLE__
#pragma GCC diagnostic pop
#endif
#endif

#if defined __GNUC__
#define Subdiv2D_Func __func__
#elif defined _MSC_VER
#define Subdiv2D_Func __FUNCTION__
#else
#define Subdiv2D_Func ""
#endif

#ifdef SUBDIV2D_STATIC_ANALYSIS
// In practice, some macro are not processed correctly (noreturn is not detected).
// We need to use simplified definition for them.
#define Subdiv2D_Error(...)                                                                                            \
    do {                                                                                                               \
        abort();                                                                                                       \
    } while (0)
#define Subdiv2D_Error_(...)                                                                                           \
    do {                                                                                                               \
        abort();                                                                                                       \
    } while (0)
#define Subdiv2D_Assert(cond)                                                                                          \
    do {                                                                                                               \
        if (!(cond))                                                                                                   \
            abort();                                                                                                   \
    } while (0)
#define Subdiv2D_ErrorNoReturn(...)                                                                                    \
    do {                                                                                                               \
        abort();                                                                                                       \
    } while (0)
#define Subdiv2D_ErrorNoReturn_(...)                                                                                   \
    do {                                                                                                               \
        abort();                                                                                                       \
    } while (0)

#else // SUBDIV2D_STATIC_ANALYSIS

    namespace detail {
        class MessageAccumulator {
          public:
            MessageAccumulator(const char* fmtString) { os_ << fmtString; }
            std::string get() const { return os_.str(); }

            template <typename T, typename... Args> void format(T&& arg, Args&&... a) {
                if (hasFirstArg_) {
                    os_ << ";";
                } else {
                    os_ << " [";
                    hasFirstArg_ = true;
                }
                os_ << (std::forward<T>(arg));
                format(std::forward<Args>(a)...);
            }
            /// base case
            void format() {
                if (hasFirstArg_) {
                    os_ << "]";
                }
            }

          private:
            bool hasFirstArg_ = false;
            std::ostringstream os_;
        };
        /// doesn't actually do printf-style formatting, but better than nothing.
        template <typename... Args> inline std::string format(const char* fmtString, Args&&... a) {
            MessageAccumulator msg(fmtString);
            msg.format(std::forward<Args>(a)...);
            return msg.get();
        }
    } // namespace detail

/** @brief Call the error handler.

Currently, the error handler prints the error code and the error message to the standard
error stream `stderr`. In the Debug configuration, it then provokes memory access violation, so that
the execution stack and all the parameters can be analyzed by the debugger. In the Release
configuration, the exception is thrown.

@param code one of Error::Code
@param msg error message
*/
#define Subdiv2D_Error(code, msg) ::sensics::subdiv2d::error(code, msg, Subdiv2D_Func, __FILE__, __LINE__)

/**  @brief Call the error handler.

This macro can be used to construct an error message on-fly to include some dynamic information,
for example:
@code
    // note the extra parentheses around the formatted text message
    Subdiv2D_Error_( ::sensics::subdiv2d::Error::StsOutOfRange,
    ("the value at (%d, %d)=%g is out of range", badPt.x, badPt.y, badValue));
@endcode
@param code one of Error::Code
@param args printf-like formatted error message in parentheses
*/
#define Subdiv2D_Error_(code, args)                                                                                    \
    ::sensics::subdiv2d::error(code, ::sensics::subdiv2d::detail::format args, Subdiv2D_Func, __FILE__, __LINE__)

/** @brief Checks a condition at runtime and throws exception if it fails

The macros Subdiv2D_Assert (and Subdiv2D_DbgAssert(expr)) evaluate the specified expression. If it is 0, the macros
raise an error (see cv::error). The macro Subdiv2D_Assert checks the condition in both Debug and Release
configurations while Subdiv2D_DbgAssert is only retained in the Debug configuration.
*/

#define SUBDIV2DAUX_CONCAT_EXP(a, b) a##b
#define SUBDIV2DAUX_CONCAT(a, b) SUBDIV2DAUX_CONCAT_EXP(a, b)

#define SUBDIV2D_VA_NUM_ARGS_HELPER(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define SUBDIV2D_VA_NUM_ARGS(...) SUBDIV2D_VA_NUM_ARGS_HELPER(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

#define Subdiv2D_Assert_1(expr)                                                                                        \
    if (!!(expr))                                                                                                      \
        ;                                                                                                              \
    else                                                                                                               \
        ::sensics::subdiv2d::error(::sensics::subdiv2d::Error::StsAssert, #expr, Subdiv2D_Func, __FILE__, __LINE__)
#define Subdiv2D_Assert_2(expr1, expr2)                                                                                \
    Subdiv2D_Assert_1(expr1);                                                                                          \
    Subdiv2D_Assert_1(expr2)
#define Subdiv2D_Assert_3(expr1, expr2, expr3)                                                                         \
    Subdiv2D_Assert_2(expr1, expr2);                                                                                   \
    Subdiv2D_Assert_1(expr3)
#define Subdiv2D_Assert_4(expr1, expr2, expr3, expr4)                                                                  \
    Subdiv2D_Assert_3(expr1, expr2, expr3);                                                                            \
    Subdiv2D_Assert_1(expr4)
#define Subdiv2D_Assert_5(expr1, expr2, expr3, expr4, expr5)                                                           \
    Subdiv2D_Assert_4(expr1, expr2, expr3, expr4);                                                                     \
    Subdiv2D_Assert_1(expr5)
#define Subdiv2D_Assert_6(expr1, expr2, expr3, expr4, expr5, expr6)                                                    \
    Subdiv2D_Assert_5(expr1, expr2, expr3, expr4, expr5);                                                              \
    Subdiv2D_Assert_1(expr6)
#define Subdiv2D_Assert_7(expr1, expr2, expr3, expr4, expr5, expr6, expr7)                                             \
    Subdiv2D_Assert_6(expr1, expr2, expr3, expr4, expr5, expr6);                                                       \
    Subdiv2D_Assert_1(expr7)
#define Subdiv2D_Assert_8(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8)                                      \
    Subdiv2D_Assert_7(expr1, expr2, expr3, expr4, expr5, expr6, expr7);                                                \
    Subdiv2D_Assert_1(expr8)
#define Subdiv2D_Assert_9(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9)                               \
    Subdiv2D_Assert_8(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8);                                         \
    Subdiv2D_Assert_1(expr9)
#define Subdiv2D_Assert_10(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9, expr10)                      \
    Subdiv2D_Assert_9(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9);                                  \
    Subdiv2D_Assert_1(expr10)

#define Subdiv2D_Assert(...) SUBDIV2DAUX_CONCAT(Subdiv2D_Assert_, SUBDIV2D_VA_NUM_ARGS(__VA_ARGS__))(__VA_ARGS__)

/** same as Subdiv2D_Error(code,msg), but does not return */
#define Subdiv2D_ErrorNoReturn(code, msg)                                                                              \
    ::sensics::subdiv2d::errorNoReturn(code, msg, Subdiv2D_Func, __FILE__, __LINE__)

/** same as Subdiv2D_Error_(code,args), but does not return */
#define Subdiv2D_ErrorNoReturn_(code, args)                                                                            \
    ::sensics::subdiv2d::errorNoReturn(code, cv::format args, Subdiv2D_Func, __FILE__, __LINE__)

#endif // SUBDIV2D_STATIC_ANALYSIS

/** replaced with Subdiv2D_Assert(expr) in Debug configuration */
#ifdef _DEBUG
#define Subdiv2D_DbgAssert(expr) Subdiv2D_Assert(expr)
#else
#define Subdiv2D_DbgAssert(expr)
#endif

} // namespace subdiv2d
} // namespace sensics

#endif // INCLUDED_AssertAndError_h_GUID_DCA9D718_0346_4088_A1AC_7EE0A243A538

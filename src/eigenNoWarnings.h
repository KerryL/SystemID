// File:  eigenNoWarnings.cpp
// Date:  9/26/2014
// Auth:  K. Loux
// Copy:  (c) Kerry Loux 2014
// Desc:  Wrapper for including Eigen headers that disables compiler warnings.
//        Allows using higher warning level for our own code without being
//        spammed with messages from 3rd party libraries.

// Disable warnings only for files included here
#if defined(__GCC__)
// Eigen headers included with -isystem option (not -I), so nothing required here
#elif defined(_MSC_VER)// MSVC++
#pragma warning(push)
#pragma warning(disable:4018)// signed/unsigned mismatch
#pragma warning(disable:4456)// declaration hides previous local declaration
#pragma warning(disable:4714)// function marked as __forceinline not inlined
#pragma warning(disable:4189)// local variable is initialized but not referenced
#pragma warning(disable:4127)// conditional expression is constant
#pragma warning(disable:4459)// declaration hides global declaration
#endif

// Eigen headers
#include <Eigen/Eigen>

// Remove our disables
#if defined(__GCC__)
// Nothing required here
#elif defined(_MSC_VER)// MSVC++
#pragma warning(pop)
#endif

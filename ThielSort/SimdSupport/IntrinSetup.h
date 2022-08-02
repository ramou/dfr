/* $Header:  $ */ 

#ifndef IntrinSetupH
#define IntrinSetupH

/********************************************************************
 *
 *                            IntrinSetup.H
 *
 * Defines, cache line len stuff, and some basic types.
 ********************************************************************/


#include <immintrin.h>

namespace SimdSupport {
  
typedef unsigned short UShort;
typedef unsigned int UInt;

#define SIMD_VEC(size) __attribute__ ((vector_size (size), may_alias))
#define SIMD_ALIGN(size) __attribute__ ((aligned(size)))
#define SIMD_ATTR(size) SIMD_VEC(size) SIMD_ALIGN(size)

#define FORCE_INLINE(RTN_VAL) inline RTN_VAL\
  __attribute__((__gnu_inline__, __always_inline__, __artificial__))
#define FORCE_INLINE_2 inline						\
  __attribute__((__gnu_inline__, __always_inline__, __artificial__))

// This seems as good a place as any for this.
#ifndef CACHE_LINE_LEN

#define CACHE_LINE_LEN 64
#define LOG_CACHE_LINE_LEN 6

static const int CacheLineLen = CACHE_LINE_LEN;
/**
 * Return true if addr is aligned on a CL.
 **/
inline bool IsClAligned(void *addr) 
{ return ((long)addr & (CACHE_LINE_LEN-1)) == 0; }

#endif  // ifndef CACHE_LINE_LEN

};  // namespace SimdSupport

#endif

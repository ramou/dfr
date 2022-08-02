/* $Header:  $ */ 

#ifndef IlElemH
#define IlElemH
#include <emmintrin.h>
#include "Prefetchers.h"

/********************************************************************
 *
 *                            IlElem.H
 *
 * Wrapper for dealing data type Elem_T into interleaved
 * ladles. Provides transparent operations for storing into the ladle
 * and copying out of it. If IL_N is 1, there is no interleaving
 * required. IL_N is a power of 2.
 *
 * Elem_T : type of data element
 * IL_N : degree of interleaving. If IL_N is 1, there is no interleaving
 *        required.
 ********************************************************************/


// IlElem                                                       IlElem
template<typename Elem_T, int IL_N>
class IlElem {
  Elem_T Val[IL_N];  // Only [0] is used. rest is undefined.
  static constexpr int IlNSHelper(int i) {
    if (i & 1) return 0;
    else return 1 + IlNSHelper(i >> 1);
  }
 public:
  static constexpr int IlN = IL_N;
  // IlN = 2^^IlNShift. Return IlNShift.
  static constexpr int IlNShift() { return IlNSHelper(IlN); }
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  IlElem() { Val[0] = 0;}
  IlElem(Elem_T v) { Val[0] = v; }
  //~IlElem();
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.
  IlElem & operator=(Elem_T v) { Val[0] = v; return *this; }
  // This should replace normal assignment. I hope.
  template<int IL_LEN> IlElem & operator=(IlElem<Elem_T, IL_LEN> &v) {
    Val[0] = (Elem_T)v;
    return *this;
  }
  operator Elem_T() const { return Val[0]; }
  Elem_T *GetEP() { return &Val[0]; }  // Compiles to noop (I hope) or lea??
};  // IlElem

// CopyElements                                                 CopyElements
// On titan, DEST_LOCALITY 0 is good, but SOURCE_LOCALITY -1 beats 3
// (slightly) Would doing it 1 ladle in advance help? This is with
// IL_N = 2. We should only prefetch the first time in a ladle
// set. This decision belongs closer to the logic!!
#define DEST_LOCALITY 0
#define SOURCE_LOCALITY -1
#define DO_PREFETCHING false
// Copy from interleaved ladle to a bucket with interleave == 1.
template<typename ELEM, int IL_N>
NO_INLINE void CopyElements(ELEM *dest_v, const IlElem<ELEM, IL_N> *src_v, long n) {
  //if (DEST_LOCALITY >= 0) PrefetchVector<1, DEST_LOCALITY>(dest_v, n);
  //if (SOURCE_LOCALITY >= 0) PrefetchVector<0, SOURCE_LOCALITY>(src_v, n);
  GdbHook(3);
  constexpr int UnwindN = 8;
  constexpr int NEPerCl = (CACHE_LINE_LEN + (sizeof(ELEM) - 1)) / sizeof(ELEM);
  constexpr int NILEPerCl =
    (CACHE_LINE_LEN + (sizeof(ELEM)* IL_N - 1)) / (sizeof(ELEM) * IL_N);
  constexpr int LookAhead = 2;
  if (DO_PREFETCHING) {
    // prefetch ahead 1 CL (dest only): +0.6% or so (sample size = 1).
    // But next run is much worse:-( . Third is better again, but not so much.
    for (int i = 0; i < LookAhead; i++) 
      for (int ii = 0; ii < UnwindN; ii += NEPerCl)
	__builtin_prefetch(dest_v + ii + i * UnwindN, 1, 0);
  }
  int i = 0;
  int loopEnd = n - (UnwindN - 1);
  for (; i < loopEnd; i += UnwindN) {
    if (DO_PREFETCHING) {
      for (int ii = 0; ii < UnwindN; ii += NEPerCl)
	__builtin_prefetch(dest_v + ii + LookAhead * UnwindN, 1, 0);
    }
#if 1    
    for (int ii = 0; ii < UnwindN; ii++) {
      dest_v[ii] = src_v[ii];
    }
#elif 0
    // Assuming UnwindN == 8, try streaming reads / writes
    // Almost identical with above.
    ELEM a = src_v[0]; ELEM b = src_v[1]; ELEM c = src_v[2]; ELEM d = src_v[3];
    dest_v[0] = a; dest_v[1] = b; dest_v[2] = c; dest_v[3] = d;
    a = src_v[4]; b = src_v[5]; c = src_v[6]; d = src_v[7];
    dest_v[4] = a; dest_v[5] = b; dest_v[6] = c; dest_v[7] = d;
#elif 0
    // This version looked good in windows on my desktop, but sucks
    // badly on titan! No idea why.
    long long int *p = reinterpret_cast<long long int*>(dest_v);
    // Assuming UnwindN == 8, try 64 bit NT writes
    ELEM a = src_v[0]; ELEM b = src_v[1]; ELEM c = src_v[2]; ELEM d = src_v[3];
    _mm_stream_si64(p + 0, a);
    _mm_stream_si64(p + 1, b);
    _mm_stream_si64(p + 2, c);
    _mm_stream_si64(p + 3, d);
    a = src_v[4]; b = src_v[5]; c = src_v[6]; d = src_v[7];
    _mm_stream_si64(p + 4, a);
    _mm_stream_si64(p + 5, b);
    _mm_stream_si64(p + 6, c);
    _mm_stream_si64(p + 7, d);
#else
    // This version 2% faster than above, but still sucks compared to
    // versions that don't use NT! Haven't tried the simd stuff.
    long long int *p = reinterpret_cast<long long int*>(dest_v);
    _mm_stream_si64(p + 0, src_v[0]);
    _mm_stream_si64(p + 1, src_v[1]);
    _mm_stream_si64(p + 2, src_v[2]);
    _mm_stream_si64(p + 3, src_v[3]);
    _mm_stream_si64(p + 4, src_v[4]);
    _mm_stream_si64(p + 5, src_v[5]);
    _mm_stream_si64(p + 6, src_v[6]);
    _mm_stream_si64(p + 7, src_v[7]);
#endif    
    dest_v += UnwindN;
    src_v += UnwindN;
  }
  // Now pick up leftovers
  n = n - i;
#if 0
  for (int ii = 0; ii < n; ii++) {
    dest_v[ii] = src_v[ii];
  }
#else
  // This version noticably (but not wonderfully) better: 0.5%.
  constexpr int bit0 = 1;
  if (bit0 > n) return;
  if (bit0 & n) {
    dest_v[0] = src_v[0]; dest_v += bit0; src_v += bit0;
  }
  constexpr int bit1 = bit0 << 1;
  if (bit1 > n) return;
  if (bit1 & n) {
    dest_v[0] = src_v[0];
    dest_v[1] = src_v[1];
    dest_v += bit1; src_v += bit1;
  }
  constexpr int bit2 = bit1 << 1;
  if (bit2 > n) return;
  if (bit2 & n) {
    dest_v[0] = src_v[0];
    dest_v[1] = src_v[1];
    dest_v[2] = src_v[2];    
    dest_v[3] = src_v[3];    
    dest_v += bit2; src_v += bit2;
  }
  constexpr int bit3 = bit2 << 1;
  if (bit3 > n) return;
  if (bit3 & n) {
    dest_v[0] = src_v[0];
    dest_v[1] = src_v[1];
    dest_v[2] = src_v[2];    
    dest_v[3] = src_v[3];    
    //dest_v += bit3; src_v += bit3;
  }
  static_assert(UnwindN <= bit3);

#endif

}
// Following versions all have the same code.

// Version to copy from a source with interleave == 1 to a dest
// with interleave == 1. Allows bulk move.
template<typename ELEM>
void CopyElements(ELEM *dest_v, const IlElem<ELEM, 1> *src_v, long n) {
  if (DEST_LOCALITY >= 0) PrefetchVector<1, DEST_LOCALITY>(dest_v, n);
  if (SOURCE_LOCALITY >= 0) PrefetchVector<0, SOURCE_LOCALITY>(src_v, n);
  std::memcpy(dest_v, src_v, n*sizeof(ELEM));
}
// Overloaded version to copy from IlElem<ELEM, 1> to ELEM to move
// ladle to bucket. Allows bulk move.
template<typename ELEM, int IL_N>
void CopyElements(ELEM *dest_v, const IlElem<ELEM, 1> *src_v,
                  long n) {
  if (DEST_LOCALITY >= 0) PrefetchVector<1, DEST_LOCALITY>(dest_v, n);
  if (SOURCE_LOCALITY >= 0) PrefetchVector<0, SOURCE_LOCALITY>(src_v, n);
  std::memcpy(dest_v, src_v, n*sizeof(ELEM));
}
// Overloaded version to copy from ELEM to Elem.
template<typename ELEM>
void CopyElements(ELEM *dest_v, const ELEM *src_v, long n) {
  if (DEST_LOCALITY >= 0) PrefetchVector<1, DEST_LOCALITY>(dest_v, n);
  if (SOURCE_LOCALITY >= 0) PrefetchVector<0, SOURCE_LOCALITY>(src_v, n);
  std::memcpy(dest_v, src_v, n*sizeof(ELEM));
}

#endif

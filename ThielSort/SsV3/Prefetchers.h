/* $Header:  $ */
#ifndef PREFETCHERS_H
#define PREFETCHERS_H

/**
 * Make unwound nested loop of prefetch requests to prefetch a
 * vector. We unwind 8 prefetches in the inner loop. 
 *
 * Template parameters:
 * R_W : 0 (read) or 1 (write) parameter to prefetch fn. Default : 0.
 *       I'm not sure if there is much difference!
 * Locality_Value : 0 (no locality), 1 (keep around for a short while),
 *                  2 (keep around for a longer while), 3 (KEEP).
 *                  Default : 0.
 * ELEM : Data type.
 *
 * Function parameters:
 * v_p : pointer to vector to be prefetched.
 * v_len : len of vector to prefetch (in ELEMs).
 *
 * CACHE_LINE_LEN must be visible. It need not be a constant, but why
 * wouldn't it be??
 *
 * Note: function call overhead is more instrs that expected. Should
 * lose the NO_INLINE for stable version.
 **/
template<int R_W = 0, int Locality_Value = 0, typename ELEM>
void PrefetchVector(const ELEM *v_p, int v_len) {
  constexpr int UnwindN = 1 << 3;  // 8. NOTE: change here implies change below
  constexpr int InnerLoopProtection = (UnwindN-1) * CACHE_LINE_LEN;
  const char *v = reinterpret_cast<const char*>(v_p);
  v_len = v_len * sizeof(ELEM);  // now in bytes
  const char *vEnd = v + (v_len - InnerLoopProtection);
  GdbHook(2);
  while (v < vEnd) {  
    for (int ii = 0; ii < UnwindN; ii++) {
      __builtin_prefetch(v + ii * CACHE_LINE_LEN, R_W, Locality_Value);
    }
    v += UnwindN * CACHE_LINE_LEN;
  }
  // Now prefetch any leftover CLs (n < UnwindN). This works as expected.
  int n = v_len & (UnwindN * CACHE_LINE_LEN - 1);
  constexpr int testBit0 = (UnwindN * CACHE_LINE_LEN) >> 1;
  if (testBit0 & n) {
    for (int ii = 0; ii < (UnwindN >> 1); ii++)
      __builtin_prefetch(v + ii * CACHE_LINE_LEN, R_W, Locality_Value);
    v += testBit0 * CACHE_LINE_LEN;
  }
  constexpr int testBit1 = testBit0 >> 1;
  if ((testBit1 >> 1) & n) {
    for (int ii = 0; ii < (UnwindN >> 2); ii++)
      __builtin_prefetch(v + ii * CACHE_LINE_LEN, R_W, Locality_Value);
    v += testBit1 * CACHE_LINE_LEN;
  }
  if ((testBit0 >> 2) & n) {
    for (int ii = 0; ii < (UnwindN >> 3); ii++)
      __builtin_prefetch(v + ii * CACHE_LINE_LEN, R_W, Locality_Value);
  }
}
/**
 * Similar to the above, but instead of prefetching a single vector,
 * we prefetch the next few cache lines for each vector in a vector of
 * vectors (i.e., pointers). This supports prefetching the next few
 * slots in a vector of buckets. Since the number of slots per vector
 * should be small, we don't worry about ELEMs split between cache
 * lines, and just step in terms of ELEMs.
 *
 * Template parameters R_W and Locality_Value are as above.
 * ELEM is now the ultimate pointed type, and v_p has type **ELEM.
 *
 * Function parameters:
 * v_p : pointer to vector of vectors to be prefetched.
 * v_len : len of vector (in buckets (*ELEM)).
 *
 * lead_n : we prefetch starting at v_p[i] + lead_n (in ELEMs).  if we
 *          are partway through a deal, the cache line containing
 *          *v_p[i] is probably already in cache, so we want to get
 *          the NEXT cache line, unless we are near the beginning of
 *          the current one, when requesting the current 1 again will
 *          avoid poluting the cache with something that won't be used
 *          for a while. Use lead_n to tune this behaviour.
 *
 * cnt_n : within each vector, prefetch cnt_n cache lines, not just 1.
 *         Default = 1. cnt_n<= 0 is a noop!
 **/
template<int R_W = 0, int Locality_Value = 0, typename ELEM>
void PrefetchBuckets( ELEM **v_p, int v_len, int lead_n, int cnt_n=1) {
  const int UnwindN = 16;
  const int NElemPerCL = CACHE_LINE_LEN / sizeof(ELEM);
  // In case v_p isn't a multiple of UnwindN :
  ELEM **vEnd = v_p + v_len - (UnwindN-1);

  for (; v_p < vEnd; v_p += UnwindN) {
    for (int ii = 0; ii < UnwindN; ii++) {
      const ELEM *p = v_p[ii];
      for (int iii = 0; iii < cnt_n; iii++) {
	__builtin_prefetch
	  (p + (lead_n + iii * NElemPerCL), R_W, Locality_Value);
      }
    }
  }
  vEnd += (UnwindN-1);
  for (; v_p < vEnd; v_p++) {
    const ELEM *p = *v_p;
    for (int iii = 0; iii < cnt_n; iii++) {
      __builtin_prefetch
	(p + (lead_n + iii * NElemPerCL), R_W, Locality_Value);
    }
  }

}


#endif

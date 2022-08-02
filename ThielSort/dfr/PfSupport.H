/* $Header:  $ */ 

#ifndef PfSupportH
#define PfSupportH

/********************************************************************
 *
 *                            PfSupport.H
 *
 * A few stand alone functions to support prefetching.
 ********************************************************************/

/**
 * Round p down to a CL addr (0 mod CacheLineLen) and return it as a char *. 
 * MakeCLAddr() caters to the case where neither the start nor the end
 * of the vec is CL aligned, so just skipping by CacheLineLen may
 * leave a partial CL at the end out. If start_p > end_p, a single CL
 * is prefetched if both addresses are in the same CL, else nothing is
 * prefetched.
 **/
char *MakeCLAddr(void *p) {
  if (sizeof(p) != sizeof(long long int)) MyAbort("Unexpected sizes!!");
  long long int pInt = reinterpret_cast<long long int>(p);
  pInt &= ~(CacheLineLen - 1);
  return reinterpret_cast<char *>(pInt);
}
const char *MakeCLAddr(const void *p) {
  if (sizeof(p) != sizeof(long long int)) MyAbort("Unexpected sizes!!");
  long long int pInt = reinterpret_cast<long long int>(p);
  pInt &= ~(CacheLineLen - 1);
  return reinterpret_cast<const char *>(pInt);
}

// PrefetchVec                                                  PrefetchVec
/**
 * Prefetch CLs for the range start_p : end_p (inclusive) using template
 * parameters R_W and Locality_Id. Locality_Id is :
 * 0 (don't keep in cache), 
 * 1 (keep for a little while), 
 * 2 (keep for a longer while), 
 * 3 (keep in all levels of cache a long time)
 *
 * R_W is : 0 (read) or 1 (write)
 * Pointers are const void* so caller needn't worry about type.
 *
 * Names are [Rev]PrefetchVec: 
 * if Rev is present, prefetch in reverse order, else forward order.
 *
 **/
template<int R_W, int Locality_Id>
  void PrefetchVec(const void *start_p, const void *end_p) {
  const char *lowClP = MakeCLAddr(start_p);
  const char *topClP =  reinterpret_cast<const char *>(end_p);  
  while (lowClP <= topClP) {
    __builtin_prefetch(lowClP, R_W, Locality_Id);
    lowClP += CacheLineLen;
  }
}  // PrefetchVec
template<int R_W, int Locality_Id>
  void RevPrefetchVec(const void *start_p, const void *end_p) {
  const char *lowClP = MakeCLAddr(start_p);
  const char *topClP =  reinterpret_cast<const char *>(end_p);
  while (lowClP <= topClP) {
    __builtin_prefetch(topClP, R_W, Locality_Id);
    topClP -= CacheLineLen;
  }
}  // RevPrefetchVec

#endif

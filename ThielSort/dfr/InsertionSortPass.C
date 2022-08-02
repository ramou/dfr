/* $Header:  $ */ 

#ifndef InsertionSortPassC
#define InsertionSortPassC

/********************************************************************
 *
 *                            InsertionSortPass.C
 *
 * Implimentation for InsertionSortPass<> class.
 ********************************************************************/

// DoBracketedSort_imp                                 DoBracketedSort_imp
NO_INLINE void DoBracketedSort_imp(Elem_T *src_p, Elem_T *src_e) {
  Elem_T *p = src_p;
  while (p < src_e) {
#if 0
    // gcc sees the reuse of p[1] as p[0] next iteration and prevents
    // the reload by moving p[1] to the p[0] temp and loading the new
    // p[1] into the p[1] temp. Not bad.
    while (p[0] <= p[1]) p++;
#else
    // For this one, gcc also sees the reuse of p[1] as p[0] and
    // prevents the reload. Here, it makes use of several registers so
    // it doesn't have to move, and uses indexed reads (p[1], p[2],
    // p[3],...) to save the inc. But that means that each break goes
    // to a separate (distant) location to fix up p. Also unwinds the
    // loop by 2 so we actually get 16 compares in the loop. Doubt
    // that pays!
    while (p[0] <= p[1]) {
      p++; if (p[0] > p[1]) break;
      p++; if (p[0] > p[1]) break;
      p++; if (p[0] > p[1]) break;
      p++; if (p[0] > p[1]) break;
      p++; if (p[0] > p[1]) break;
      p++; if (p[0] > p[1]) break;
      p++; if (p[0] > p[1]) break;
      p++; 
    }
#endif    
    // p[0] > p[1]. Thus p[0] not in sentinel, but p[1] may be.
    // Using ISNS_SHORT_RUN to test for short run.
    if (__builtin_expect(TopMask(p[ISNS_SHORT_RUN]) > TopMask(p), true)) {
      /**
       * Just do an IS. The worst case seems to be bad here because if
       * there are 1000 sorted elements in the SeBlock before we find
       * 1 out of order (the "sorted preface"), we may have to sift
       * each of Nf elements down 1000 times. If Nf is 10, the cost
       * will be 9000 compares (and moves). But a closer look belies
       * that this is a "high" cost: if we sort all 1009 elements
       * with an n*logn sort, and note that log2(1000) is 10, the cost
       * of the sort will be O(10000), but O() calcualtions ignore
       * various constant terms and factors, so the sift will
       * likely be cheaper. As the size of the sorted preface grows,
       * the cost of the sort increases in an nlogn fashion, while the
       * cost increase of the sift is strictly linear.
       *
       * One wants to do an average case analysis to get a good value
       * for Nf. My guess would be that numbers in the 5:10 range
       * would give good results. Remember that the worst case is not
       * just that we have a bunch of unsorted values at the end of a
       * large sorted preface, but that each of those unsorted numbers
       * belongs at the beginning of the SeBlock.
       **/ 
      //GdbHook(43);
      
      // Assuming only a single out of place saves 14% on IS at
      // N=50000000 compared to calling DoIs here!
      BC::LogCnts(0);
      SiftDown(p[1], p);
      p++;
    } else if (__builtin_expect(TopMask(p[Nf]) > TopMask(p), true)) {
      // short run
      DoIs(p, p + Nf);
      p = p + Nf - 1;  // Skipping to +Nf fails if p+Nf not an SeBlock bndry
    } else {
      if (__builtin_expect(TopMask(p + Nf) < TopMask(p), false)) {
	// p+Nf is past end. Do IS to end and done.
	int n = src_e - p;
	if (n > 1) {
	  DoIs(p, src_e);
	}  // Else p is just before sentinel, so nothing to do (or count).
	return;
      }
      // Found a longish SeBlock. Find its boundries and sort it.
      // +1 since FindNextSeb() returns first outside the seb.
      Elem_T *startP = BC:: template FindNextSeb<-1>(p, src_p-1) +1;
      // But the src_e is supposed to be outside the end, so good.
      Elem_T *endP = BC:: template FindNextSeb<1>(p, src_e);
      SortSeBlock(startP, endP, p);
      p = endP - 1;
    }
	  
  }
}
// CpySrt                                                  CpySrt
/**
 * Sort the semi-sorted non-empty chunk while copying the result to
 * NewDestV[]. There is a bottom sentinel in NewDestV.  The first
 * SeBlock is not an extension of a SeBlock started in the previous
 * chunk.
 **/
void CpySrt(Elem_T *src_p, long src_n) {
  if (TopMask(src_p) == TopMask(src_p[src_n - 1])) {
    // Single SeBlock.  NOTE: May not be large! Check if copy IS will work!
    BC::SortLargeSeBlock(src_p, src_p + src_n);
    std::memcpy(NewDestV, src_p, src_n * sizeof(Elem_T));
  } else if (src_n < IsThreshold) {
    DoCpyIs(src_p, src_p + src_n, NewDestV);
  } else {
    Elem_T *newEnd = BC::SortTopSentinelArea(src_p, src_p + src_n);
    // copy sorted tail to dest.
    long newN = newEnd - src_p;
    std::memcpy(NewDestV + newN, newEnd, (src_n - newN) * sizeof(Elem_T));
    // Create sentinel. 
    for (int i = 0; i < Nf; i++) newEnd[i] = BC::SmallestRec;
    BracketedCopySort(src_p, newEnd, NewDestV);
  }
  NewDestV += src_n;  // Update for next chunk.
}
/**
 * As CpySrt() but this is an overflow chunk so it may start with data
 * in the same SeBlock as the immediately previous chunk. If so,
 * combine the 2 pieces of that SeBlock efficiently into a sorted
 * whole, then call CpySrt() to deal with any further SeBlocks.
 **/
NO_INLINE void BucketExtend(Elem_T *src_p, long src_n) {
#if 1
  BC::IsNear(NewDestV, src_n);
#endif
  static const int SmallN = 4;
  CountType nToAdd = 0;  // Becomes size of partial SeBlock at start.
  Elem_T *endP = src_p;  // Becomes end of partial SeBlock at start.
  // Not a real loop. Just a scope to break out of.
  while (TopMask(src_p) == TopMask(NewDestV - 1)) {
    // Extending prior SeBlock.
    endP = BC:: template FindNextSeb<1>(src_p, src_p + src_n);
    nToAdd = endP - src_p;
    if (nToAdd <= SmallN) {
      DoCpyIs(src_p, endP, NewDestV);
      break;
    }
    if (nToAdd <= IsThreshold) {
      // DoIS() assumes src_p[1] out of order:-( . So find that condition
      Elem_T *p = src_p;
      for (; p < endP; p++) {
	if (p[1] < *p) break;
      }
      if (p != endP) {  // Found out of order condition.
	BC::NoSentIsSort(src_p, endP);
      }  // Else all recs already in order. 
    } else {  // Largish sort
      BC::SortLargeSeBlock(src_p, endP);
    }
    // Now merge the 2 sorted blocks into place.
    DoMerge(NewDestV + nToAdd - 1, NewDestV - 1, endP - 1);
    break;
  }
  NewDestV += nToAdd;
  if (src_n == nToAdd) return;
  CpySrt(endP, src_n - nToAdd);
}
/**
 * baseP and addonP point to the last elements of the primary and
 * overflow parts of an SeBlock, respectively. storeP > baseP is where
 * to store the biggest element from the 2 sorted list. When storeP
 * catches up (down?) with baseP, the merge is complete.
 **/
NO_INLINE void DoMerge(Elem_T *storeP, Elem_T *baseP, Elem_T *addonP) {
  Elem_T tAddon = *addonP;  // always == *addonP
  Elem_T tBase = *baseP;  // tBase is always == *baseP
  do {
    if (tBase <= tAddon) {
      *storeP-- = tAddon;
      if (storeP == baseP) break;
      tAddon = *--addonP;
    } else { *storeP-- = tBase; tBase = *--baseP; }
  } while (true);
}
/**
 * Do a modified insertion pass over src_p <= p < src_e.
 *
 * guaranttees at call:
 * - src_p < src_e.
 * - There is a bottom sentinel before the start of dest_v.
 *   NOTE: It may not be immediately before.
 * - There is a top sentinel (of size Nf) at src_e.
 **/
void BracketedCopySort(Elem_T *src_p, Elem_T *src_e, Elem_T *dest_v) {
  Elem_T *p = src_p;
  while (p < src_e) {
    // t = dest_v[-1] is top of sorted part. *dest_v slot to store next.
    // *p is next item to store into *dest_v unless it is < t.
    Elem_T t = dest_v[-1];
    while (t <= p[0]) {
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; if (t > p[0]) break;
      *dest_v++ = t = p[0];
      p++; 
    }
#ifdef IS_FINDING_INDEX
    BC::IsNear(dest_v, 20);
#endif
    // t == p[-1] == dest_v[-1] > p[0]. It may be top sentinel
    // if p[0] is first sentinel, p[Nf] is past sentinels!
    if (__builtin_expect(TopMask(p[Nf-1]) > TopMask(t), true)) {
      dest_v = DoCpyIs(p, p + Nf, dest_v);
      p = p + Nf;
    } else {
      if (__builtin_expect(TopMask(p[Nf-1]) < TopMask(t), false)) {
	// p+Nf is past end. Do IS to end and done.
	if (p < src_e) {
	  DoCpyIs(p, src_e, dest_v);
	  //SiftDown(*src_p++, dest_v++);
	  //if (src_p < src_e) DoCpyIs(src_p, src_e, dest_v);
	}  // Else *p == sentinel
	return;
      }
      // Found a longish SeBlock. Find its boundries and sort it.
      // startP is start of SeBlock.
      Elem_T *startP = BC:: template FindNextSeb<-1>(p, src_p - 1) +1;
      Elem_T *endP = BC:: template FindNextSeb<1>(p, src_e);
      long n = p - startP;
      if (n < IsThreshold) {
	dest_v = DoCpyIs(p, endP, dest_v);
      } else {
	BC::SortLargeSeBlock(startP, endP);
	dest_v -= n;
	while (startP < endP) *dest_v++ = *startP++;
      }
      p = endP;
    }
	  
  }
}

#endif

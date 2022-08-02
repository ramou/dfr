/* $Header:  $ */ 

#ifndef InsertionSortPassH
#define InsertionSortPassH

/********************************************************************
 *
 *                            InsertionSortPass.H
 *
 * A class of functions to complete the sort after an FR sort on some
 * high order portion of the key that may leave some unsorted
 * SeBlocks. See InsertionSortPass.Notes.
 *
 * NOTE: First version assumes Elem_T is both both the record and the
 * key.
 ********************************************************************/

#include <algorithm>

//#define IS_COUNTING

#include "IsConstBase.H"

// InsertionSortPass                                    InsertionSortPass
template<typename Elem_T>
class InsertionSortPass : private IsConstBase<Elem_T> {
  typedef IsConstBase<Elem_T> BC;

  using BC::Nf;
  using BC::TopMask;
  using BC::SiftDown;
  
  static const int IsThreshold = BC::IsThreshold;

  Elem_T *OriginalDest;  // Original start of dest. Needed to search down.
  Elem_T *NewDestV;  // After call to CopySort() where to store next chunk.
public:
  /**
   * Complete the sort by using an insertion sort to fix any
   * (supposedly) small out of order chunks after the diversion. If a
   * chunk is too large for that, we just use std::sort() or, in
   * principle, call DFR on that chunk.
   **/  
  void CompleteSort(Elem_T *src_p, Elem_T *src_e, Elem_T tm_mask) {
    src_p = BC::BaseInit(src_p, src_e, tm_mask, src_p);
    if (src_p != src_e) {
      // Bottom sentinel in place
      CountType n = src_e - src_p;
      if (n < IsThreshold) {  // Handle small case
	BC::NoSentIsSort(src_p, src_e);
	return;
      }
      Elem_T *newEnd = BC::SortTopSentinelArea(src_p, src_e);
      DoBracketedSort(src_p, newEnd);
    }  // Else all 1 SeBlock. Sorted in BaseInit()
  }
  /**
   * Similar to the above, but we need to copy the inorder values to
   * dest_v[] and the source is passed in 1 bucket (primary in src_p
   * (of size src_n), overflow in of_p (of size of_n)) at a
   * time. Either the primary source or the overflow (or both) may be
   * empty.
   *
   * Since the initial (non empty) call provides the address of the
   * destination vector and tm_mask, and must set up the bottom
   * sentinel in dest_v, there are 2 functions:
   *
   * 1 : CopySort() : skip any leading empty buckets and process the
   *     first non empty bucket after setting up the bottom
   *     sentinel. It returns true on a totally empty bucket, meaning
   *     it needs to be called again with the next bucket, or false if
   *     the initialization has been completed. The tm_mask and dest_v
   *     parameters should be the same for all calls.
   * 
   * 2 : ExtendCopySort() : Called for each bucket after CopySort()
   *     returns true.
   *
   * NOTE: The bucket sources will be changed. In particular, we
   * will install a sentinel near its end and we don't fix it
   * afterwards!
   *
   * Parameters (both functions):
   *
   * src_p and src_n delimit the bucket's primary source. If src_n ==
   * 0, it is empty.
   *
   * of_p and of_n delimit the bucket's overflow source. If of_n ==
   * 0, it is empty.
   *
   * If both sources are empty the functions do nothing.
   **/

  /**
   * Process all buckets up to and including the first non empty
   * bucket. Extra parameters :
   *
   * tm_mask is the top mask for the sort. It is 1 for the bits in the
   * key which were included in the radix sort and 0 elsewhere.
   *
   * dest_v is the destination where the sorted data will be
   * placed. It is assumed to be big enough.
   *
   * return true to request another bucket (since this one was empty)
   * or false to refuse further buckets, which should be passed to
   * ExtendCopySort().
   **/  
  bool CopySort(Elem_T *src_p, long src_n, Elem_T *of_p, long of_n,
		Elem_T tm_mask, Elem_T *dest_v) {
    if (src_n) {  // We have a primary source.
      BucketInit(src_p, src_n, tm_mask, dest_v);
      if (of_n) BucketExtend(of_p, of_n);
      return false;  // Don't call me again.
    } else if (of_n) {  // We only have an overflow sourc
      BucketInit(of_p, of_n, tm_mask, dest_v);
      return false;  // Don't call me again.
    } 
    return true;  // Else hit me again
  }
  /**
   * Add a new chunk to a dest previously passed to CopySort() while
   * sorting it. Context is in place.
   **/
  void ExtendCopySort(Elem_T *src_p, long src_n, Elem_T *of_p, long of_n) {
    if (src_n) {  // We have a primary source.
      CpySrt(src_p, src_n);
      if (of_n) BucketExtend(of_p, of_n);
    } else if (of_n) {  // We only have an overflow source
      CpySrt(of_p, of_n);
    }
  }
  
private:
  /**
   * Do a modified insertion pass over src_p <= p < src_e. 
   *
   * guaranttees at call:
   * - src_p < src_e.
   * - src_p[0] is smallest element in vector.
   * - src_e is addr of sorted tail of size >= Nf. This space can be used 
   *   as a sentinel.
   **/
  void DoBracketedSort(Elem_T *src_p, Elem_T *src_e) {
    // Save Nf values from end of vector here so we can replace them by
    // sentinels.
    Elem_T SentinelBackup[Nf];
    for (int i = 0; i < Nf; i++) {
      SentinelBackup[i] = src_e[i];
      src_e[i] = BC::SmallestRec;
    }

    DoBracketedSort_imp(src_p, src_e);

    // Now fix up sentinel space.
    for (int i = 0; i < Nf; i++) {
      src_e[i] = SentinelBackup[i];
    }
  }
  /**
   * Do a sort over the SeBlock src_p <= p < src_e. *ooo_p is the
   * first out-of-order rec (w.r.t. ooo_p[1]). We choose a sort
   * depending on the size of the block and the out-of_order
   * segment. We assume a bottom sentinel is in place.
   **/
  NO_INLINE void SortSeBlock(Elem_T *src_p, Elem_T *src_e, Elem_T *ooo_p) {
    GdbHook(10);
    int oooN = src_e - ooo_p;
    if (oooN < IsThreshold) {
      DoIs(ooo_p, src_e);
    } else {
      BC::SortLargeSeBlock(src_p, src_e);
    }
  }
  /**
   * Like above, but we create the 3rd parameter with the required
   * condition (op[0] > op[1]) before calling it.  Note that all is in
   * order up to src_e, we just return.
   **/
  void SortSeBlock(Elem_T *src_p, Elem_T *src_e) {
    Elem_T* op = src_p;
    Elem_T* oEnd = src_e-1;
    while ((op < oEnd) && (op[0] <= op[1])) op++;
    if (op < src_e) SortSeBlock(src_p, src_e, op);
  }
  /**
   * Do an insertion sort in place.
   * Assumes src_p[0] is in place vis a vis earlier records in src_p.
   * Assumes src_p[0] > src_p[1]. So src_e > src_p+1. 
   * Assumes top sentinel doesn't start before src_e.
   * Assumes bottom sentinel is in place.
   **/
  void DoIs(Elem_T *src_p, Elem_T *src_e) {
    BC::LogCnts(src_e - src_p);
    src_e--;  // Since we look at src_p[1]!
    // Set up local sentinel.
    Elem_T tmp = src_e[1];
    src_e[1] = BC::SmallestRec;  // local top sentinel    
    do {
      SiftDown(src_p[1], src_p);
      while (src_p[0] <= src_p[1]) src_p++;
    } while (src_p < src_e);
    src_e[1] = tmp;
  }
  /**
   * Do a copy insertion sort.
   * Assumes bottom sentinel is in place (may sift down below src_p).
   * Assumes src_p < src_e (n > 0). Store the first copied item at dest_v,
   * and return the addr of the first free spot in dest_v.
   **/
  Elem_T *DoCpyIs(Elem_T *src_p, Elem_T *src_e, Elem_T *dest_v) {
    BC::LogCnts(src_e - src_p);
    do {
      if (src_p[0] < dest_v[-1]) {
	SiftDown(src_p[0], dest_v - 1);
      } else *dest_v = src_p[0];
      dest_v++;
    } while (++src_p < src_e);
    return dest_v;
  }

  /**
   * src_p, src_n are the first non empty chunk. Init, set up bottom
   * sentinel, sort the chunk, and move it to dest.
   **/  
  //////////////  Made changes here for new BaseInit(). Not tested!!
  void BucketInit(Elem_T *src_p, long src_n, Elem_T tm_mask, Elem_T *dest_v) {
    Elem_T *srcE = src_p + src_n;
    if (src_n < 12) {  // My guess at reasonable. Talk to theorists!!

      this->NoSentIsSort(src_p, srcE);
      std::memcpy(dest_v, src_p, src_n * sizeof(Elem_T));
      NewDestV = dest_v + src_n;
    } else {
      Elem_T newSrcP = BC::BaseInit(src_p, srcE, tm_mask, dest_v);
      long sebN = srcE - src_p;
      std::memcpy(dest_v, src_p, sebN * sizeof(Elem_T));
      NewDestV = dest_v + sebN;
      if (newSrcP != srcE) {  // Still work to do.
	CpySrt(newSrcP, src_n - sebN);
      }  // else sorted in BaseInit().
    }
  }
#include "InsertionSortPass.C"
  
};  // InsertionSortPass

#endif

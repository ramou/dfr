/* $Header:  $ */ 

#ifndef SimdMoveH
#define SimdMoveH

/********************************************************************
 *
 *                            SimdMove.H
 *
 * Holds the functions that move simd units (and single Uts). It extends
 * SrtSimdTypes<>, where the actual intrinsics are called. No data
 * except constants to describe the simd units and associated
 * alignments and block sizes (see below). All functions are static.
 * All address parameters to public functions are Rec_Type*.
 *
 * Doesn't depend on or call anything except from Simd_Type.
 ********************************************************************/

//#include "SrtType.H"
#include "SimdUnitMover.h"

namespace SimdUnitMover {

  // SimdMove                                                       SimdMove
  /**
   * Template parameter Simd_Log_Size is log of size of simd unit in bytes.
   **/
  template<class Rec_Type, int Simd_Log_Size = 5> struct SimdMove {
    //typedef Rec_Type SortRecType;
    
    typedef SimdUnitMover::RemoteSdu<Simd_Log_Size> RemoteSdu;
    typedef SimdUnitMover::LocalSdu<Simd_Log_Size> LocalSdu;

    static const int SimdLogSize = Simd_Log_Size;
    static const int ByteSize = 1 << Simd_Log_Size;
    static const int UtsPerSimd =  ByteSize / sizeof(Rec_Type);
    static const int NUtsPerCL = CACHE_LINE_LEN / sizeof(Rec_Type);
  
    // Describe a block of simd units that we try to move efficiently.
    // Number of simd units in a simd block (an unwound loop)
    static const int LogSimdsPerBlock = 2;
    static const int SimdsPerBlock = 1 << LogSimdsPerBlock;
    // Number of items in a simd block.
    static const int UtsPerBlock = UtsPerSimd * SimdsPerBlock;

    /**
     * Functions to round n down and round up to a multiple of
     * mul, which must be a power of 2. 
     **/
    static constexpr int RoundDown(int n, int mul)
    { return n & ~(mul - 1); }
    static constexpr int RoundDownSimd(int n)
    { return RoundDown(n, UtsPerSimd); }
    static constexpr int RoundDownBlk(int n)
    { return RoundDown(n, UtsPerBlock); }
    // Return n mod mul, where mul defaults to UtsPerSimd and must be 2**x.
    static constexpr int ModSimd(int n, int mul = UtsPerSimd)
    { return n & (mul - 1); }
    // Return n mod UtsPerBlock.
    static constexpr int ModSimdBlk(int n) { return ModSimd(n, UtsPerBlock); }
    // Return number of sdus in n_recs recs. May be 0.
    static constexpr int GetNSimds(int n_recs) {
      return n_recs / UtsPerSimd;
    }
    // Return number of blocks in n_recs recs. May be 0.
    static constexpr int GetNBlks(int n_recs) {
      return n_recs / UtsPerBlock;
    }
    // Return true iff p is aligned on a Simd boundry.
    static bool IsAligned(const void *p) {
      ULong Mask = (UtsPerSimd * sizeof(Rec_Type)) - 1;
      return ((ULong)p & Mask) == 0;
    }
    /**
     * Prefetching the source CLs before copying them is probably always
     * good. Doing the copy in reverse order is also probably good,
     * since ladles are filled from the bottom. Those filled later are
     * more likely still in cache. Select the loop order and the
     * prefetch option using this enum. 
     **/
    enum SimdMoveOptions {
      ReverseLoop = 1,  
      PrefetchOption = 2
    };
    static const int DefaultMoveOption = ReverseLoop + PrefetchOption;
    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor : None

    //////////////////////////////////////////////////////////////////////
    //                   Other functions.

    // Primary functions. No logic, just move homogeneous
    // stuff. 
  
    /**
     * Prefetch ladle in prep for copying it to the bucket.
     * src_v - addr of start of area to copy.
     * end_v - addr just past area to copy. 
     * move_option
     **/
    INLINE_ATT static
    void DoPrefetch(const Rec_Type *src_v, const Rec_Type *end_v,
		    int move_option = DefaultMoveOption) {
      if (move_option & PrefetchOption) {  
	if (move_option & ReverseLoop) {
	  for (end_v -= NUtsPerCL; end_v >= src_v; end_v -= NUtsPerCL) {
	    __builtin_prefetch(end_v, 0, 3);
	  }
	} else {
	  for ( ; end_v >= src_v; src_v += NUtsPerCL) {
	    __builtin_prefetch(src_v, 0, 3);
	  }
	}
      }
    }
    /**
     * Use reverse move to copy contents of src_v : src_end to the space
     * immediately before dest_end. Uses ReverseLoop and PrefetchOption.
     * Assumes at least 1 block must be moved. src_v : src_end is some
     * number of complete blocks.
     **/
    //INLINE_ATT static
    NO_INLINE static
    void RevBulkMove(Rec_Type *dest_end, const Rec_Type *src_v,
		     const Rec_Type *src_end) {
      GdbHook(2);
      //DoPrefetch(src_v, src_end, ReverseLoop | PrefetchOption);
      LocalSdu t1, t2, t3, t4;
      RemoteSdu *destEP = RemoteSdu::SimdPtr(dest_end);
      const LocalSdu *srcV = LocalSdu::SimdPtr(src_v);
      const LocalSdu *srcEP = LocalSdu::SimdPtr(src_end);
      const int NWay = 4;
      t1 = srcEP[-1];
      t2 = srcEP[-2];
      t3 = srcEP[-3];
      t4 = srcEP[-4];
      srcEP -= NWay;
      while (srcEP > srcV) {
	destEP[-1] = t1;
	destEP[-2] = t2;
	destEP[-3] = t3;
	destEP[-4] = t4;
	destEP -= NWay;
	t1 = srcEP[-1];
	t2 = srcEP[-2];
	t3 = srcEP[-3];
	t4 = srcEP[-4];
	srcEP -= NWay;
      } 
      destEP[-1] = t1;
      destEP[-2] = t2;
      destEP[-3] = t3;
      destEP[-4] = t4;
    }
    /**
     * Move n_uts (>=0) records from src_v to dest_v in units of
     * UtsPerSimd. n_uts must be a multiple of UtsPerSimd. Not super
     * efficient for large n_uts, but does work.
     *
     * Note: May break CL alignment (unless enforced by caller) but
     * preserves simd alignment.
     **/
    //INLINE_ATT static
    NO_INLINE static
    void MoveSimds(Rec_Type *dest_v, const Rec_Type *src_v, int n_uts) {
      GdbHook(3);
      const int MaxToUnwind = 3;  // Will want more later!!
      if (n_uts == 0) return;
      const LocalSdu *src = LocalSdu::SimdPtr(src_v);
      RemoteSdu *dest = RemoteSdu::SimdPtr(dest_v);
      int nSimds = n_uts / UtsPerSimd;
      if (nSimds > MaxToUnwind) {
	for (int i = 0; i < nSimds; i++) dest[i] = src[i];
      } else if (nSimds == 3) {
	LocalSdu t0 = src[0];
	LocalSdu t1 = src[1];
	LocalSdu t2 = src[2];
	dest[0] = t0;
	dest[1] = t1;
	dest[2] = t2;
      } else if (nSimds == 2) {
	LocalSdu t0 = src[0];
	LocalSdu t1 = src[1];
	dest[0] = t0;
	dest[1] = t1;
      } else if (nSimds == 1) {
	LocalSdu t0 = src[0];
	dest[0] = t0;
      } else MyAbort("MaxToUnwind mismatch in MoveSimds");
    }
    /**
     * Move n_uts (>=0) Uts from src_v to dest_v. Expected n_uts is
     * smaller than a simd. May break simd alignment at dest_v. 
     **/  
    INLINE_ATT static void MoveUts(Rec_Type *dest_v, const Rec_Type *src_v,
				   int n_uts) {
      if (n_uts == 0) return;
      if (n_uts == 1) { dest_v[0] = src_v[0]; return; }
      for (int i = 0; i < n_uts; i++) dest_v[i] = src_v[i];
      // if (n_uts == 3) use simd move??
    }
    // Secondary functions                                 Secondary functions

    // Secondary functions reflect higher level, aggregate units,
    // which they break up and send in pieces to the primary functions.
    /**
     * Move n_uts records from src_v to dest_v in units of UtsPerBlock.
     * n_uts must be a multiple of UtsPerBlock.
     *
     * NOTE: Assumption that when we need to flush, at least 1 block
     * must be moved adds complications when we must split the flush to
     * different locations as 1 buffer gets filled. In this case, 1 or
     * both of parts may or may not satisfy that assumption. So I am
     * dropping that idea (for now).
     **/
    INLINE_ATT static
    void RevMoveBlock(Rec_Type *dest_v, const Rec_Type *src_v,
		      int blk_n, int simd_n) {
      MoveSimds(dest_v + blk_n, src_v + blk_n, simd_n);
      RevBulkMove(dest_v + blk_n, src_v, src_v + blk_n);
    }

    /**
     * The paired sizes above allow optimization when simd_n is
     * literal. If only the total length is known, it is just a pain to
     * call. This version takes the total n. It still assumes n is an
     * integral number of simd units but doesn't require >= a simd
     * block. Thus it is usable for cases where an originally large
     * enough chunk of data has to be split into smaller chunks and
     * stored in multiple places.
     **/
    INLINE_ATT static void MoveBlock
    (Rec_Type *dest_v, const Rec_Type *src_v, int n) {
      if (n >= UtsPerBlock) {
	RevMoveBlock(dest_v, src_v, RoundDownBlk(n), ModSimdBlk(n));
      } else {
	MoveSimds(dest_v, src_v, n);
      }
    }
  private:

  };  // SimdMove

};  // SimdUnitMover
#endif

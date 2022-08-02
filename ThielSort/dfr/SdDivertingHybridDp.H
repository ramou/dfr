/* $Header:  $ */

#ifndef SdSdDivertingHybridDpH
#define SdSdDivertingHybridDpH

/**
 *                     SdDivertingHybridDp.H 
 *
 * This Sort descriptor file declares a templated sort descriptor
 * class which defines constants and types for a set of related sort
 * versions. A sort is done by calling DoSort(). Variations within the
 * set are specified with template parameter:
 *
 * SRec_T : the type of the sort record.
 *
 **/

#include "BcSortDesc.H"

template<typename SRec_T>
class SdClass<SRec_T, FR_DIVERTING_WITH_HYBRID_FLUSHER_SD> :
  protected BcSortDesc {

  typedef SRec_T ElemT;
  typedef IlElem<SRec_T, IlN> IlElemT;
  typedef SRec_T BucketT;

  using DealerCore = DealerCoreTemplate<SRec_T>;
  using FlushMoverT = FlushMovers::FMoverT<SRec_T>;
  using FlusherT = FlusherTemplate<SRec_T>;
  using AlignerT = typename FlushMoverT::AlignerT;

  static constexpr int RtPerMu = FlushMoverT::RtPerMu;
  
  static void ExitDoSortMsg(const char *msg = "") {
#if 0
    std::cout << "Exit DoSort() at RecursiveLevel = " << RecursiveLevel <<
      ' ' << msg << std::endl;
#endif
  }

public:
  
  /**
   * This file implements the following sort:
   **/
  static constexpr const char *SdSpec =
    "Diverting direct ladle. Fast radix,  n < 8 passes.(SdDivertingHybridDp.H)";

  // Define the bucket sizes using 8 bit digits
  typedef BFDefinition<0, 8, 8, 8, 8, 8, 8, 8, 8> BFDefT;
  template<int Fld_Id>
  using BitFld = BitFieldAccess<Fld_Id, BFDefT>;

  static constexpr int NFields = BFDefT::NFields;

  /////////////////////////////////////////////////////////////////////////
  // Bug (fixed) when a lot of identical values showed with (flatter) rng spec:
  // -r  Normal 4253852969 9.223372E+18 2.305843E+18  1000000000 
  // Pointier test has rng spec:
  // -r Normal 1 4.3E9 1.07E9 1000000000
  // -r Normal 368392917 4.3E9 1.07E9 1000000000 loops somewhere
  /////////////////////////////////////////////////////////////////////////
  void DoSort(ElemT *src_v, CountType src_n,
	      int n_effective_flds = NFields, ElemT top_mask = 0) {
#if 0
    std::cout << "In DoSort() at RecursiveLevel = "<< RecursiveLevel <<
      std::endl;
     std::cout << "Sorting vector at " << src_v << " with length = " <<
      src_n << '\n';      
#endif
    // See if recursion through insertion sort has shown that the
    // whole block is identical records. Return if so.
    if (n_effective_flds == 0) {
      ExitDoSortMsg("n_effective_flds == 0");
      return;
    } else {  // This is a recursive sort.
    }
    // Use these values to define space available for estimated buckets.
    // src_v and src_n still used in most places.
    ElemT *AlignedBktSpace = AlignerT::AlignMuTUp(src_v);
    CountType AlignedBktSpaceN = src_n - (src_v - AlignedBktSpace);
	
    DataPoolPair<SRec_T, BFDefT> PoolPair(SdSpec, src_v, src_n);
    // Suppress dumpting config info when recursing.
    if (DumpConfig && (n_effective_flds==NFields)) PrintDiversionThresholds();
    
    // If counting a pass, the counts are stored here.
    CountType BktCntr[MaxNBuckets+2];  
    // If PreferNewBuffer, don't use original source as a bkt bfr.
    // With variable number of passes, must be careful with buffers!!
    ElemT *BktBfr[2] = {nullptr, AlignedBktSpace};
    // Select the BktBfr to use for this pass. 0 selects temp space, 1
    // selects the original source vector (usually).
    int BktBfrNdx = 0;
    if (PreferNewBuffer) {
      // PreferNewBuffer for this Sd seems bad on titan. But good when
      // bypassing ladles?? Don't understand:-((
      // Don't use original src as bkt buffer on even numbered passes?
      BktBfr[1] = nullptr;
    }

    LbDataInfoT<ElemT> LbDataInfo(src_v[0]);
    NoopDataInfoT<ElemT> NoopDataInfo;
    CountDataInfoT<ElemT> CountDataInfo(0, BktCntr);
    // If looking for high order unused bits with sampling, and it
    // nPasses is 2, will need to use CntLbDataInfo on first pass.
    CntLbDataInfoT<ElemT> CntLbDataInfo(src_v[0], BktCntr);
    BktCntr[MaxNBuckets] = src_n;
    
    FlusherT Flusher(src_n, &PoolPair.Timers);
    PoolPair.DppInit(&Flusher, BatchSize);
    
    int NBuckets = MaxNBuckets;  // NBuckets const in this SD.
    
    int BucketLens[MaxNBuckets + 1];  // ~[NBuckets] becomes total len
    GuessBkts(BucketLens, AlignedBktSpaceN, RtPerMu, NBuckets);
    int nPasses = FindNDeals(src_n);
    ElemT topMask = top_mask;
    // Select digit for first pass based on total number of passes
    // desired.
    
    // NOTE: Fix if dropping unused high order digit: use
    //       n_effective_flds as number of digits, not NFields.
    int lowOrderDigit = n_effective_flds + 1 - nPasses;
    // Check for < nPasses digits left
    while (lowOrderDigit < 1) {
      lowOrderDigit++;
      nPasses--;
    }
    //std::cout<<"lowOrderDigit="<<lowOrderDigit<<", nPasses= "<<nPasses<<'\n';
    PoolPair.StopSetupTimer();
    // First pass needs Lb check if dropping high order bits on sampling.
    // I left first pass separate so this can be added when needed.
    // Loop ends need fixing if dropping high order bits on sampling.
    if ((nPasses & 1) == 1) {  // Odd number of passes. Needs copy at end.
      if (nPasses == 1) {  // 1 pass not counted, then copy.
        topMask |= PoolPair.template DealAPass
          (lowOrderDigit,nullptr, BucketLens, &NoopDataInfo);
      } else {  
	// Since we will need a copy at the end, no need to count
	// for last pass. Copy will move data back home.
        topMask |= PoolPair.DealAPass
          (lowOrderDigit,nullptr, BucketLens, &NoopDataInfo);
        BktBfrNdx ^= 1;  // Didn't use this above, but swap anyway.
	for (int curDigit = lowOrderDigit+1; curDigit<=n_effective_flds; curDigit++) {
	  topMask |= PoolPair.DealAPass
	    (curDigit, BktBfr[BktBfrNdx], BucketLens, &NoopDataInfo);
	  BktBfrNdx ^= 1;
	}
      }
      // Now copy back to original source and do insertion sort.
#if USE_COPY_INSERTION_SORT == 1
      PoolPair.RunCopyInsertionSort(src_v, topMask);
#else      
      PoolPair.CopyDpToVec(src_v);
      MergedOfMem.MarkFree();  // Release OF area.
      PoolPair.RunInsertionSort(topMask, src_v, src_v + src_n);
#endif
    } else {  // Even number of passes. Count during penultimate pass,
      // and last pass stores data in place.
      // I think still counting because need to allocate and create
      // a second temp vector puts overhead through the roof!!
      if (nPasses == 2) {
	topMask |= PoolPair.template DealAPass
	  (lowOrderDigit,nullptr, BucketLens, &CountDataInfo);
	Flusher.SetCountedPass(true);
	topMask |= PoolPair.template DealAPass
          (lowOrderDigit + 1, src_v,BktCntr,&NoopDataInfo);
      } else {  // > 2 passes needed.
	topMask |= PoolPair.DealAPass
	  (lowOrderDigit,nullptr, BucketLens, &NoopDataInfo);
	BktBfrNdx ^= 1;  // Didn't use this above, but swap anyway.
	int curDigit = lowOrderDigit + 1;
	for (; curDigit < (n_effective_flds - 1); curDigit++) {
	  topMask |= PoolPair.DealAPass
	    (curDigit, BktBfr[BktBfrNdx], BucketLens, &NoopDataInfo);
	  BktBfrNdx ^= 1;
	}
	// Collect counts last pass so data will end up back home.
	topMask |= PoolPair.DealAPass
	  (curDigit, BktBfr[BktBfrNdx], BucketLens, &CountDataInfo);
	// Use exact counts (BktCntr) to put data in final home.
	// NOTE: Can't use BktBfr[BktBfrNdx] for last pass if
	// PreferNewBuffer!!
	Flusher.SetCountedPass(true);	
	topMask |= PoolPair.DealAPass
	  (curDigit + 1, src_v, BktCntr, &NoopDataInfo);
      }  // nPasses == 2
      // Shouldn't be needed after counted pass, but it is. Why??
      MergedOfMem.MarkFree();  // Release OF area.
      if (nPasses < n_effective_flds) {  // Still unsorted digits
	PoolPair.RunInsertionSort(topMask, src_v, src_v + src_n);
      }  // else radix sort includes loworder (and thus all) digits.
    }

    PoolPair.PrintInfo();
    ExitDoSortMsg();

  }  // DoSort

};  // SdDivertingHybridDp

#endif


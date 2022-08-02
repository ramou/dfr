/* $Header:  $ */ 

#ifndef DataPoolPairH
#define DataPoolPairH

#include "TimeWrapper.H"

/********************************************************************
 *
 *                            DataPoolPair.H
 *
 * The DataPoolPair class contains 2 DataPools, the objects that they
 * share, and the functions that configure and deal passes. It swaps
 * the source and target each pass.
 *********************************************************************/

// DataPoolPair                                                DataPoolPair
template<typename SRec_T, typename BFld_Defn_T> class DataPoolPair {

  typedef SRec_T ElemT;
  using DataPoolType = DataPoolTemplate<SRec_T>;

  using FlushMoverT = FlushMovers::FMoverT<SRec_T>;
  using FlusherT = FlusherTemplate<SRec_T>;

  template<int Fld_Id>
  using BitFld = BitFieldAccess<Fld_Id, BFld_Defn_T>;

  static constexpr int NFields = BFld_Defn_T::NFields;
  
  void Swap() {
    SrcDpNdx = 1 - SrcDpNdx;
    SrcDp = DataDefn + SrcDpNdx;
    TargetDp = DataDefn + (1 - SrcDpNdx);
  }
  // Called by both constructors to do what is common to them
  void InitConstruct(const char *sd_name, long n_to_sort) {
    if (RecursiveLevel == 0) {
      PrintConfig(n_to_sort, FlusherT::WhoIsLM(), FlusherT::WhoAmI(),
		  DealerCoreTemplate<SRec_T>::WhoAmI(),
		  FlushMoverT::WhoAmI(), sd_name);
    } else {
      //std::cout << "Recursing to level " << RecursiveLevel << std::endl;
    }
    SetupTimerId = Timers.TimeWInit(n_to_sort);  // Start timing here.
    for (int i = 0; i <= NFields; i++) TimerId[i] = -1;
  }
 public:

  // Public so Sd classes can use it directly.
  TimeWrapper Timers;
  
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  // If sorting data by sections, pass the n_to_sort for the entire
  // job. Use DppInit() to initialize the various sections.
  //
  // NOTE: When using segments, each final result must be extracted
  // from it's DataPool before starting the next segment.
  DataPoolPair(const char *sd_name, long n_to_sort) {
    InitConstruct(sd_name, n_to_sort);
  }
  // Older call when not using segments.
  // Use DppInit(FlusherT *flusher_p, int batch_size).
  DataPoolPair(const char *sd_name, ElemT *src_vec, long n_to_sort) {
    InitConstruct(sd_name, n_to_sort);
    //std::cout << "DataPoolPair: " << src_vec << ", n_to_sort= " << n_to_sort << '\n';
    SourceV = src_vec;
    SortN = n_to_sort;
  }
  // Init the sort for a segment of the total job.

  // Init DataPoolPair when entire source vector was passed to constructor.
  void DppInit(FlusherT *flusher_p, int batch_size) {
    SortSegment++;
    /**
     * Looks like smaller ladles don't help after about 100000. No
     * surprise there. We see gain at 7500. (my I5 desktop).
     **/    
    static const int MAX_LADLE_THRESHOLD = 100000;  
    if (!USE_EFFECTIVE_BATCH_SIZE || (SortN >= MAX_LADLE_THRESHOLD)) {
      BatchSize = batch_size;
    } else if (SortN < batch_size)  {  // We can shorten the ladle
      BatchSize = SortN;
    } else {
      BatchSize = batch_size / 2;  // Can be tuned more precisely, I'm sure!!
    }
    FlusherP = flusher_p;
    SrcDpNdx = 0;
    SrcDp = DataDefn + SrcDpNdx;
    SrcDp->WrapSourceVec(SourceV, SortN);
    TargetDp = DataDefn + (1 - SrcDpNdx);
  }
  // Init DataPoolPair for each segment when using
  // DataPoolPair(const char *sd_name, long n_to_sort).
  // src_vec and n_to_sort describe the current (or only) segment to sort.
  void DppInit(FlusherT *flusher_p, int batch_size, ElemT *src_vec,
	       long n_to_sort) {
    SourceV = src_vec;
    SortN = n_to_sort;
    DppInit(flusher_p, batch_size);
  }
  /**
   * Specify the expected temp space vector size needed to be
   * allocated when using a segmented data vector. Should be done
   * before first call to DppInit() to avoid repeatedly reallocating
   * for a slightly larger size.
   **/
  void SetTargetBfrSize(CountType size) {
    DataDefn[1].InitBuffer(size);
    BucketBuffer.MarkFree();    
  }
  //~DataPoolPair();
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  void StopSetupTimer() { Timers.StopSetup(); }
  void ReStartSetupTimer() { Timers.ReStartTimer(SetupTimerId); }
  void PrintInfo() { Timers.PrintInfo(); }  
  /**
   * These functions redirect to the appropriate dp, and swap source and
   * target dp after a deal.
   **/
  
  /**
   * Functions to deal a pass from a DataPool. Parameters:
   *
   * - ElemT *bkt_buffer : Buffer to deal to. If NULL, use locally
   *   allocated temp buffer. Its length >= bkt_lens[numBuckets],
   *   where numBuckets comes from DBA.
   *
   * - int *bkt_lens : bkt_lens[b] contains the exact length to be
   *   allocated to bucket[b] before overflow. May come from a real
   *   exact count or an estimate.
   *
   * ** NOTE: bkt_lens[numBuckets] = sum(bkt_lens[0:n_buckets-1]) ==
   *          size of DataPool to deal to.
   *
   * - data_info_p : Carries, Data_Info type and allows transfer of
   *   DataInfo values between copy visible to caller and copy used by
   *   Dealer.
   *
   **/

  /**
   * Select a digit based on a simple int parameter and call
   * DealFromDP() based on that. Means that the Sd versions don't
   * have to worry about matching the pass with the digit.
   **/
  template<typename Data_Info> NO_INLINE
  ElemT DealAPass(int digit_id, ElemT *bkt_buffer, int *bkt_lens,
		 Data_Info *data_info_p) {

    switch (digit_id) {
    case 1: {
      DealFromDP<BitFld<1> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<1>::InPlaceMask();
    }
    case 2: {
      DealFromDP<BitFld<2> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<2>::InPlaceMask();
    }
    case 3: {
      DealFromDP<BitFld<3> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<3>::InPlaceMask();
    }
    case 4: {
      DealFromDP<BitFld<4> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<4>::InPlaceMask();
    }
    case 5: {
      DealFromDP<BitFld<5> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<5>::InPlaceMask();
    }
    case 6: {
      DealFromDP<BitFld<6> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<6>::InPlaceMask();
    }
    case 7: {
      DealFromDP<BitFld<7> >(bkt_buffer, bkt_lens, data_info_p);
      return  BitFld<7>::InPlaceMask();
    }
    case 8: {
      DealFromDP<BitFld<8> >(bkt_buffer, bkt_lens, data_info_p);
      return BitFld<8>::InPlaceMask();
    }
    default: {
      MyAbort("Invalid digit id");
      return bkt_buffer[0];  // No return after abort. Keep compiler happy
    }
    }
  }  // DealAPass
// Do a sort pass from SrcDp to TargetDp. Swap dps after the deal.
  // Return the digit mask for this digit.
  template<typename DBA, typename Data_Info_T>
  void DealFromDP(ElemT *bkt_buffer, int *bkt_lens, 
		  Data_Info_T *data_info_p) {
    const int numBuckets = DBA::NumBuckets();
    TimerId[DBA::FldId] = Timers.ReStartTimer
      (TimerId[DBA::FldId], DBA::FldId);
    
    typedef DealerCoreTemplate<ElemT, Data_Info_T, DBA> DealerCoreT;
    DealerExec<SRec_T, DealerCoreT> Dealer(FlusherP, BatchSize);
    DealerCoreT &DlrCore = Dealer;
    Data_Info_T &DtaInf = DlrCore;
    
    DlrCore.InitDataInfo(*data_info_p, numBuckets);

    bkt_buffer = TargetDp->InitBuffer(bkt_buffer, bkt_lens[numBuckets],
    				      numBuckets);
    FlusherP->InitFlusher(&Dealer, BatchSize, SrcDp->GetBktBuffer(),
			  numBuckets,
			  Timers.GetFlushTimer(TimerId[DBA::FldId]),
			  TargetDp);
    // Requires InitFlusher() before InitBkts()!
    FlusherP->InitBkts(&Dealer, bkt_buffer, bkt_lens);
    
    SrcDp->DealFromDP(&Dealer);
    FlusherP->template Close<DBA>(TargetDp);
    data_info_p->ExtractInfo(Dealer);

#ifdef OCD_TESTING
    TargetDp->template ConfirmTargetBkts<DBA>(SortN);
#endif
    Swap();
    Timers.StopTimer(TimerId[DBA::FldId]);
  }
// Do a sort pass from SrcDp to TargetDp. Swap dps after the deal.
  // Return the digit mask for this digit. Deals directly into counted
  // buckets. bkt_lens contains the exact bkt lens for the deal. Uses
  // CountedNullFlusher flusher class, which essentially has nothing
  // to do.
  template<typename DBA, typename Data_Info_T>
  void DealCountedBktsFromDP(ElemT *bkt_buffer, int *bkt_lens, 
		  Data_Info_T *data_info_p) {
    const int numBuckets = DBA::NumBuckets();
    TimerId[DBA::FldId] = Timers.ReStartTimer
      (TimerId[DBA::FldId], DBA::FldId);
    
    typedef DealerCoreTemplate<ElemT, Data_Info_T, DBA> DealerCoreT;
    DealerExec<SRec_T, DealerCoreT> Dealer(FlusherP, BatchSize);
    DealerCoreT &DlrCore = Dealer;
    Data_Info_T &DtaInf = DlrCore;
    
    DlrCore.InitDataInfo(*data_info_p, numBuckets);

    bkt_buffer = TargetDp->InitBuffer(bkt_buffer, bkt_lens[numBuckets],
    				      numBuckets);
    FlusherP->InitFlusher(&Dealer, BatchSize, SrcDp->GetBktBuffer(),
			  numBuckets,
			  Timers.GetFlushTimer(TimerId[DBA::FldId]),
			  TargetDp);
    // Requires InitFlusher() before InitBkts()!
    FlusherP->InitBkts(&Dealer, bkt_buffer, bkt_lens);
    
    SrcDp->DealFromDP(&Dealer);
    FlusherP->template Close<DBA>(TargetDp);
    data_info_p->ExtractInfo(Dealer);

#ifdef OCD_TESTING
    TargetDp->template ConfirmTargetBkts<DBA>(SortN);
#endif
    Swap();
    Timers.StopTimer(TimerId[DBA::FldId]);
  }
  // Copy the data from a dp to the contiguous target vector. dest_vec
  // is assumed to be large enough to hold the data.
  void CopyDpToVec(ElemT *dest_vec) {
    CopyTimerId = Timers.ReStartTimer(CopyTimerId, "   Copy time = ");
    //std::cout << "Copy data from Dp " << SrcDp << '\n';
    SrcDp->CopyDpToVec(dest_vec);
    Timers.StopTimer(CopyTimerId);
  }
  /**
   * Run the insertion sort on SourceV.  *** Deprecated ***
   **/
  void RunInsertionSort(ElemT top_mask) {
    RunInsertionSort(top_mask, SourceV, SourceV + SortN);
  }
  // Need this version for segmented or recursive sort, where IS may
  // be run on a portion of the total vector.
  void RunInsertionSort(ElemT top_mask, ElemT *src_v, ElemT *src_e) {
    InsertSortTimerId =
      Timers.ReStartTimer(InsertSortTimerId, "\n   Insertion sort time = ");
    //std::cout << "Enter CompleteSort(" << src_v << ", " << src_e << ", top_mask)\n";
    IsObject.CompleteSort(src_v, src_e, top_mask);
    //std::cout << "Exit CompleteSort()\n";
    Timers.StopTimer(InsertSortTimerId);
  }
  void RunCopyInsertionSort(ElemT *dest_vec, ElemT top_mask) { 
   InsertSortTimerId =
      Timers.ReStartTimer(InsertSortTimerId, "   Copy insertion sort time = ");
    // I assume an empty sort vector was sensed long since!
    SrcDp->InsertionSortWithCopy(&IsObject, top_mask, dest_vec);
    Timers.StopTimer(InsertSortTimerId);
  }

private:
  
  ElemT *SourceV;
  long SortN;  // Set to number to sort by DealSrcVector().
  int SrcDpNdx;  // Index the source DataPool. 1 - SrcDpNdx is the destination
  // Required info for sorting.
  DataPoolType DataDefn[2];
  DataPoolType *SrcDp = nullptr;
  DataPoolType *TargetDp = nullptr;
  FlusherT *FlusherP = nullptr;  
  InsertionSortPass<ElemT> IsObject;
  int BatchSize = 0;
  int SetupTimerId = -1;
  int CopyTimerId = -1;
  int InsertSortTimerId = -1;
  
  // If SortSegment <= 1, we're sorting first (or only) segment.
  // DealFromDP() will save timer ids in TimerId[].  Otherwise, we're
  // past the first, and will restart timers with saved timer ids.
  int TimerId[NFields + 1];
  int SortSegment = -1;  
};  // DataPoolPair

#endif

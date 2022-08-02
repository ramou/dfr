/* $Header:  $ */ 

#ifndef DealerExecH
#define DealerExecH

#include "DealerBase.H"
#include "DataInfo.H"

/********************************************************************
 *
 *                            DealerExec.H
 *
 * The dealer gathers the various components that execute a deal.
 *
 * DealerExec::ProcessBlk() is called by the source DataPool once for
 * each contiguous vector of ELEMs in the source. It divides the
 * source stream into batches and calls the Ladler and DataInfo
 * classes to process individual ELEMs. At the end of each batch, it
 * calls the Flusher. It passes the DBA class as a template parameter
 * to the various places that need it.
 *
 * DealerExec::Close() is called by the DataPool. It tells the Flusher
 * to empty all the ladles and have the sink DataPool finish up.
 *
 * A separate dealer is generated for each pass. Template parameters:
 *
 * DataInfo_T : a class to collect info about each element as it is
 *            dealt. Calls nothing.
 *
 * Flusher_t : the class that moves ELEMs from the ladles to the 
 *           DataPool sink. Extends LadleMgr, interacts with DataPool.
 *
 * DBA : (Deal Byte Access) : the class that contains the byte numbers
 *       (wrt the key) of the deal byte and (when needed) the count
 *       byte. Calls nothing. Only static members.
 *
 * DealerExec<> is allocated in DataPoolPair
 ********************************************************************/

// DealerExec                                                    DealerExec

template<typename SRec_T, typename DealerCore_T>
class DealerExec : public DealerCore_T,
//public DataInfo_T,
		   public DealerBase<SRec_T> {

  typedef SRec_T ElemT;
  
  static const int UnwindN = 4;  /////////// Was 8. Needs 4 for slicing
  
public:
  using FlusherT = FlusherTemplate<SRec_T>;
  using DBA = typename DealerCore_T::DBA;
  using DataInfoT = typename DealerCore_T::DataInfoT;
  using DealerCore_T::InitDc;

  static const int DealByte = DBA::DealByte;
  static const int CountByte = DBA::CountByte;

  // No real diff on my desktop (I5).
  static const int RwFlag = 1;  // PF for writing
  //static const int RwFlag = 0;  // PF for reading
  static const int Locality = 3;  // Make it sticky. Far better my desktop
  //static const int Locality = 2;  // Make it somewhat sticky.
  
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  DealerExec(FlusherT *flusher_p, int batch_size) :
    FlusherP(flusher_p), BatchSize(batch_size), RunningBS(batch_size) {}

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  constexpr int GetDealByte() { return DBA::DealByte; }
  
  // Deal a contiguous vector of ELEMs. 
  void ProcessBlk(const ElemT *e_vec, int n) {
    int nTilFlush = RunningBS;
    if (n >= nTilFlush) {  // Complete current batch and continue
      e_vec = UnwindDealer(e_vec, nTilFlush);
      FlusherP->Flush();
      n -= nTilFlush;
      nTilFlush = BatchSize;
      while (n >= nTilFlush) {  
	e_vec = UnwindDealer(e_vec, nTilFlush);
	FlusherP->Flush();
	n -= nTilFlush;
      }
    }
    if (n > 0) {  // Room for all, don't flush yet
      UnwindDealer(e_vec, n);
      nTilFlush -= n;
    }
    RunningBS = nTilFlush;
  }

private:

  // Deal n ELEMs using an unwound loop where possible. n <= batch
  // size, so no flushing here. We do pause to prefetch ladles if
  // enabled.
  NO_INLINE const ElemT *UnwindDealer (const ElemT *e_vec, int n) {
    if (SELECT_DEALER_CORE >= SLICING_DEALER_CORE) {
      return DealerCore_T::UnwindDealer(e_vec, n);
    } else {  // Not SLICING_DEALER_CORE
      GdbHook(0);
      // Make these local so they stay in registers.
      ElemT live = DealerCore_T::DataInfoT::LiveBits;
      const ElemT bm = DealerCore_T::DataInfoT::BitMask;
      const ElemT *eEnd = e_vec + (n - UnwindN);
      const ElemT *dbgEnd = e_vec + n;

      while (e_vec <= eEnd) {
	InnerLoop(e_vec, live, bm);
	e_vec += UnwindN;
      }
      // Now the last few
      eEnd += UnwindN;
      while (e_vec < eEnd) {
	Deal1<0>(e_vec, live, bm);
	e_vec++;
      }
      DealerCore_T::DataInfoT::LiveBits = live;
      return e_vec;
    }
  }

  // Deal UnwindN ELEMs in an unwound loop.
  void InnerLoop(const ElemT *p, ElemT &live, ElemT bm) {
    if (ReadPfAdvance >= 0)
      __builtin_prefetch(p + ReadPfAdvance, 0, 0);
    if ((DealerPfLead == 0) || (DealerPfFreq == 0)) {
      // No pf.
      for (int i = 0; i < UnwindN; i++) Deal1<0>(p + i, live, bm);
    } else if (DealerPfFreq == 1) {
      // Always do pf
      for (int i = 0; i < UnwindN; i++) Deal1<DealerPfLead>(p + i, live, bm);
    } else {  // do pf every 2nd store
      for (int i = 0; i < UnwindN; i+=2) {
	Deal1<DealerPfLead>(p + i, live, bm);
	Deal1<0>(p + i+1, live, bm);
      }
    }
    if (DEAL_PF_FAR > 0) {  // PF with a bigger lead at end of each loop.
      int bktId = FarBktID;
      __builtin_prefetch(DealerCore_T::GetStoreAddr(bktId) + DEAL_PF_FAR,
			 RwFlag, Locality);
      FarBktID = (++bktId) & (MaxNBuckets-1);  //*** Assume NBuckets is 2**n!!
    }
  }

  // Deal a single ELEM. Dlr_Pf_Lead allows us to omit pf
  // sometimes. Ignored if DealerCore_T::Deal1Template == 0.
  template<int Dlr_Pf_Lead>
  void Deal1(const ElemT *p, ElemT &live, ElemT bm) {
    ElemT t = *p;
    DealerCore_T:: template Deal1<Dlr_Pf_Lead>(t, DBA::GetDealB(p));
    DataInfoT::template Collect<DBA>(t, live, bm);
  }

  FlusherT *FlusherP = nullptr;
  int BatchSize = 0;  
  int RunningBS = 0;  // N left to fill current batch.
  int FarBktID = 0;
  
};  // DealerExec

#endif

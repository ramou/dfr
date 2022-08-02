/* $Header:  $ */ 

#ifndef OverflowH
#define OverflowH

#include "TimeWrapper.H"
#include "AlignedBlock.H"
#include "MuAligner.H"
#include "BigOverflow.H"

/********************************************************************
 *
 *                            Overflow.H
 *
 * Define a class to manage the 3 overflow areas and initialize the
 * bucket overflow area vectors that use them. 
 *
 * Like the data pools, the buffers are CL aligned and all accesses
 * can make assumptions about block alignment as long as requested
 * blocks have appropriate lengths.
 ********************************************************************/

// Overflow                                                       Overflow
template<typename Elem_T, int N_Buckets>
class Overflow {
public:
  typedef FlushMovers::MoverClass<Elem_T, SELECT_FLUSH_MOVER> Mover_T;
  typedef typename Mover_T::MuT MuT;
  typedef MuAligner<Elem_T, MuT> Aligner;
  
  static const int RtPerMu = Mover_T::RtPerMu;
  static const int NBuckets = N_Buckets;
#ifndef OF_COMPRESS_THRESHOLD
#define OF_COMPRESS_THRESHOLD 0.4
#endif
  
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  // src_n just gives a limit that can't be overrun by overflows when
  // we have moved from the preface area to the primary area (which is
  // the source vector).
  Overflow(CountType src_n, Mover_T *mover_p = nullptr) :
    SafeN(src_n), MoverP(mover_p) {}
  //~Overflow();

  /**
   * over_flow_safe points to an area which can't overflow, to be used
   * as the primary overflow buffer this pass. safety_margin is the
   * amount of excess to allocate to the preface in addition to the
   * minimum requirement. When it is batch size, there is always
   * enough space to hold the entire batch's overflow. Thus we need
   * only check for end of the preface once per batch, and once in
   * Close() before flushing leftovers.
   *
   * We round over_flow_safe up to a CL boundry so we can safely use
   * aligned simd loads and stores, and add 2 * RtPerMu to
   * safety_margin to make up for shortening the Over_Flow_Safe buffer
   * at both ends. If safety_margin < (256 * NBuckets), it is
   * increased to that value since we may allocate incomplete chunks
   * of up to RtPerMu during Close().
   **/
  void OverflowInit(Elem_T *over_flow_safe, int safety_margin) {
    if (safety_margin < (256 * NBuckets)) safety_margin = 256 * NBuckets;
    safety_margin += 2 * RtPerMu;
    over_flow_safe = Aligner::AlignMuTUp(over_flow_safe);
    OverflowInit_imp(over_flow_safe, safety_margin);
  }
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  /**
   * Interface defns used in several functions:
   *
   * Elem_T** bkt_addr : bkt_addr[b] points to primary space for bkt[b].
   * int *bkt_n : bkt_n[b] is # recs dealt to primary bkt[b].
   * Elem_T** ofbkt_addr : ofbkt_addr[b] points to overflow for bkt[b].
   * int* ofbkt_n : ofbkt_n[b] is # overflowed recs for bkt[b].
   * int* of_bkt_list : list of bkt numbers of overflowing bkts
   * int n_ofbkts : number of overflowing buckets (==size of of_bkt_list).
   * pool_p : pointer to the target DataPool (type is template parm).
   * int n_buckets : number of buckets being dealt.
   **/
  /**
   * Ensure that there is still a safety margin in the current
   * overflow area.
   **/
  void CheckOfArea() {
    if (CurAreaEnd >= CurAreaNext) return;
    // Must be preface, and is full.
    PrefaceEnd = CurAreaNext;
    CurAreaNext = SafeAreaOrigin;
    CurAreaEnd = SafeAreaOrigin + SafeN;
    UsingPreface = false;
  }
  /**
   * Allocate a block of size n_elems and return its tail address. 
   **/
  Elem_T *GetOfBlockEnd(int n_elems) {
#if 0
    std::cout<< ", GetOfBlockEnd: N=" << n_elems << " start =" <<
      CurAreaNext << "\n";
#endif
    return CurAreaNext += n_elems;
  }
  /**
   * Allocate a block of size n_elems and return its addr
   **/
  Elem_T *GetOfBlockP(int n_elems) {
#if 0
    std::cout<< ", GetOfBlockP: N=" << n_elems << " start =" <<
      CurAreaNext << "\n";
#endif
    Elem_T *p = CurAreaNext;
    CurAreaNext = p + n_elems;
    return p;
  }
  void Store1(Elem_T &v) { *CurAreaNext++ = v; }
  /**
   * Dealing for this pass is done, and overflow count for each bucket
   * is available. Allocate the bucket overflow areas as needed,
   * compute the addrs of the overflow chunks, deal the overflowing
   * recs into them, and create the primary and overflow chunks in
   * DataPool pool_p.
   *
   * of_bkt_list : a vector of bucket numbers for all overflowing buckets
   * n_ofbkts : number of overflowing buckets.
   * bkt_addr: Addr of start of primary buckets.
   * bkt_n:    # of items dealt to bkt[b] before OF.
   * of_counts : ~[b] is number overflowing for bucket[b]. Undefined
   *   when b is not in of_bkt_list.
   * n_ofe : total number of overflowing Elems.
   * align_size : the alignment requirement.
   * pool_p : pointer to the target DataPool
   **/
  template<typename DBA, typename Pool_T>
  void CreateChunks(int *of_bkt_list, int n_ofbkts,
		    Elem_T** bkt_addr, CountType *bkt_n,
		    int *of_counts, long n_ofe, int align_size,
		    Pool_T *pool_p) {
    AllocOfBkts(of_bkt_list, n_ofbkts, bkt_addr, bkt_n,
		of_counts, n_ofe, align_size);
    // Must do before copy, since FillOfBuckets() changes OFBktV.
    CreateDpChunks(bkt_addr, bkt_n, OFBktV, of_counts, pool_p);
    FillOfBuckets<DBA>(align_size);
    //pool_p->CountRecs(SafeN, "At end of CreateChunks: ");  // Counts good here.
  }
private:

  /**
   * Move the data from the ungrouped preface and primary overflow
   * areas to the per bucket areas. If align_size > 1, all blocks are
   * assumed to be multiples of align_size. This can speed the
   * scanning and copying operations.  Parameters:
   *
   * DBA: (template) Deal Byte Access class.
   * align_size: the alignment requirement.
   **/
  template<typename DBA>
  void FillOfBuckets(int align_size) {
    if (PrefaceEnd > PrefaceOrigin) {
      Move2Bkts<DBA>(PrefaceOrigin, PrefaceEnd, align_size);
    }
    if (SafeAreaEnd > SafeAreaOrigin) {
      Move2Bkts<DBA>(SafeAreaOrigin, SafeAreaEnd, align_size);
    }
  }
  //////////////////////// Data members
  // Initialized at the start of the pass. 
  
  Elem_T *MergedOrigin;
  long MergedUsed { 0 };
  // Total n in compressed bkt space. Used when sliding primary bucket
  // to one end of the buffer to make room for the OF recs in the main
  // buffer.
  long BktTotalN=0;

  Elem_T *SafeAreaOrigin;
  Elem_T *SafeAreaEnd; // End of used part of primary area

  Elem_T *PrefaceOrigin;
  Elem_T *PrefaceEnd;  // End of used part of preface area

  Elem_T *CurAreaNext;  // Where to drop next item(s).
  Elem_T *CurAreaEnd;   // Nominal end of currently active OF area.

  CountType SafeN;  // Yes. is src_n.
  bool UsingPreface;

  // Allow copy in Move2Bkts to reference the flush mover if this is
  // set. Init at construction.
  Mover_T *MoverP;

  BigOverflow<Elem_T, N_Buckets, RtPerMu> ReorderBkts;
  
  /**
   * When overflows exist they must be grouped by bucket and moved to
   * MergedOrigin at the end of the pass. OFBktV[b] will be set to
   * the addr of each overflow bucket. If bucket[b] hasn't overflowed,
   * it is nullptr. We need a vector of addresses of the starts, to be
   * returned to the LadleManager, and a running index within each
   * bucket(??).
   **/
  Elem_T *OFBktV[N_Buckets];

#include "Overflow.C"
  
};  // Overflow

#endif

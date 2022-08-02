/* $Header:  $ */ 

#ifndef FlushStatusH
#define FlushStatusH

/********************************************************************
 *
 *                            FlushStatus.H
 *
 * Gather the values that define the current status for a possibly
 * flushing ladle and its bucket. Functions to use that data for
 * prefetching and flushing.
 *
 * NOTE: Parallel vectors may be better, since LadleStart is constant
 * after init for each pass, while the rest vary. Also precludes any
 * simd ops.
 ********************************************************************/

// FlushStatus                                                   FlushStatus
template<typename SRec_T, typename Bkt_Mu_T, int Rt_Per_Mu>
class FlushStatus {
  // Addr of start of this ladle. Const after init. 
  SRec_T *LadleStart;
  // Addr of next free slot in bucket. If overflowing, addr of
  // overflow slot to use. Incremented on call to flush mover. 
  union {
    SRec_T *BktAvail;
    Bkt_Mu_T* MuBktAvail;
  };
  short IlLength;  // of the ladle
  short IlLeftovers;

public:
  typedef IlElem<SRec_T, IlN> IlElemT;
  static constexpr int LeftoverMask = (Rt_Per_Mu - 1) * IlN;
  
  /////////////////////////////////////////////////////////////////
  //                Initialization.

  void InitLdlAddr(SRec_T *ldl_p, SRec_T *bkt_p) {
    LadleStart = ldl_p;
    BktAvail = bkt_p;
  }

  //////////////////////////////////////////////////////////////////
  //               Readers / writers

  // This "init" is really used to store target overflow address.
  void InitBktAvail(SRec_T *p) { BktAvail = p; }

  // Get addr of the ladle as a ptr to record, to an Il record, or an Mu.
  SRec_T *Ladle() const { return LadleStart; }
  IlElemT* IlLadle() const { return reinterpret_cast<IlElemT*>(LadleStart); }
  template<typename Ldl_Mu_T>
  Ldl_Mu_T *LdlMu() const { return reinterpret_cast<Ldl_Mu_T*>(LadleStart); }

  // As above, but the addr returned is after adding il_ndx to
  // LadleStart.  To use a bare (IlN not included) index, use the
  // above functions and add it to the desired form of ladle address.
  SRec_T *Ladle(int il_ndx) const { return LadleStart + il_ndx; }
  IlElemT* IlLadle(int il_ndx) const
  { return reinterpret_cast<IlElemT*>(LadleStart + il_ndx); }
  template<typename Ldl_Mu_T>
  Ldl_Mu_T *LdlMu(int il_ndx) const
  { return reinterpret_cast<Ldl_Mu_T*>(LadleStart + il_ndx); }
  
  // Get bucket free space as a tail pointer to record or Mu.
  // BktAvail is updated.
  SRec_T* BktEnd(int len) {
    return BktAvail += len;
  }
  
  Bkt_Mu_T *MuBktEnd(int len) {
    BktEnd(len);
    return MuBktAvail;
  }

  static Bkt_Mu_T *ToBktMuT(SRec_T *p)
  { return reinterpret_cast<Bkt_Mu_T*>(p); }
  
  // Get start of unused rest of bucket.
  SRec_T* Bkt() const { return BktAvail; }
  Bkt_Mu_T *MuBkt() const { return reinterpret_cast<Bkt_Mu_T*>(BktAvail); }

  // Store the current il len, including leftover count.
  void StoreIlLen(int il_cur_n) {
    int t = il_cur_n & LeftoverMask;
    IlLeftovers = t;
    IlLength = il_cur_n - t;
  }
  // Get Il length, plain length, Il leftovers, or plain leftovers
  int IlLen() const { return IlLength; }
  int Len() const { return IlLength / IlN; }
  int IlLeftos() const  { return IlLeftovers; }
  int Leftos() const { return IlLeftovers / IlN; }
  
  //////////////////////////////////////////////////////////////////
  //               Other functions

  /**
   * Flush the ladle to BktAvail and inc BktAvail. Ignores leftovers.
   * mover_p points to the mover class.
   **/
  template<typename Mover_T>
  void FlushToMvr(Mover_T * mvr_p, int flush_n) {
    mvr_p->FlushMus(IlLadle(), IlLadle() + flush_n, MuBktEnd(flush_n),
		    flush_n);
  }
  /**
   * As above, but also slides leftover items to start of ladle. Most
   * flushes done through here.  mover_p points to the mover class.
   **/
  template<typename Mover_T>
  void FlushLdl(Mover_T * mvr_p) {
    FlushToMvr(mvr_p, Len());
    if (Rt_Per_Mu > 1)  // Else can't be any leftovers.
      SlideTail();
  }
  /**
   * Flush ladle[start_ndx : start_ndx + len] to bkt_end_p. parm's
   * don't include IlN. 
   **/
  template<typename Mover_T>
  void Shorten(Mover_T * mvr_p, int start_ndx, int len,
               Bkt_Mu_T *bkt_end_p) const {
    mvr_p->FlushMus(IlLadle() + start_ndx, IlLadle() + start_ndx + len,
                    bkt_end_p, len);
  }
  template<typename Mover_T>
  Bkt_Mu_T* ShortenLo(Mover_T * mvr_p, int start_ndx, int len,
               Bkt_Mu_T *bkt_end_p) const {
    mvr_p->template Shorten<true>(IlLadle() + start_ndx, bkt_end_p, len);
  }
  template<typename Mover_T>
  Bkt_Mu_T* ShortenHo(Mover_T * mvr_p, int start_ndx, int len,
               Bkt_Mu_T *bkt_end_p) const {
    mvr_p->template Shorten<false>(IlLadle() + start_ndx, bkt_end_p, len);
  }
  // Slide leftovers down to ladle start.
  void SlideTail() const {
    if (Rt_Per_Mu == 1) return;  // Never any tail to flush
    int nLo = IlLeftos();
    SRec_T* srcP = Ladle(IlLen());
    SRec_T *destP = Ladle();
    // for written to allow unwinding (with internal exit). Does it???
    for (int i = 0; i < Rt_Per_Mu * IlN; i+=IlN) {
      if (i >= nLo) return;
      destP[i] = srcP[i];
    }
  }

  ///////////////// Prefetchers
  /**
   * Prefetch bucket space to flush this ladle.
   **/
  void PfBktSpace() {
    if (BktNToPrime >= 0) {
      const SRec_T *p = BktAvail;
      RevPrefetchVec<1, BktLocality>(p, p + Len() - 1);
    }
  }
  /**
   * Prefetch space holding il_len recs. We don't need the exact ladle addr.
   **/
  void PfLadle(int il_len) {
    if (LdlNToPrime >= 0) {
      const SRec_T *p = Ladle();
      RevPrefetchVec<0, 3>(p, p + il_len - 1);
    }
  }
};  // FlushStatus

#endif

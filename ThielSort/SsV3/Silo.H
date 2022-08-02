/* $Header:  $ */ 

#ifndef SiloH
#define SiloH

/********************************************************************
 *
 *                            Silo.H
 *
 * Doing a read prefetch ahead of processing the source data when it
 * is scattered (i.e., LeapFrog), means following the chain
 * twice. This adds both complexity and an additional load on the
 * cache. This object helps this problem by storing the block address
 * and length and doing the prefetch separately.
 *
 * The silo is just a simple circular vector that manages its own
 * subscripting so the processing loop doesn't need to. It is
 * templated on data type and silo size, which must be > the look
 * ahead distance. The size of the silo is 2**N_Power.
 ********************************************************************/

// Silo                                                       Silo
template<typename ELEM_T, int N_Power>class Silo {

public:

  static const int Size = 1 << N_Power;
  
  struct BlockSpec {
    ELEM_T *Addr;
    long Len;
    void Store(ELEM_T *p, int len) { Addr = p; Len = len; }
  };
  
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  Silo() { CurStore = CurRead = 0; BaseAddr = nullptr; }
  Silo(ELEM_T *base_addr, int look_ahead) { Reset(base_addr, look_ahead); }
  //~Silo();

  /**
   * Reset silo to empty. base_addr != nullptr is the address of the
   * vector that the ints passed to PfAndStore() index. look_ahead is
   * the look ahead distance for the prefetch. Abort if >= Size.
   **/
  void Reset(ELEM_T *base_addr, int look_ahead) {
    CurStore = CurRead = 0;
    BaseAddr = base_addr;
    if (base_addr == nullptr) MyAbort("Bad base_addr in Silo::Reset()");
    if (look_ahead >= Size) MyAbort("Bad look_ahead in Silo::Reset()");    
  }
  
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  /**
   * Prefetch (locality 0) the block. Save the addr and n so GetNext()
   * can return them. n > 0.
   **/  
  void PfAndStore(int start_ndx, int len) {
    ELEM_T *p = BaseAddr + start_ndx;
    SiloV[CurStore].Store(p, len);
    CurStore = (CurStore + 1) & (Size - 1);
    __builtin_prefetch(p, 0, 0);
    char *curCL = reinterpret_cast<char *>(p) + CacheLineLen;
    // PrefetchVector here seems much worse. Not in lined?
    //PrefetchVector<1, 0>(p, len);
#if 1
    char *endCL = curCL + sizeof(ELEM_T) * len;
    for (; curCL < endCL; curCL += CacheLineLen) {
      __builtin_prefetch(curCL, 0, 0);
    }
#endif
  }
  /**
   * If CurRead == CurStore, return nullptr, indicating end of
   * data. Otherwise prefetch (locality 0) each CL in the next block and
   * return a pointer to its BlockSpec. If do_pf is false, we skip the
   * prefetch.
   **/
  BlockSpec *GetNext(bool do_pf = true) {
    if (__builtin_expect(CurStore == CurRead, 0)) return nullptr;
    BlockSpec *p = SiloV + CurRead;
    CurRead = (CurRead + 1) & (Size - 1);
    return p;
  }
  
private:
  BlockSpec SiloV[Size];
  ELEM_T *BaseAddr;
  // Indices into Vec.
  int CurStore = 0;
  int CurRead = 0;
  // 
};  // Silo

#endif

/* $Header:  $ */ 

#ifndef DealerCoresH
#define DealerCoresH

/********************************************************************
 *
 *                            DealerCores.H
 *
 * Several cases of dealer cores (tiny classes that deal a single item
 * at a time). The correct one is picked by Config.
 *
 * Later versions (added with #includes at end of this file) provide
 * DealerCores with more complex loop unwinding to try to optimize
 * cache and register use to reduce pipeline stalls.
 ********************************************************************/

//////////// NOTE: Inc of store addr/index after use seems to add an
////////////       extra mov instruction: Local value copied to
////////////       another register before inc and store. Only a
////////////       register move, but...

//////////// If DC_POST_INC == 1, inc after use. If == 0, inc before use.
#ifndef DC_POST_INC
#define DC_POST_INC 1
#endif

// Generic DealerCore class template. Doesn't exist.
template<typename SRec_T, typename DataInfo_T, typename DBA_T,
	 int Selector_Id> class DealerCoreClass;

/**
 * Since we have distinct init calls for dealercore variations, we use
 * a base class to send a generic error message if the wrong init is
 * ever called. It won't be, because a compile time if prevents it,
 * but the compiler still needs to see the target:-(( .
 **/
class DealerCoreError {
protected:
  [[noreturn]]
  void Error() { MyAbort("Invalid call for this DealerCore variant."); }

public:
  void InitDc(void **ldl_addrs, void *ldl_space, int *ldl_ndx) 
  { Error(); }
  //template<typename void> 
  void InitDc(void **bkt_addrs, Overflow<void, MaxNBuckets> *of_area,
	      const int *bkt_lens) { Error(); }
  template<typename SRec_T>
  const SRec_T *UnwindDealer (const SRec_T *e_vec, int n)
  { Error(); return nullptr; }  
  
 };
// DirectStore                                                 DirectStore
//////////////////////////////////////////////////////////////////
// Declare DealerCore for direct store ladles, no prefetch.
template<typename SRec_T, typename DataInfo_T, typename DBA_T>
class DealerCoreClass<SRec_T, DataInfo_T, DBA_T,
		      DIRECT_STORE_DEALER_CORE> :
  public DealerCoreError,
  public DataInfo_T,
  public DBA_T {
  
  // DirectStore
  SRec_T *LadleSpace;
  // Use LdlAddr if WideIndices, LdlNdx otherwise.
  SRec_T **LdlAddr;  // Addr to store next SRec_T
  int *LdlNdx;      // Index to store next SRec_T
    
public:
  using DBA = DBA_T;
  using DataInfoT = DataInfo_T;

  SRec_T *GetStoreAddr(int b) {
    if (WideIndices) return LdlAddr[b];
    else return LadleSpace + LdlNdx[b];
  }
  void InitDc(SRec_T **ldl_addrs, SRec_T *ldl_space,
	      int *ldl_ndx) {
    LdlAddr = ldl_addrs; LdlNdx = ldl_ndx; LadleSpace = ldl_space;
    static_assert((SELECT_DEALER_CORE != DIRECT_STORE_DEALER_CORE) ||
		  (DEALER_PF_LEAD == 0),
		  "Config requests DEALER_PF_LEAD > 0 with "
		  "DIRECT_STORE_DEALER_CORE");
  }
  // Basic Deal1(). No prefetch.
  template<int Dlr_Pf_Lead>void Deal1(SRec_T v, int bkt_id) {
    if (WideIndices) { *LdlAddr[bkt_id] = v; LdlAddr[bkt_id] += IlN; }
    else { LadleSpace[LdlNdx[bkt_id]]  = v; LdlNdx[bkt_id] += IlN; }
  }
  static const char *WhoAmI() {
    return WideIndices ? "WDirect" : "NDirect";
  }
  
};  // DirectStore
//////////////////////////////////////////////////////////////////////////
// DirectStoreWithPf                                     DirectStoreWithPf
//Declare DealerCore for direct store ladles, with prefetch option.
template<typename SRec_T, typename DataInfo_T, typename DBA_T>
class DealerCoreClass<SRec_T, DataInfo_T, DBA_T,
		      DIRECT_STORE_DEALER_CORE_WITH_PF> :
  public DealerCoreError,
  public DataInfo_T,
  public DBA_T {
  
  // DirectStoreWithPf
  SRec_T *LadleSpace;
  // Use LdlAddr if WideIndices, LdlNdx otherwise.
  SRec_T **LdlAddr;  // Addr to store next SRec_T
  int *LdlNdx;      // Index to store next SRec_T
  // No real diff on my desktop (I5).
  static const int RwFlag = 1;  // PF for writing
  //static const int RwFlag = 0;  // PF for reading
  static const int Locality = 3;  // Make it sticky. Far better my desktop
  //static const int Locality = 2;  // Make it somewhat sticky.
    
public:
  using DataInfoT = DataInfo_T;
  using DBA = DBA_T;

  SRec_T *GetStoreAddr(int b) {
    if (WideIndices) return LdlAddr[b];
    else return LadleSpace + LdlNdx[b];
  }

  void InitDc(SRec_T **ldl_addrs, SRec_T *ldl_space,
	      int *ldl_ndx) {
    LdlAddr = ldl_addrs; LdlNdx = ldl_ndx; 
#if (DC_POST_INC == 1)
    LadleSpace = ldl_space;
#else
    LadleSpace = ldl_space - IlN;
#endif
  }
  // Deal1() that does prefetch Dlr_Pf_Lead SRec_Ts past initial
  // store point.  0 == NO prefetch.
  template<int Dlr_Pf_Lead>void Deal1(SRec_T v, int bkt_id) {
    /**
     * NOTE: At some point DEALER_PF_LEAD was set to default to
     * IL_VALUE, which was wrong when IL_VALUE > 1 was also used to
     * adjust RealPfLead here. Some timings will reflect this
     * bug. Fixed 2021-11-18.
     **/
    /**
     * NOTE: The value of RealPfLead must reflect DC_POST_INC, since
     * the pf is done using an inc'ed addr/index. Thus, if the inc is
     * done after the store, the pf addr is already leading by 1, and
     * RealPfLead is added on top of that, effectively adding 1 (in
     * the bucket, not the bin) to RealPfLead. But if the inc is done
     * before the store, RealPfLead is the real lead. So must subtract
     * 1 from Dlr_Pf_Lead when we comput RealPfLead if DC_POST_INC is
     * 0. This was done at the pf before today, meaning the
     * Dlr_Pf_Lead was effectively reduced by 1 when DC_POST_INC was 1.
     **/
    
    constexpr int RealPfLead = (Dlr_Pf_Lead - DC_POST_INC)  * IlN; 
    if (WideIndices) {
#if 0
      // This make a miserably slow deal. I'm using it wrong or it
      // just doesn't work well. It is also not flexable: Only good
      // for long long int. Doesn't compile even for unsigned!
      // Seems unusable.
      if (NtDeal) {
	_mm_stream_si64(LdlAddr[bkt_id], v);
      } else {
	*LdlAddr[bkt_id] = v;
	if (Dlr_Pf_Lead > 0)
	  __builtin_prefetch(LdlAddr[bkt_id]+RealPfLead, RwFlag, Locality);
      }
#else
      *LdlAddr[bkt_id] = v;
      if (Dlr_Pf_Lead > 0)
	__builtin_prefetch(LdlAddr[bkt_id]+RealPfLead, RwFlag, Locality);
#endif
      LdlAddr[bkt_id] += IlN;
    } else {
#if (DC_POST_INC == 1)
      LadleSpace[LdlNdx[bkt_id]]  = v;
      LdlNdx[bkt_id] += IlN;
#else
      long t = LdlNdx[bkt_id] + IlN;
      LdlNdx[bkt_id] = t;
      LadleSpace[t]  = v;
#endif      
      if (Dlr_Pf_Lead > 0)
	// Note that there is no hidden overhead: Only the prefetcht0 instr
	// is added. (Checked code with pf on alternating stores.)
	__builtin_prefetch(LadleSpace+LdlNdx[bkt_id]+RealPfLead,
			   RwFlag, Locality);
    }
  }
  static const char *WhoAmI() {
    return WideIndices ? "WDirect, PF" : "NDirect, PF";
  }
};  // DirectStoreWithPf
//////////////////////////////////////////////////////////////////

// Exclude variants
//#include "UnwindDC.H"
//#include "DcPfMeNow.H"

// NOTE: defaults provide access to WhoAmI(), which is dependent only
// on SELECT_DEALER_CORE.
template<typename SRec_T, typename DataInfo_T = NoopDataInfoT<SRec_T>,
	 typename DBA_T = DealByteAccess<SRec_T, 1, 2>>
using DealerCoreTemplate = DealerCoreClass
  <SRec_T, DataInfo_T, DBA_T, SELECT_DEALER_CORE>;

#endif

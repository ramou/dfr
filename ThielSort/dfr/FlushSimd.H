/* $Header:  $ */ 

#ifndef FlushSimdH
#define FlushSimdH
  
/********************************************************************
 *
 *                            FlushSimd.H
 *
 * Move data using simd NT instructions. 
 * See FlushMovers.H for specs.
 ********************************************************************/
#include <immintrin.h>

#include "SimdMove.H"

// FlushSimd                                                       FlushSimd
template<typename SRec_T>
class MoverClass<SRec_T, FLUSH_SIMD_MOVER> {  // FlushSimd
public:
  // Following 2 types are the same. Name diff emphasizes function.
  typedef SimdUnitMover::RemoteSdu<SIMD_MOVE_SIZE> MuT;
  typedef SimdUnitMover::RemoteSdu<4> MuT128; 
  typedef SimdUnitMover::LocalSdu<SIMD_MOVE_SIZE> LocalSdu;
  typedef LocalSdu LdlMuT;
  typedef SimdUnitMover::LocalSdu<4> LocalSdu128;
  
  typedef IlElem<SRec_T, IlN> IlElemT;
    
  // Assumes sizeof(SRec_T) plays nicely with sizeof(simd)!
  static constexpr int RtPerSimd = MuT::SimdLen / sizeof(SRec_T);
  static constexpr int RtPerSimd128 = MuT128::SimdLen / sizeof(SRec_T);
  
  static constexpr int RtPerMu = FindRtPerMu(RtPerSimd);
  static constexpr int LogRtPerMu = FindPower(RtPerMu);
  
  // Unwind loops by UnwindN MUs (which may not mean UnwindN simds!)
  static constexpr int UnwindN = 4;

  static void staticassert() {
    // Just a place to put these
    static_assert(RtPerMu >= RtPerSimd, "Config error with RtPerMu");
  }

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  //FlushSimd();
  //~FlushSimd();
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  void FlushMus(IlElemT *ldl_p, IlElemT * ldl_end_p, MuT* block_end_p,
                int n_e) {
    using AlignerT = MuAligner<SRec_T, MuT>;
#ifdef OCD_TESTING
    if (!AlignerT::MuAligned(block_end_p))
	ErrTag("FlushSimd alignment error", true);
#endif
    if (IlN == 1)
      FlushMusNoIl(reinterpret_cast<LdlMuT*>(ldl_p),
		   reinterpret_cast<LdlMuT*>(ldl_end_p),
		   block_end_p, n_e / RtPerSimd);
    else {
      if (LoadSimdReg < 4) FlushMusIl(ldl_p, ldl_end_p, block_end_p, n_e);
      else if (LoadSimdReg == 4)
	FlushMusIlSimdGthr(ldl_p, ldl_end_p, block_end_p, n_e);
      else if (LoadSimdReg == 5)
	FlushMusIlNT16(ldl_p, ldl_end_p, block_end_p, n_e);

    }
  }
  /**
   * Copy n_mus MUs from src_p to dest_p. Pointers must respect NT
   * simd alignment requirements.
   **/  
  static void CopyMus(const SRec_T *src_p, SRec_T *dest_p, int n_mus) {
    // Cause I can't be bothered to propagate the const through a
    // library hierarchy now. I should fix that eventually.
    SRec_T *p = const_cast<SRec_T *>(src_p);
    FlushMusNoIl(reinterpret_cast<LdlMuT*>(p),
		 reinterpret_cast<LdlMuT*>(p) + n_mus,
		 reinterpret_cast<MuT*>(dest_p) + n_mus, n_mus);
  }
  // As above, but n_e is number of elements, not MUs. We convert to
  // MUs and **ignore** leftovers.
  static void CopyRecs(const SRec_T *src_p, SRec_T *dest_p, int n_e) {
    CopyMus(src_p, dest_p, n_e >> LogRtPerMu);
  }
#if 1
  typedef __m128i GthrNdxSduT;  // Type of gather index simd unit.
  typedef int GthrNdxT;         // Type of gather index.

  typedef long long int IntelLL;  // Intel/g++ decl for simd intrinsic
  typedef const long long int IntelCLL;  // and const version
  
  static const int IntPerGthrNdx = sizeof(GthrNdxSduT) / sizeof(GthrNdxT);

  static const char *WhoAmI() {
    const char* VersionId[] =
      {"FlushSimd(v3:no n_e)SimdLd1",
       "FlushSimd(v3:no n_e)SimdLd2",
       "FlushSimd(v3:no n_e)SimdLd3",
       "FlushSimd(v3:no n_e)SimdLd4",
       "FlushSimd(v3:no n_e)SimdLd5",
       "FlushSimd(v3:no n_e)SimdLd6"};
    
    if (LoadSimdReg > 5) ErrTag("LOAD_SIMD_REG config error", true);
    if ((LoadSimdReg >= 4) && (sizeof(int64_t) != sizeof(SRec_T)))
      ErrTag("Size of SRec_T != 8 in FlushSimd", true);
    return VersionId[LoadSimdReg - 1];
  }

  private:

// Version for IlN > 1, using 16 byte NT stores, since so expensive to load 
  // scattered data into a 32 bit register. Is this any better??
  NO_INLINE
  void FlushMusIlNT16(IlElemT *ldl_p, IlElemT * ldl_end_p, MuT* block_end_p,
		      int n_e) {
    MuT128* blockEndP = reinterpret_cast<MuT128*>(block_end_p);
    FlushMusByMum128(ldl_p + RtPerMu * (UnwindN-1), ldl_end_p, blockEndP);
    // 0:(UnwindN-1) Mus left.
    IlElemT* p = ldl_end_p - 2 * RtPerMu;
    // Assumes UnwindN == 4!
    if (ldl_p <= p) {  // Move 2 more Mus. 
      for (int i = 0; i < 2 * RtPerMu / RtPerSimd128; i++) {
	ldl_end_p -= RtPerSimd128;
	LocalSdu128 t; t. template Set(ldl_end_p[1], ldl_end_p[0]);
	blockEndP--;
	*blockEndP = t;
      }
    }
    if (ldl_p < ldl_end_p) {  // 1 more to go
      for (int i = 0; i < RtPerMu / RtPerSimd128; i++) {
        ldl_end_p -= RtPerSimd128;
        LocalSdu128 t; t. template Set(ldl_end_p[1], ldl_end_p[0]);
        blockEndP--;
        *blockEndP = t;
      }
    }

  }

  
// Move sets of UnwindN Mus from ldl_end_p[-i] to block_end_p[-i]
  // until ldl_end_p <= ldl_p.
  NO_INLINE void FlushMusByMum128(IlElemT *ldl_p, IlElemT *&ldl_end_p,
                        MuT128* &block_end_p) {
    static constexpr int UnwSimdN = UnwindN * RtPerMu / RtPerSimd128;
    GdbHook(1);
    while (ldl_p < ldl_end_p) {
      for (int i = -1; i >= -UnwSimdN; i--) {
	ldl_end_p -= RtPerSimd128;
	LocalSdu128 t; t. template Set(ldl_end_p[1], ldl_end_p[0]);;
	block_end_p[i] = t;
      }
      block_end_p -= UnwSimdN;
    }
  }

// Version for IlN == 1, so just load the simd units.
#if 0
  // This version (NT loads from ladles) seems worse than normal loads.
  NO_INLINE
  void FlushMusNoIl(LdlMuT *ldl_mu_p, LdlMuT * ldl_mu_end_p, MuT* block_end_p,
		    int n_mu) {
    // This approach gets rid of the mucking around to exit loop with
    // a test for eq / neq.
    static_assert(UnwindN == 4, "UnwindN mismatch in Flush1By1");
    if (n_mu & 2) {
      if (n_mu & 1) {  // 3 to do
	LdlMuT t4; t4.NtMoveIn(ldl_mu_end_p-1);
	LdlMuT t5; t5.NtMoveIn(ldl_mu_end_p-2);
	LdlMuT t6; t6.NtMoveIn(ldl_mu_end_p-3);
	block_end_p[-1] = t4;
	block_end_p[-2] = t5;
	block_end_p[-3] = t6;
	block_end_p -= 3;
	ldl_mu_end_p -= 3;
      } else {  // only 2
	LdlMuT t4; t4.NtMoveIn(ldl_mu_end_p-1);
	LdlMuT t5; t5.NtMoveIn(ldl_mu_end_p-2);
	block_end_p[-1] = t4;
	block_end_p[-2] = t5;
	block_end_p -= 2;
	ldl_mu_end_p -= 2;
      }
    } else if (n_mu & 1) {  // 1 more to go
      LdlMuT t6; t6.NtMoveIn(ldl_mu_end_p-1);
      block_end_p[-1] = t6;
      block_end_p--;
      ldl_mu_end_p--;
    }
    if (ldl_mu_p == ldl_mu_end_p) return;
    //FlushMusByMuNoIl(ldl_mu_p, ldl_mu_end_p, block_end_p);
    GdbHook(1);
    ldl_mu_end_p -= UnwindN;
    LdlMuT t0; t0.NtMoveIn(ldl_mu_end_p+3);
    LdlMuT t1; t1.NtMoveIn(ldl_mu_end_p+2);
    LdlMuT t2; t2.NtMoveIn(ldl_mu_end_p+1);
    LdlMuT t3; t3.NtMoveIn(ldl_mu_end_p+0);
    while (ldl_mu_p != ldl_mu_end_p) {
      ldl_mu_end_p -= UnwindN;
      block_end_p -= UnwindN;
      block_end_p[3] = t0;
      t0.NtMoveIn(ldl_mu_end_p+3);
      block_end_p[2] = t1;
      t1.NtMoveIn(ldl_mu_end_p+2);
      block_end_p[1] = t2;
      t2.NtMoveIn(ldl_mu_end_p+1);
      block_end_p[0] = t3;
      t3.NtMoveIn(ldl_mu_end_p+0);
    }
    block_end_p -= UnwindN;
    block_end_p[3] = t0;
    block_end_p[2] = t1;
    block_end_p[1] = t2;
    block_end_p[0] = t3;
  }
#else
  // Changing this version back to normal loads.
  //NO_INLINE
  static void FlushMusNoIl(LdlMuT *ldl_mu_p, LdlMuT * ldl_mu_end_p,
			   MuT* block_end_p, int n_mu) {
    // This approach gets rid of the mucking around to exit loop with
    // a test for eq / neq.
    static_assert(UnwindN == 4, "UnwindN mismatch in Flush1By1");
    // Make n_mu a multiple of UnwindN by moving any leftovers.
    if (n_mu & 2) {
      LdlMuT t4 = ldl_mu_end_p[-1];
      LdlMuT t5 = ldl_mu_end_p[-2];
      if (n_mu & 1) {  // 3 to do
	LdlMuT t6 = ldl_mu_end_p[-3];
	block_end_p[-1] = t4;
	block_end_p[-2] = t5;
	block_end_p[-3] = t6;
	block_end_p -= 3;
	ldl_mu_end_p -= 3;
      } else {  // only 2
	block_end_p[-1] = t4;
	block_end_p[-2] = t5;
	block_end_p -= 2;
	ldl_mu_end_p -= 2;
      }
    } else if (n_mu & 1) {  // 1 more to go
      LdlMuT t6 = ldl_mu_end_p[-1];
      block_end_p[-1] = t6;
      block_end_p--;
      ldl_mu_end_p--;
    }
    if (ldl_mu_p == ldl_mu_end_p) return;
    //FlushMusByMuNoIl(ldl_mu_p, ldl_mu_end_p, block_end_p);
    GdbHook(1);
    // N MUs to move is x * UnwindN.
    ldl_mu_end_p -= UnwindN;
    LdlMuT t0 = ldl_mu_end_p[+3];
    LdlMuT t1 = ldl_mu_end_p[+2];
    LdlMuT t2 = ldl_mu_end_p[+1];
    LdlMuT t3 = ldl_mu_end_p[+0];
    while (ldl_mu_p != ldl_mu_end_p) {
      ldl_mu_end_p -= UnwindN;
      block_end_p -= UnwindN;
      block_end_p[3] = t0;
      block_end_p[2] = t1;
      t0 = ldl_mu_end_p[+3];
      t1 = ldl_mu_end_p[+2];
      block_end_p[1] = t2;
      block_end_p[0] = t3;
      t2 = ldl_mu_end_p[+1];
      t3 = ldl_mu_end_p[+0];
    }
    block_end_p -= UnwindN;
    block_end_p[3] = t0;
    block_end_p[2] = t1;
    block_end_p[1] = t2;
    block_end_p[0] = t3;
  }

#endif  
  void FlushMusByMuNoIl(LdlMuT *ldl_mu_p, LdlMuT *ldl_mu_end_p,
			MuT* block_end_p) {
    GdbHook(1);
    ldl_mu_end_p -= UnwindN;
    LdlMuT t0 = ldl_mu_end_p[3];
    LdlMuT t1 = ldl_mu_end_p[2];
    LdlMuT t2 = ldl_mu_end_p[1];
    LdlMuT t3 = ldl_mu_end_p[0];
    while (ldl_mu_p != ldl_mu_end_p) {
      ldl_mu_end_p -= UnwindN;
      block_end_p -= UnwindN;
      block_end_p[3] = t0;
      t0 = ldl_mu_end_p[3];
      block_end_p[2] = t1;
      t1 = ldl_mu_end_p[2];
      block_end_p[1] = t2;
      t2 = ldl_mu_end_p[1];
      block_end_p[0] = t3;
      t3 = ldl_mu_end_p[0];
    }
    block_end_p -= UnwindN;
    block_end_p[3] = t0;
    block_end_p[2] = t1;
    block_end_p[1] = t2;
    block_end_p[0] = t3;
    
  }

// Version for interleaved ladle (IlN > 1). Must gather recs into simd unit.
  /**
   * Basically, we must gather the il data into a simd register and
   * then NT store the register to remote memory. The best way to gather
   * the data is not clear, and may vary by architecture. The docs
   * indicate that the intrinsics generate sequences of instrs and are
   * not necessarily very good!
   *
   * So I provide 4 ways to do the load:
   * - use an insert intrinsic to insert each item into place in the simd
   * - use a set intrinsic to fill the simd register with 1 "instr"
   * - gather the IL data to an aligned buffer and load the simd from there
   * - use the simd gather instr to fill the simd register.
   *
   * One of these is selected by the LOAD_SIMD_REG = {1, 2, 3, 4} parameter.
   **/
  
  NO_INLINE
  void FlushMusIl(IlElemT *ldl_p, IlElemT * ldl_end_p, MuT* block_end_p,
		  int n_e) {
    GdbHook(2);
    FlushMusByMuIl(ldl_p + RtPerSimd * (UnwindN-1), ldl_end_p, block_end_p);
    // Does this unwind?  Yes.
    for (int i = -1; i >= -UnwindN; i--) {
      LdlMuT t;
      ldl_end_p -= RtPerSimd;
      if (ldl_p > ldl_end_p) return;
      LoadSimd(ldl_end_p, t);
      block_end_p[i] = t;
    }
  }

  NO_INLINE void FlushMusByMuIl(IlElemT *ldl_p, IlElemT *&ldl_end_p,
                        MuT* &block_end_p) {
    GdbHook(1);
    while (ldl_p < ldl_end_p) {
      for (int i = -1; i >= -UnwindN; i--) {
	LdlMuT t;
	ldl_end_p -= RtPerSimd;
	LoadSimd(ldl_end_p, t);
	block_end_p[i] = t;
      }
      block_end_p -= UnwindN;
    }
  }

  // fill dest_reg with ldl_p[0, 1, ..., RtPerSimd-1].
  void LoadSimd(const IlElemT *ldl_p, LocalSdu &dest_reg) {
    if (LoadSimdReg == 1) {  // Use insert instrs
      // pretty good results (my I5 desktop).
      dest_reg. template Insert<0>(ldl_p[0]); 
      dest_reg. template Insert<1>(ldl_p[1]); 
      dest_reg. template Insert<2>(ldl_p[2]); 
      dest_reg. template Insert<3>(ldl_p[3]);
#if 0
      // Invalid calls to simd because byte 32 is invalid in 32 byte
      // register.  Compiler throws an error on the template number,
      // even though the if is compile time.
      if (RtPerSimd == 4) return;
      dest_reg. template Insert<4>(ldl_p[4]); 
      dest_reg. template Insert<5>(ldl_p[5]); 
      dest_reg. template Insert<6>(ldl_p[6]); 
      dest_reg. template Insert<7>(ldl_p[7]);
#endif
      return;
    } else if (LoadSimdReg == 2) {  // Use set instr
      // slightly better results than == 1 above
      if (RtPerSimd == 4)
	dest_reg. template Set(ldl_p[3], ldl_p[2], ldl_p[1], ldl_p[0]);
#if 0
      // Cause even though not callable, g++ looks for it:-(
      else if (RtPerSimd == 8) 
	dest_reg. template Set(ldl_p[7], ldl_p[6], ldl_p[5], ldl_p[4],
			       ldl_p[3], ldl_p[2], ldl_p[1], ldl_p[0]);
#endif
      return;
    } else if (LoadSimdReg == 3) {  // gather to Bfr
      // Awful results (worse than no simd).
      for (int i = 0; i < RtPerSimd; i++) Rec[i] = ldl_p[i];
      dest_reg = Bfr;
    } else if (LoadSimdReg == 4) {  // Use simd gather instr.
      IntelCLL* ldlP = reinterpret_cast<IntelCLL*>(ldl_p);
      dest_reg.MoveIn
	(_mm256_i32gather_epi64(ldlP, GatherNdx.SdU, sizeof(int64_t)));
    } else ErrTag("Bad LoadSimdReg config: die", true);
  }
// FlushMusIlSimdGthr                                    FlushMusIlSimdGthr

  NO_INLINE void FlushMusIlSimdGthr
  (IlElemT *ldl_p, IlElemT * ldl_end_p, MuT* block_end_p, int n_e) {
    GdbHook(3);
    GthrNdxSduT gNdx = GatherNdx.SdU;
    GatherByMus(ldl_p + RtPerSimd * (UnwindN-1), ldl_end_p, block_end_p, gNdx);
    for (int i = -1; i >= -UnwindN; i--) {
      ldl_end_p -= RtPerSimd;
      if (ldl_p > ldl_end_p) return;
      IntelCLL* ldlEndP = reinterpret_cast<IntelCLL*>(ldl_end_p);
      block_end_p[i].MoveIn
	(_mm256_i32gather_epi64(ldlEndP, gNdx, sizeof(int64_t)));
    }
  }  // FlushMusIlSimdGthr
  // GatherByMus                                                  GatherByMus
  /**
   * Gather just sucks. It might improve some when GatherByMus() is
   * inlined, but seems not likely to beat the simd set
   * sequence. Intel has designed the instruction to use mask, and
   * that overhead is there all the time, including creating a mask
   * copy and copying it fresh to a new register for each gather.
   *
   * This version seems to be about the same as copying the values 1
   * by 1 to a buffer and loading from there. Part of the cost will be
   * the need to prefetch the ladle space (not being done, but wasn't
   * done for the set instruction either).
   **/
  
  NO_INLINE
  void GatherByMus(IlElemT *ldl_p, IlElemT* &ldl_end_p, MuT* &block_end_p,
		   GthrNdxSduT g_ndx) {
    GdbHook(4);
    while (ldl_p < ldl_end_p) {
      for (int i = -1; i >= -UnwindN; i--) {
	ldl_end_p -= RtPerSimd;
	IntelCLL* ldlEndP = reinterpret_cast<IntelCLL*>(ldl_end_p);
	block_end_p[i].MoveIn
        (_mm256_i32gather_epi64(ldlEndP, g_ndx, sizeof(int64_t)));
      }
      block_end_p -= UnwindN;
    }

}  // GatherByMus

// Bfr is used to gather il data for loading to a simd register. 
union {
    LocalSdu Bfr;
    SRec_T Rec[RtPerSimd];
  } __attribute__ ((aligned (1 << MuT::SimdLogLen)));

  struct {
    union {
      GthrNdxSduT SdU;  // parameter to gather intrinsic
      GthrNdxT IV[IntPerGthrNdx] = { 0*IlN, 1*IlN, 2*IlN, 3*IlN};
    };
  } GatherNdx;
  
#else
  static const char *WhoAmI() { return "FlushSimd(V2:PartN)"; }

  // - INcrease SimdUnwindN to 8.

  // - Use (unwound) loop for partial blocks: Looks like it isn't
  //   the compares that are expensive, but breaking the pipe (Duh!). All
  //   the things I tried to not use a loop with lots of compares
  //   were not helpful.
  NO_INLINE void FlushMus(IlElemT *ldl_p, IlElemT * ldl_end_p,
			  MuT* block_end_p, int n_e) {

    int nMu = n_e >> LogRtPerMu;
    constexpr int SimdUnwindN = 8;  // N of simd ops to unwind
    // Move by blocks of SimdUnwindN explicitly unwound.
    // First, move any partial blocks
    int partN = nMu & (SimdUnwindN - 1);
    nMu -= partN;
    LdlMuT* ldlEndP = reinterpret_cast<LdlMuT*>(ldl_end_p) - partN;
      
    block_end_p -= partN;
    GdbHook(1);
    static_assert(SimdUnwindN == 8, "SimdUnwindN mismatch");

    for (int i = 0; i < partN; i++) block_end_p[i].NtMoveIn(ldlEndP[i]);

    // Now complete blocks
    while (nMu > 0) {
      // NOTE: code nice and compact, but always uses ymm0, so a
      // block stalls cpu.
      block_end_p[-1].NtMoveIn(ldlEndP[-1]);
      block_end_p[-2].NtMoveIn(ldlEndP[-2]);
      block_end_p[-3].NtMoveIn(ldlEndP[-3]);
      block_end_p[-4].NtMoveIn(ldlEndP[-4]);
      block_end_p[-5].NtMoveIn(ldlEndP[-5]);
      block_end_p[-6].NtMoveIn(ldlEndP[-6]);
      block_end_p[-7].NtMoveIn(ldlEndP[-7]);
      block_end_p[-8].NtMoveIn(ldlEndP[-8]);
      ldlEndP -= SimdUnwindN;
      block_end_p -= SimdUnwindN;
      nMu -= SimdUnwindN;
    }
  }  // FlushMus
#endif
} __attribute__ ((aligned (1 << MuT::SimdLogLen)));  // MoverClass

#endif

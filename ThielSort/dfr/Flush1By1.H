/* $Header:  $ */ 

#ifndef Flush1By1H
#define Flush1By1H

#include "PfSupport.H"

/********************************************************************
 *
 *                            Flush1By1.H
 *
 * Move data 1 element at a time. No simd, no non temporal instructions.
 * See FlushMovers.H for specs.
 ********************************************************************/

// Flush1By1                                                       Flush1By1
template<typename SRec_T>
class MoverClass<SRec_T, FLUSH_1BY1_MOVER> {  // Flush1By1
  static const int NativeRtPerMu = 1;

public:
  typedef SRec_T MuT;  // MU is just a record
  typedef IlElem<SRec_T, IlN> IlElemT;
    
  static constexpr int RtPerMu = FindRtPerMu(NativeRtPerMu);
  static constexpr int LogRtPerMu = FindLogRtPerMu(1, RtPerMu);
  static constexpr int UnwindN = 4;
  static constexpr int UWN = (RtPerMu > 1) ? RtPerMu : UnwindN;
  
  static const char *WhoAmI() {
    // Cause its a convenient place to put it.
    static_assert(RtPerMu > 0, "Config error with RtPerMu");
    const char* VersionId[] =
      {"Flush1By1(Unwind4 within loop)",
       "Flush1By1(Unwind4 with read ahead(n_mu test))"};

    return VersionId[FLUSH_MOVER_VARIANT];
  }
  
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  void FlushMus(IlElemT *ldl_p, IlElemT * ldl_end_p, MuT* block_end_p,
		int n_e) {
    /**
     * g++ doesn't like testing addresses for < or for >. It wants to
     * test for = or !=, and mucks around a lot to make sure that it
     * gets correct code. It looks like we have to make sure that we
     * have N%UnwindN = 0 before calling FlushMusByMu() to avoid a lot
     * of extra code. So if RtPerMu is < UnwindN, we move any n =
     * N%UnwindN here, before the call.
     *
     * NOTE: Even when I request a test for != in the FlushMusByMu() loop,
     *       g++ doesn't trust me and mucks around with the address to
     *       (presumably) make sure its right!
     *
     * NOTE: Assumes that if RtPerMu > UnwindN, it is a multiple. 
     **/
    if (RtPerMu < UnwindN) {
      static_assert(UnwindN == 4, "UnwindN mismatch in Flush1By1");
      // Move first few til n_e remaining is 0 mod UnwindN.
      if (int nTail = n_e % UnwindN) {
	if (nTail & 2) {
	  block_end_p[-1] = ldl_end_p[-1];
	  block_end_p[-2] = ldl_end_p[-2];
	  ldl_end_p -= 2;
	  block_end_p -= 2;
	}
	if (nTail & 1) {
	  block_end_p[-1] = ldl_end_p[-1];
	  block_end_p--;
	  ldl_end_p--;
	}
	n_e -= nTail;
      }
    }      
    // N to move is multiple of UnwindN.
    if (__builtin_expect(n_e > 0, true)) {
      FlushMusByMu(ldl_end_p, block_end_p, n_e / UnwindN);
    }
    PrefetchVec<0,3>(ldl_end_p - n_e, ldl_end_p);
#if 0
      // Simple loop
      while (ldl_end_p >= ldl_p) *--block_end_p = *--ldl_end_p;
#endif
  }

  private:
  // n_e % UnwindN == 0, n_e > 0.
  NO_INLINE void FlushMusByMu(IlElemT *ldl_end_p, MuT* block_end_p,
			      int n_mu) {
    GdbHook(1);
#if (0)
    while (ldl_p < ldl_end_p) {
      for (int i = -1; i >= -UWN; i--) block_end_p[i] = ldl_end_p[i];
      ldl_end_p -= UWN;
      block_end_p -= UWN;
    }
#else
    // We need N > 0, N to move is multiple of UnwindN, MuT is 1 record.
    // Read ahead.
    MuT t1 = ldl_end_p[-1];
    MuT t2 = ldl_end_p[-2];
    MuT t3 = ldl_end_p[-3];
    MuT t4 = ldl_end_p[-4];
    
    for (n_mu--; n_mu; n_mu--) {
      ldl_end_p -= UnwindN;
      block_end_p[-1] = t1; t1 = ldl_end_p[-1];
      block_end_p[-2] = t2; t2 = ldl_end_p[-2];
      block_end_p[-3] = t3; t3 = ldl_end_p[-3];
      block_end_p[-4] = t4; t4 = ldl_end_p[-4];
      block_end_p -= UnwindN;
    }
    // Store items that were read ahead.
    block_end_p[-1] = t1; 
    block_end_p[-2] = t2; 
    block_end_p[-3] = t3; 
    block_end_p[-4] = t4;
#endif
  }
  // Original (shotgun) version
  NO_INLINE void FlushMus1By1(IlElemT *ldl_p, IlElemT * ldl_end_p,
			      MuT* block_end_p, int n_e) {
    /////////// Changing SimdUnwindN may require fix below!!
    constexpr int SimdUnwindN = 8;  // N of moves to unwind
    // Move by blocks of SimdUnwindN explicitly unwound.
    // First, move any partial blocks
    int partN = n_e & (SimdUnwindN - 1);
    GdbHook(1);
    //// Assumes SimdUnwindN == 8!!!
    constexpr int move4 = 4;
    if (partN & move4) {
      ldl_end_p -= move4;
      block_end_p -= move4;
      for (int i = 0; i < move4; i++) block_end_p[i] = ldl_end_p[i];
    }
    constexpr int move2 = 2;
    if (partN & move2) {
      ldl_end_p -= move2;
      block_end_p -= move2;
      for (int i = 0; i < move2; i++) block_end_p[i] = ldl_end_p[i];
    }
    if (partN & 1) {
      ldl_end_p--;
      block_end_p--;
      block_end_p[0] = ldl_end_p[0];
    }
    // Now complete blocks
    while (ldl_p < ldl_end_p) {
      block_end_p[-1] = ldl_end_p[-1];
      block_end_p[-2] = ldl_end_p[-2];
      block_end_p[-3] = ldl_end_p[-3];
      block_end_p[-4] = ldl_end_p[-4];
      block_end_p[-5] = ldl_end_p[-5];
      block_end_p[-6] = ldl_end_p[-6];
      block_end_p[-7] = ldl_end_p[-7];
      block_end_p[-8] = ldl_end_p[-8];
      ldl_end_p -= SimdUnwindN;
      block_end_p -= SimdUnwindN;
    }
#if 0
    for (int i = 1; i <= n_e; i++) {
      if (ldlEndP[i] == 0xbdf1b1deb8b4a) ErrTag("Found suspect");
    }
#endif
  }

};  // Flush1By1

#endif

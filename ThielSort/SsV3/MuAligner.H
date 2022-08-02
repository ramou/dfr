/* $Header:  $ */ 

#ifndef MuAlignerH
#define MuAlignerH

/********************************************************************
 *
 *                            MuAligner.H
 *
 * Round addrs up or down so they are Mu aligned. Assumes:
 * - an MU is an integral number of Rec_Ts.
 * - an MU occupies an integral number of CLs.
 * - if an MU starts on a CL boundry, no other Rec_T in the MU will 
 *   be on a CL boundry.
 *
 * Testing for MU alignment is really testing for CL alignment.
 * If the size of the MU in bytes is a multiple of CL length, rounding 
 * an addr up or down to a CL boundry is fast.
 * Otherwise, it will require a small loop.
 ********************************************************************/

// MuAligner                                                       MuAligner
template<typename Rec_T, typename Mu_T> class MuAligner {
  using UUL = long long unsigned;
  union UUnion {
    UUnion(Rec_T *ptr) : p{ptr} {
      static_assert(sizeof(Rec_T) == sizeof(UUL), "sizeof(pointer is ??");
    }
    UUnion(Rec_T *, Mu_T *ptr) : muP{ptr} {
      static_assert(sizeof(Rec_T) == sizeof(UUL), "sizeof(pointer is ??");
    }

    bool Aligned() { return (i & (sizeof(Mu_T) - 1)) == 0; }

    Mu_T *muP;
    Rec_T *p;
    UUL i;
  };

public:
  constexpr static  int RtPerMu = sizeof(Mu_T) / sizeof(Rec_T);
  // IsPowerOf2 is true iff sizeof(Mu_T) is a power of 2.
  constexpr static bool IsPowerOf2 = (sizeof(Mu_T) & (sizeof(Mu_T) - 1)) == 0;

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  /**
   * Round p up to the next MU. We require that sizeof(MU) is a power of 2.
   **/
  static Rec_T *AlignMuTUp(Rec_T *p) {
    if (IsPowerOf2) {
      UUnion u(p);
      u.i = (u.i + (sizeof(Mu_T) - 1)) & ~(sizeof(Mu_T) - 1);
      return u.p;
    } else {
      while (!RtAligned(p)) p++;
      return p;
    }
  }
  /**
   * Round p down to the next MU. We require that sizeof(MU) is a power of 2.
   **/
  static Rec_T *AlignMuTDown(Rec_T *p) {
    if (IsPowerOf2) {
      UUnion u(p);
      u.i = u.i & ~(sizeof(Mu_T) - 1);
      return u.p;
    } else {
      while (!RtAligned(p)) p--;
      return p;
    }
  }
  /**
   * return true iff p is aligned on a CL boundry. If true, p can be
   * the addr of an aligned MU. They need different names when a MU is
   * just a RT.
   **/
  static bool RtAligned(Rec_T *p) {
    UUnion u(p);
    return u.Aligned();
  }
  static bool MuAligned(Mu_T *p) {
    UUnion u(nullptr, p);
    return u.Aligned();
  }
};  // MuAligner

#endif

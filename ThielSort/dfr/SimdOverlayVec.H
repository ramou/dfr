/* $Header:  $ */ 

#ifndef SimdOverlayVecH
#define SimdOverlayVecH

/********************************************************************
 *
 *                            SimdOverlayVec.H
 *
 * SimdOverlayVec overlays a simd vector and 2 base vectors of
 * types Base_T and Aux_T using an anonymous union. If we declare Name:
 *
 * SimdOverlayVec<...> Name;
 *
 * then Name.BV[], Name.XV[], and Name.SdV[] are the 3 vectors.
 *
 * template parameters:
 *
 * Base_T : the base (non-simd) type being manipulated
 *
 * Base_Vec_Len : the size of the Base_T vector (defaults to a vector
 *   long enough to be overlaid by a single simd value).  Note that if
 *   Base_Vec_Len isn't divisible by the simd length (NPerSimd) the
 *   last simd element will extend past the end of the Base_T vector.
 *
 * Log_Su_Len : the log2 of the size of the simd value (in bytes).
 *   5 => 32 byte values.  
 *
 * Aux_T : allows a Base_T to be treated as some other type using a
 *   union to overlay the bit pattern. Defaults to
 *   Base_T. sizeof(Base_T) == sizeof(Aux_T) is required.
 *
 * When Base_T is a pointer it can be useful to do arithmetic on the
 * elements with or without the simd instrs. Aux_T supports this sort
 * of thing by essentially providing a reinterpret_cast<Aux_T>().
 *
 * Vector sizes and types are provided. This is a struct, so all
 * values are public. 
 *
 * In principle, there can be no default initialization of unions,
 * since the value types cannot be assumed to be compatable. C++ thus
 * forbids a default constructor when a type contains an anonymous
 * union. A constructor with parameters must be used. We provide 2
 * constructors:
 *
 * - SimdOverlayVec() : does nothing. Space is not initialized.
 * - SimdOverlayVec(Base_T val) : inits each BV element to val.
 *
 ********************************************************************/

/**
 ****** I haven't built the infrastructure for AVX2 yet!!
 ****** Not even for AVX with Base_T = long!
 * 
 * So I will fool the compiler by making SimdT be just a wrapper for
 * Base_T. The code will work as expected, but won't use simd.
 **/
template<typename Base_T, int Base_Vec_Len = 0> struct FakeSimd {
  using BaseT = Base_T;
  static const int NPerSimd = 1;
  //FakeSimd(Base_T v) { Val = v; }
  FakeSimd& operator=(Base_T v) { Val = v;  return *this; }
  operator Base_T() const { return Val; }
    
  Base_T Val;
};

// SimdOverlayVec                                            SimdOverlayVec

template<typename Base_T, int Base_Vec_Len = 0,
	 int Log_Su_Len = SIMD_SHIFT_SIZE,
	 typename Aux_T = Base_T> 
struct SimdOverlayVec {
  using BaseT = Base_T;
  using AuxT = Aux_T;
  using SimdT = FakeSimd<Base_T, Base_Vec_Len>;
  // Should be SimdSupport::VectorType<Base_T, Log_Su_Len>;
  static const int NPerSimd = SimdT::NPerSimd;
  // size of base vector
  static const int BSize = Base_Vec_Len ? Base_Vec_Len : NPerSimd;
  // And size of simd vector.
  static const int SdSize = (Base_Vec_Len + NPerSimd - 1) / NPerSimd;  
  union {
    Base_T BV[BSize];  // Base Vector
    Aux_T XV[BSize];
    /////////////////////  No real simd yet!!!!!
    SimdT SdV[SdSize];  // Simd vector
  };
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  SimdOverlayVec() {
    static_assert(sizeof(Base_T) == sizeof(Aux_T),
		  "Type mismatch in SimdOverlayVec");
  }
  SimdOverlayVec(Base_T val) {
    for (int i = 0; i < BSize; i++) BV[i] = val;
  }
};  // SimdOverlayVec

#endif

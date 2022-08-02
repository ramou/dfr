/* $Header:  $ */ 

#ifndef VTBaseH
#define VTBaseH

/********************************************************************
 *
 *                            VTBase.H
 *
 * Base class for the vector types. Declares types and fields and
 * handles "conversions" between vectors of the same length.
 ********************************************************************/

namespace SimdSupport {

  
// VTBase                                                       VTBase
template<int Log_Size>
  class VTBase : private SpecificVecTypes<Log_Size> {
  typedef SpecificVecTypes<Log_Size> SpecTypes;
  
public:

  static const int SuLen = 1<<Log_Size;  // size of simd unit in bytes
  
  // Following are for use when the vector type matters to the user
  typedef typename SpecTypes::IVecType IVecType;
  typedef typename SpecTypes::UIVecType UIVecType;
  typedef typename SpecTypes::LVecType LVecType;
  typedef typename SpecTypes::ULVecType ULVecType;
  typedef typename SpecTypes::SVecType SVecType;
  typedef typename SpecTypes::USVecType USVecType;
  typedef typename SpecTypes::BVecType BVecType;
  typedef typename SpecTypes::UBVecType UBVecType;
  
  protected:
  
  typedef float SIMD_VEC(SuLen) FSVecType;  // for internal use at this time

#if 0
  // hope i don't need these!
  // Following are for use when the vector type matters to the intrinsics
  typedef __m128i AnonIntVec;  // Anonymous int data.
  typedef __m128d DoubleVec;   // Double precision
  typedef __m128  FloatVec;    // Float
#endif
  
  union {
      IVecType IVec;
      UIVecType UIVec;
      LVecType LVec;
      ULVecType ULVec;
      SVecType SVec;
      USVecType USVec;
      BVecType BVec;
      UBVecType UBVec;
      FSVecType FSVec;
    };
    // Store from VTBase
    inline void Store(const VTBase&src) { LVec = src.LVec; }
    // Store from individual vec types
    inline void Store(const IVecType &src) { IVec = src; }
    inline void Store(const UIVecType &src) { UIVec = src; }
    inline void Store(const LVecType &src) { LVec = src; }
    inline void Store(const ULVecType &src) { ULVec = src; }
    inline void Store(const SVecType &src) { SVec = src; }
    inline void Store(const USVecType &src) { USVec = src; }
    inline void Store(const BVecType &src) { BVec = src; }
    inline void Store(const UBVecType &src) { UBVec = src; }
    inline void Store(const FSVecType &src) { FSVec = src; }
    
  public:
    VTBase() : LVec{0,0} {}
    // Copy constructors
    VTBase(IVecType v) : IVec(v) {}
    VTBase(UIVecType v) : UIVec(v) {}
    VTBase(LVecType v) : LVec(v) {}
    VTBase(ULVecType v) : ULVec(v) {}
    VTBase(SVecType v) : SVec(v) {}
    VTBase(USVecType v) : USVec(v) {}
    VTBase(BVecType v) : BVec(v) {}
    VTBase(UBVecType v) : UBVec(v) {}
    VTBase(FSVecType v) : FSVec(v) {}
    // Casting operators
    operator IVecType() const { return IVec; }
    INLINE_ATT operator const UIVecType&() const { return UIVec; }
    operator LVecType() const { return LVec; }
    operator ULVecType() const { return ULVec; }
    operator SVecType() const { return SVec; }
    operator USVecType() const { return USVec; }
    operator BVecType() const { return BVec; }
    operator UBVecType() const { return UBVec; }
    operator FSVecType() const { return FSVec; }

  } SIMD_ALIGN(1<<Log_Size);  // VTBase

};  // namespace SimdSupport

#endif

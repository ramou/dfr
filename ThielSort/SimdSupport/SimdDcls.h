/* $Header:  $ */ 

#ifndef SimdDclsH
#define SimdDclsH

/********************************************************************
 *
 *                            SimdDcls.H
 *
 * A bunch of declarations that the implementation needs. Only the
 * first set should be needed to use the package.
 ********************************************************************/

namespace SimdSupport {

  // Names used to create app types

  // Select simd unit lengths. Values is log base 2 of the length in bytes:
  static const int LogLen16 = 4;
  static const int LogLen32 = 5;
  //static const int LogLen64 = 6;  // Not here yet

  // Select vector element type:
  typedef int SimdInt;
  typedef unsigned int SimdUInt;
  typedef long long int SimdLong;
  typedef long long unsigned int SimdULong;
  typedef short SimdShort;
  typedef unsigned short SimdUShort;
  typedef signed char SimdByte;
  typedef unsigned char SimdUByte;

  // Template used to create app types

  /**
   * Use this template to create and application type.
   *
   * Log_Su_Len must be one of the simd unit length tags above.
   * Elem_T must be one of the vector element types above
   *
   * template<typename Elem_T, int Log_Su_Len> class VectorType;
   **/  
  
  // Following shouldn't be used directly.

  // The vector types specific to the simd unit size. 
  template<int Log_Su_Len> class SpecificVecTypes;
  
  // SuLen == 16
  template<>class SpecificVecTypes<LogLen16> {
  public:
    typedef SimdInt    SIMD_VEC(16) IVecType;
    typedef SimdUInt   SIMD_VEC(16) UIVecType;
    typedef SimdLong   SIMD_VEC(16) LVecType;
    typedef SimdULong  SIMD_VEC(16) ULVecType;
    typedef SimdShort  SIMD_VEC(16) SVecType;
    typedef SimdUShort SIMD_VEC(16) USVecType;
    typedef SimdByte   SIMD_VEC(16) BVecType;
    typedef SimdUByte  SIMD_VEC(16) UBVecType;
  };
  // SuLen == 32
  template<>class SpecificVecTypes<LogLen32> {
  public:
    typedef SimdInt    SIMD_VEC(32) IVecType;
    typedef SimdUInt   SIMD_VEC(32) UIVecType;
    typedef SimdLong   SIMD_VEC(32) LVecType;
    typedef SimdULong  SIMD_VEC(32) ULVecType;
    typedef SimdShort  SIMD_VEC(32) SVecType;
    typedef SimdUShort SIMD_VEC(32) USVecType;
    typedef SimdByte   SIMD_VEC(32) BVecType;
    typedef SimdUByte  SIMD_VEC(32) UBVecType;
  };
};  // namespace SimdSupport

#endif

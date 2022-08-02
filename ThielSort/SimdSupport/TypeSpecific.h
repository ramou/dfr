/* $Header:  $ */ 

#ifndef TypeSpecificH
#define TypeSpecificH

/********************************************************************
 *
 *                            TypeSpecific.h
 *
 * Provides a templated class declaration to contain the code and
 * declarations specific to each specialization of VectorType<>. This
 * generic template is not defined, but specializations for each vec
 * length / element type provide the required implementations and
 * ensure type compatibility between the generically typed intrinsics
 * and the app specific types.
 *
 * Note that these classes are not intended to be instantiated even
 * though they are base classes of VectorType<>. All values and
 * functions are static. I intend them all to be identical except for
 * the template parameters that define the specialization and the
 * actual calls to the intrinsics.
 ********************************************************************/
#include <immintrin.h>

namespace SimdSupport {

  /**
   * Primary template for classes that implement the intrinsics. It is
   * not instantiated itself, but the various specializations are
   * implemented below. Template parameters:
   *
   * Elem_T : type of vector element, from the list of vector element
   *          types in SimdDcls.H.
   * Log_Su_Len : indicates size of simd unit, from the list of 
   *              of unit lengths in SimdDcls.H.
   **/
  template<class Elem_T, int Log_Su_Len> class ElemBase;

#if 0
  // useful??
  typedef VTBase<4> VTBase16;
  typedef VTBase<5> VTBase32;
#endif
  
// SimdInt                                                      SimdInt
  template<>class ElemBase<SimdInt, LogLen16> {
  protected:
    // Rename template parameters. Flag as such.
    static const int Log_Len = LogLen16;
    typedef SimdInt E_Type;

    static const int SuLen = 1 << Log_Len;
    typedef VTBase<Log_Len> AnonTyp;
    typedef E_Type SIMD_VEC(SuLen) VecType;

      // E_Type oriented
    inline static AnonTyp BroadCast(E_Type v)
    { return  _mm_set1_epi32(v); }
    inline static E_Type Extract(VecType &vec, int ndx)
    { return _mm_extract_epi32(vec, ndx);}
    inline static VecType Insert(VecType vec, int ndx, E_Type v) {
      switch (ndx) {
      case 0: return reinterpret_cast<VecType>(_mm_insert_epi32(vec, v, 0));
      case 1: return reinterpret_cast<VecType>(_mm_insert_epi32(vec, v, 1));
      case 2: return reinterpret_cast<VecType>(_mm_insert_epi32(vec, v, 2));
      case 3: return reinterpret_cast<VecType>(_mm_insert_epi32(vec, v, 3));
      }
      return vec;  // Never gets here! But the compiler needs it:-((
    }
    // Shifts and logical ops
    inline static VecType ShiftLeft(const VecType &vec, int bit_cnt)
    { return _mm_slli_epi32(vec, bit_cnt);}
    inline static VecType ShiftRight(const VecType &vec, int bit_cnt)
    { return _mm_srai_epi32(vec, bit_cnt);}
    inline VecType AndInstr(const VecType op1, const VecType op2)
    { return _mm_and_si128(op1, op2); }
    inline VecType AndNotInstr(const VecType op1, const VecType op2)
    { return _mm_andnot_si128(op1, op2); }
    inline VecType OrInstr(const VecType op1, const VecType op2)
    { return _mm_or_si128(op1, op2) ; }
    inline VecType XorInstr(const VecType op1, const VecType op2)
    { return _mm_xor_si128(op1, op2) ; }
    // Arithmetic
    inline static VecType IVMult
    ( VecType left_op, VecType right_op)
    { return _mm_mullo_epi32(left_op, right_op);}
    inline static VecType IVAdd
    (const VecType &left_op, const VecType right_op)
    { return _mm_add_epi32(left_op, right_op);}
    inline static VecType IVSub
    (const VecType left_op, const VecType right_op)
    { return _mm_sub_epi32(left_op, right_op);}
    // Miscellanious
    static void DumpVec(AnonTyp vec, const char *title) {
      VecType v = vec;
      const int nElem = sizeof(VecType) / sizeof(E_Type);
      std::cout << title << " = ";
      for (int i = 0; i < nElem - 1; i++) {
	std::cout << v[i] << ", ";
      }
      std::cout << v[nElem - 1] << "\n";
    }
  };  // ElemBase<SimdInt, LogLen16>

// SimdUInt                                                      SimdUInt
  template<>class ElemBase<SimdUInt, LogLen16> {
  protected:
    // Rename specialization arguments to template parameters. 
    static const int Log_Len = LogLen16;
    typedef SimdUInt E_Type;

    static const int SuLen = 1 << Log_Len;
    typedef VTBase<Log_Len> AnonTyp;
    typedef E_Type SIMD_VEC(SuLen) VecType;

    // E_Type oriented
    inline static AnonTyp BroadCast(E_Type v)
    { return  _mm_set1_epi32(v); }
    inline static E_Type Extract(const VecType &vec, int ndx)
    { return _mm_extract_epi32(vec, ndx);}
    inline static VecType Insert(VecType vec, int ndx, E_Type v) {
      switch (ndx) {
      case 0: return _mm_insert_epi32(vec, v, 0); 
      case 1: return _mm_insert_epi32(vec, v, 1); 
      case 2: return _mm_insert_epi32(vec, v, 2); 
      case 3: return _mm_insert_epi32(vec, v, 3); 
      }
      return vec;  // Never gets here! But the compiler needs it:-((
    }
    // Shifts and logical ops
    inline static VecType ShiftLeft(const VecType &vec, int bit_cnt)
    { return _mm_slli_epi32(vec, bit_cnt);}
    inline static VecType ShiftRight(const VecType vec, int bit_cnt)
    { return _mm_srli_epi32(vec, bit_cnt);}
    inline VecType AndInstr(const VecType op1, const VecType op2)   
    { return _mm_and_si128(op1, op2); }
    inline VecType AndNotInstr(const VecType op1, const VecType op2)   
    { return _mm_andnot_si128(op1, op2); }
    inline VecType OrInstr(const VecType op1, const VecType op2)   
    { return _mm_or_si128(op1, op2) ; }
    inline VecType XorInstr(const VecType op1, const VecType op2)
    { return _mm_xor_si128(op1, op2) ; }
    // Arithmetic
    inline static VecType IVMult
      (const VecType left_op, const VecType right_op)
    { return _mm_mullo_epi32(left_op, right_op);}
    inline static VecType IVAdd
      (const VecType &left_op, const VecType &right_op)
    { return _mm_add_epi32(left_op, right_op);}
    inline static VecType IVSub
      (const VecType left_op, const VecType right_op)
    { return _mm_sub_epi32(left_op, right_op);}
    // Miscellanious
    static void DumpVec(const VecType v, const char *title) {
      const int nElem = sizeof(VecType) / sizeof(E_Type);
      std::cout << title << " = ";
      for (int i = 0; i < nElem - 1; i++) {
	std::cout << v[i] << ", ";
      }
      std::cout << v[nElem - 1] << "\n";
    }
  };  // ElemBase<SimdUInt, LogLen16>

#if 0  
// SimdLong                                                      SimdLong
  // SLongVec                                                       SLongVec
  template<>class ElemBase<SimdLong, LogLen16> SLongVec {

  public:

    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    SLongVec();
    ~SLongVec();

    //////////////////////////////////////////////////////////////////////
    //                  Readers / Writers.

    //}Readers Tag for placing autogened readers.
    //////////////////////////////////////////////////////////////////////
    //                   Other functions.

  } SIMD_ALIGN(1<<LogSuLength);  // SLongVec
// SimdULong                                                      SimdULong
  template<>class ElemBase<SimdULong, LogLen16> ULongVec {

  public:

    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    ULongVec();
    ~ULongVec();

    //////////////////////////////////////////////////////////////////////
    //                  Readers / Writers.

    //}Readers Tag for placing autogened readers.
    //////////////////////////////////////////////////////////////////////
    //                   Other functions.

  } SIMD_ALIGN(1<<LogSuLength);  // ULongVec
// SimdShort                                                      SimdShort
  template<>class ElemBase<SimdShort, LogSuLength> SShortVec {

  public:
    typedef short E_Type;
    typedef SVecType NativeType;
    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    SShortVec();
    ~SShortVec();

    //////////////////////////////////////////////////////////////////////
    //                   Other functions.
    inline VecType BroadCast(short v, VecType*)
    { return  _mm_set1_epi16(v); }
    inline int Extract(IVecType vec, int ndx)
    { return _mm_extract_epi16(vec, ndx);}
    inline SVecType Insert(IVecType vec, int ndx, short v)
    { _mm_insert_epi16(vec, v, ndx); return reinterpret_cast<SVecType>(vec);}
    inline __m128i ShiftLeft(SVecType vec, int bit_cnt)
    { return _mm_slli_epi16(vec, bit_cnt);}
    inline __m128i ShiftRight(SVecType vec, int bit_cnt)
    { return _mm_srai_epi16(vec, bit_cnt);}

  } SIMD_ALIGN(1<<LogSuLength);  // SShortVec
// SimdUShort                                                      SimdUShort
  template<>class ElemBase<SimdUShort, LogLen16> UShortVec {

  public:

    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    UShortVec();
    ~UShortVec();

    //////////////////////////////////////////////////////////////////////
    //                   Other functions.
    inline VecType BroadCast(UShort v, VecType*)
    { return  _mm_set1_epi16(v); }
    inline int Extract(IVecType vec, int ndx, UShort)
    { return _mm_extract_epi16(vec, ndx);}
    inline USVecType Insert(IVecType vec, int ndx, UShort v)
    { _mm_insert_epi16(vec, v, ndx); return reinterpret_cast<USVecType>(vec);}
    inline __m128i ShiftLeft(USVecType vec, int bit_cnt)
    { return _mm_slli_epi16(vec, bit_cnt);}
    inline __m128i ShiftRight(USVecType vec, int bit_cnt)
    { return _mm_srli_epi16(vec, bit_cnt);}

  } SIMD_ALIGN(1<<LogSuLength);  // UShortVec
// SimdByte                                                      SimdByte
  template<>class ElemBase<SimdByte, LogLen16> SByteVec {

  public:

    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    SByteVec();
    ~SByteVec();

    //////////////////////////////////////////////////////////////////////
    //                  Readers / Writers.

    //}Readers Tag for placing autogened readers.
    //////////////////////////////////////////////////////////////////////
    //                   Other functions.

  } SIMD_ALIGN(1<<LogSuLength);  // SByteVec
// SimdUByte                                                      SimdUByte
  template<>class ElemBase<SimdUByte, LogLen16> UByteVec {

  public:

    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    UByteVec();
    ~UByteVec();

    //////////////////////////////////////////////////////////////////////
    //                  Readers / Writers.

    //}Readers Tag for placing autogened readers.
    //////////////////////////////////////////////////////////////////////
    //                   Other functions.

  } SIMD_ALIGN(1<<LogSuLength);  // UByteVec
#endif  

};  // namespace SimdSupport

#endif

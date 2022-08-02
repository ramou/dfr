/* $Header:  $ */ 

#ifndef SrtSimdBaseH
#define SrtSimdBaseH

/********************************************************************
 *
 *                            SrtSimdBase.H
 *
 * Contains SrtSimd base classes.  Allows SrtSimd types to come from a
 * single template, even for the 4- and 8-byte cases.
 ********************************************************************/

#include <immintrin.h>

namespace SimdUnitMover {

  /**
   * Implementation for different size SDUs. Unit_Log_Size is the log of
   * the number of bytes in the SDU (i.e., 2, 3, 4, 5, or 6). This
   * primary template has no definition. Specializations for each length
   * is given below. Note that Unit_Log_Size == 2 or 3 don't really
   * specify simd operands, so can only be used with assignment.
   **/
  template<int Unit_Log_Size> class SrtSimdBase;

// SrtSimdBase<2>                                              SrtSimdBase<2>
  template<> class SrtSimdBase<2> {
  public:
    typedef int IntelSdType;
    static constexpr int SimdLogLen = 2;
    static constexpr int SimdLen = 1 << SimdLogLen;
    
    // Move v to me.
    INLINE_ATT void NtMoveIn(const SrtSimdBase<2> v)
    { _mm_stream_si32(&SdVec, v.SdVec); }
    INLINE_ATT void MoveIn(const IntelSdType &v) { SdVec = v; }
    // Insert v into slot Ndx.
    template<int Ndx, class Rec_Type> void Insert(Rec_Type v)
    { MyAbort("No insert into i32 'simd' operands"); }
    template<class Rec_Type> void Set(Rec_Type e3, Rec_Type e2,
				      Rec_Type e1, Rec_Type e0) {
      MyAbort("No Set() instr for i32 'simd' operands");
    }
  
  protected:
    IntelSdType SdVec;
  } __attribute__ ((aligned (1 << 2)));  // SrtSimdBase<2>
  // SrtSimdBase<3>                                              SrtSimdBase<3>
  template<> class SrtSimdBase<3> {
  public:
    typedef long long int IntelSdType;
    static constexpr int SimdLogLen = 3;
    static constexpr int SimdLen = 1 << SimdLogLen;

    // Move v to me.
    INLINE_ATT void NtMoveIn(const SrtSimdBase<3> v)
    { _mm_stream_si64(&SdVec, v.SdVec); }
    INLINE_ATT void MoveIn(const IntelSdType &v) { SdVec = v; }
    operator IntelSdType() const { return SdVec; }
    // Insert v into slot Ndx.
    template<int Ndx, class Rec_Type> void Insert(Rec_Type v)
    { MyAbort("No insert into i64 'simd' operands"); }
    template<class Rec_Type> void Set(Rec_Type e3, Rec_Type e2,
				      Rec_Type e1, Rec_Type e0) {
      MyAbort("No Set() instr for i64 'simd' operands");
    }
  
  protected:
    IntelSdType SdVec;
  } __attribute__ ((aligned (1 << 3)));  // SrtSimdBase<3>
// SrtSimdBase<4>                                              SrtSimdBase<4>
  template<> class SrtSimdBase<4> {
  public:
    typedef __m128i IntelSdType;
    static constexpr int SimdLogLen = 4;
    static constexpr int SimdLen = 1 << SimdLogLen;

    // Move v to me.
    INLINE_ATT void NtMoveIn(const SrtSimdBase<SimdLogLen> &v)
    { _mm_stream_si128(&SdVec, v.SdVec); }
    INLINE_ATT void NtMoveIn(const IntelSdType &v)
    { _mm_stream_si128(&SdVec, v); }
    INLINE_ATT void MoveIn(const SrtSimdBase<SimdLogLen> &v)
    { SdVec = v.SdVec; }
    INLINE_ATT void MoveIn(const IntelSdType &v)
    { _mm_store_si128(&SdVec, v); }
    // Insert v into slot Ndx.
    template<int Ndx, class Rec_Type> void Insert(Rec_Type v) 
    { SdVec = _mm_insert_epi64(SdVec, v, Ndx); }
    template<class Rec_Type> void Set(Rec_Type e1, Rec_Type e0) {
      SdVec = _mm_set_epi64x(e1, e0);
    }

  protected:
    IntelSdType SdVec;
  } __attribute__ ((aligned (1 << 4)));  // SrtSimdBase<4>

// SrtSimdBase<5>                                              SrtSimdBase<5>
  template<> class SrtSimdBase<5> {
  public:
    typedef __m256i IntelSdType;
    static constexpr int SimdLogLen = 5;
    static constexpr int SimdLen = 1 << SimdLogLen;

    // Move v to me.
    INLINE_ATT void NtMoveIn(const SrtSimdBase<SimdLogLen> &v)
    { _mm256_stream_si256(&SdVec, v.SdVec); }
    INLINE_ATT void NtMoveIn(const IntelSdType &v)
    { _mm256_stream_si256(&SdVec, v); }
    INLINE_ATT void NtMoveIn(const SrtSimdBase<SimdLogLen>* p)
    { SdVec = _mm256_stream_load_si256(&(p->SdVec)); }
    INLINE_ATT void MoveIn(const SrtSimdBase<SimdLogLen> &v)
    { //_mm256_store_si256(&SdVec, v.SdVec);
      SdVec = v.SdVec;
    }
    INLINE_ATT void MoveIn(const IntelSdType &v)
    { _mm256_store_si256(&SdVec, v);
    }

    INLINE_ATT const IntelSdType &Get() const { return SdVec; }

    // Insert v into slot Ndx.
    template<int Ndx, class Rec_Type> void Insert(Rec_Type v)
    { SdVec = _mm256_insert_epi64(SdVec, v, Ndx); }
    template<class Rec_Type> void Set(Rec_Type e3, Rec_Type e2,
				      Rec_Type e1, Rec_Type e0) {
      SdVec = _mm256_set_epi64x(e3, e2, e1, e0);
    }
    
  protected:
    IntelSdType SdVec;
  } __attribute__ ((aligned (1 << SimdLogLen)));  // SrtSimd32Base

  using SrtSimd4Base = SrtSimdBase<2>;
  using SrtSimd8Base = SrtSimdBase<3>;
  using SrtSimd16Base = SrtSimdBase<4>;
  using SrtSimd32Base = SrtSimdBase<5>;
  //using SrtSimd64Base = SrtSimdBase<6>;

};  // SimdUnitMover

#endif

/* $Header:  $ */ 

#ifndef Simd16BaseH
#define Simd16BaseH

/********************************************************************
 *
 *                            Simd16Base.H
 *
 * Map my needs to intel intrinsics for the 128-bit simd case.
 * include IntrinUser.h, not this file!
 ********************************************************************/


// Simd16Base                                               Simd16Base
// Following are for use when the vector type matters
typedef int   SIMD_VEC(16) IVecType16;
typedef unsigned int   SIMD_VEC(16) UIVecType16;
typedef long long int  SIMD_VEC(16) LVecType16;
typedef unsigned long long int  SIMD_VEC(16) ULVecType16;
typedef short SIMD_VEC(16) SVecType16;
typedef unsigned short SIMD_VEC(16) USVecType16;
typedef SByte  SIMD_VEC(16) BVecType16;
typedef UByte SIMD_VEC(16) UBVecType16;
typedef float SIMD_VEC(16) FSVecType16;

// Following propagate the constant across the vector.
inline SdvType16 BroadCast(SByte v, SdvType16*) 
{ return  _mm_set1_epi8(v); }
inline SdvType16 BroadCast(UByte v, SdvType16*) 
{ return  _mm_set1_epi8(v); }
inline SdvType16 BroadCast(short v, SdvType16*) 
{ return  _mm_set1_epi16(v); }
inline SdvType16 BroadCast(UShort v, SdvType16*) 
{ return  _mm_set1_epi16(v); }
inline SdvType16 BroadCast(int v, SdvType16*) 
{ return  _mm_set1_epi32(v); }
inline SdvType16 BroadCast(UInt v, SdvType16*) 
{ return  _mm_set1_epi32(v); }
inline SdvType16 BroadCast(long v, SdvType16*) 
{ return  _mm_set1_epi64x(v); }
inline SdvType16 BroadCast(ULong v, SdvType16*) 
{ return  _mm_set1_epi64x(v); }

inline int GetSignBits(SdvType16 v, int word_len) {
  if (word_len == 4) {  // Get Signs from ints, etc.
    return _mm_movemask_ps(v);
  } else if (word_len == 1) {  // Get Signs from bytes.
    return _mm_movemask_epi8(v);
  }
  myabort("Invalid call to GetSignBits()");
  return 0;  // Never gets here.
}
  
// following are all combine 2 vecs of the same size and return a
// value of the same type.

 // Load / store instructions
 inline void StoreNtInstr(Vector16 *dest, const Vector16 src)  
 { _mm_stream_si128(dest, src); }
 inline SdvType16 LoadNtInstr (const Vector16 *src) 
 { return _mm_stream_load_si128((Vector16*)src); }
 inline SdvType16 LoadInstr (const Vector16 *p) 
 { return _mm_load_si128((const __m128i*)p); }
 inline void StoreInstr(Vector16 *dest, const SdvType16 src)  
 { _mm_store_si128(dest, src); }
// Memory to memory move
 inline void MoveNtInstr(Vector16 *dest, const Vector16 *src)
 { Vector16 t = _mm_load_si128(src); _mm_stream_si128(dest, t); }
 inline void NTCopy1CL(const Vector16 *src, Vector16 *dest) {
   for (int i = 0; i < 4; i++) StoreNtInstr(dest + i, src[i]);
 }

// Extract / Insert
// Note: 128 bit extracts seem to be unsigned moves. Not 256 bit versions!
// Note: 128 bit versions short seems to be much more expensive!?
//inline int Extract(IVecType16 vec, int ndx, int);
inline int Extract(IVecType16 vec, int ndx, int) {
  return _mm_extract_epi32(vec, ndx);
}
inline int Extract(IVecType16 vec, int ndx, UInt) {
  return _mm_extract_epi32(vec, ndx);
}
inline int Extract(IVecType16 vec, int ndx, short) {
  return _mm_extract_epi16(vec, ndx);
}
inline int Extract(IVecType16 vec, int ndx, UShort) {
  return _mm_extract_epi16(vec, ndx);
}
inline int Extract(IVecType16 vec, int ndx, SByte) {
  return _mm_extract_epi8(vec, ndx);
}
inline int Extract(IVecType16 vec, int ndx, UByte) {
  return _mm_extract_epi8(vec, ndx);
}
inline IVecType16 Insert(IVecType16 vec, int ndx, int v) {
  _mm_insert_epi32(vec, v, ndx);
  return vec;
}
inline UIVecType16 Insert(IVecType16 vec, int ndx, UInt v) {
  _mm_insert_epi32(vec, v, ndx);
  return reinterpret_cast<UIVecType16>(vec);
}
inline SVecType16 Insert(IVecType16 vec, int ndx, short v) {
  _mm_insert_epi16(vec, v, ndx);
  return reinterpret_cast<SVecType16>(vec);
}
inline USVecType16 Insert(IVecType16 vec, int ndx, UShort v) {
  _mm_insert_epi16(vec, v, ndx);
  return reinterpret_cast<USVecType16>(vec);
}
inline BVecType16 Insert(IVecType16 vec, int ndx, SByte v) {
  _mm_insert_epi8(vec, v, ndx);
  return reinterpret_cast<BVecType16>(vec);
}
inline UBVecType16 Insert(IVecType16 vec, int ndx, UByte v) {
  _mm_insert_epi8(vec, v, ndx);
  return reinterpret_cast<UBVecType16>(vec);
}
// Bit shift instructions. 

/////// NOTE: Can avoid passing DigitDef by template to get constant
/////// shift at cost of an xmm register to hold count. Not doing it now.
// Shift left. logical and arithmetic give same result.
inline __m128i ShiftLeft(SVecType16 vec, int bit_cnt) {
  return _mm_slli_epi16(vec, bit_cnt);
}
inline __m128i ShiftLeft(IVecType16 vec, int bit_cnt) {
  return _mm_slli_epi32(vec, bit_cnt);
}
inline __m128i ShiftLeft(LVecType16 vec, int bit_cnt) {
  return _mm_slli_epi64(vec, bit_cnt);
}
inline __m128i ShiftLeft(USVecType16 vec, int bit_cnt) {
  return _mm_slli_epi16(vec, bit_cnt);
}
inline __m128i ShiftLeft(UIVecType16 vec, int bit_cnt) {
  return _mm_slli_epi32(vec, bit_cnt);
}
inline __m128i ShiftLeft(ULVecType16 vec, int bit_cnt) {
  return _mm_slli_epi64(vec, bit_cnt);
}

// Shfit right
// Shift type (arithmetic or logical) is from parameter type.
inline __m128i ShiftRight(SVecType16 vec, int bit_cnt) {
  return _mm_srai_epi16(vec, bit_cnt);
}
inline __m128i ShiftRight(IVecType16 vec, int bit_cnt) {
  return _mm_srai_epi32(vec, bit_cnt);
}
#if 0
// Missing from docs, includes???
inline __m128i ShiftRight(LVecType16 vec, int bit_cnt) {
  return _mm_srai_epi64(vec, bit_cnt);
}
#endif
inline __m128i ShiftRight(USVecType16 vec, int bit_cnt) {
  return _mm_srli_epi16(vec, bit_cnt);
}
inline __m128i ShiftRight(UIVecType16 vec, int bit_cnt) {
  return _mm_srli_epi32(vec, bit_cnt);
}
inline __m128i ShiftRight(ULVecType16 vec, int bit_cnt) {
  return _mm_srli_epi64(vec, bit_cnt);
}

// Arithmetic instrs
inline __m128i IVMult(const IVecType16 op1, const IVecType16 op2) {
  return _mm_mullo_epi32(op1, op2);
}
template<class T>
inline __m128i IVAdd(const T op1, const T op2);
template<>
inline __m128i IVAdd(const IVecType16 op1, const IVecType16 op2) {
  return _mm_add_epi32(op1, op2);
}
template<>
inline __m128i IVAdd (const USVecType16 op1, const USVecType16 op2) {
  return _mm_adds_epu16(op1, op2);
}
#ifdef Simd16Base_home
  static const RegisterSource Simd16Base_H_Registry("$Header:  $");
#endif

#endif


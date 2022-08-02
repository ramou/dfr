inline SdvType16 BroadCast(int v, SdvType16*) { return  _mm_set1_epi32(v); }
inline SdvType16 BroadCast(UInt v, SdvType16*) { return  _mm_set1_epi32(v); }
inline int Extract(IVecType16 vec, int ndx, int) {return _mm_extract_epi32(vec, ndx);}
inline int Extract(IVecType16 vec, int ndx, UInt) {return _mm_extract_epi32(vec, ndx);}
inline IVecType16 Insert(IVecType16 vec, int ndx, int v) {_mm_insert_epi32(vec, v, ndx); return vec;}
inline UIVecType16 Insert(IVecType16 vec, int ndx, UInt v) {_mm_insert_epi32(vec, v, ndx); return reinterpret_cast<UIVecType16>(vec);}
inline __m128i ShiftLeft(IVecType16 vec, int bit_cnt) {return _mm_slli_epi32(vec, bit_cnt);}
inline __m128i ShiftLeft(UIVecType16 vec, int bit_cnt) {return _mm_slli_epi32(vec, bit_cnt);}
inline __m128i ShiftRight(IVecType16 vec, int bit_cnt) {return _mm_srai_epi32(vec, bit_cnt);}
inline __m128i ShiftRight(UIVecType16 vec, int bit_cnt) {return _mm_srli_epi32(vec, bit_cnt);}
inline __m128i IVMult(const IVecType16 op1, const IVecType16 op2) {return _mm_mullo_epi32(op1, op2);}
inline __m128i IVAdd(const IVecType16 op1, const IVecType16 op2) {return _mm_add_epi32(op1, op2);}

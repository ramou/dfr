inline SdvType16 BroadCast(short v, SdvType16*)
{ return  _mm_set1_epi16(v); }
inline SdvType16 BroadCast(UShort v, SdvType16*)
{ return  _mm_set1_epi16(v); }
inline int Extract(IVecType16 vec, int ndx, short)
{ return _mm_extract_epi16(vec, ndx);}
inline int Extract(IVecType16 vec, int ndx, UShort)
{ return _mm_extract_epi16(vec, ndx);}
inline SVecType16 Insert(IVecType16 vec, int ndx, short v)
{ _mm_insert_epi16(vec, v, ndx); return reinterpret_cast<SVecType16>(vec);}
inline USVecType16 Insert(IVecType16 vec, int ndx, UShort v)
{ _mm_insert_epi16(vec, v, ndx); return reinterpret_cast<USVecType16>(vec);}
inline __m128i ShiftLeft(SVecType16 vec, int bit_cnt)
{ return _mm_slli_epi16(vec, bit_cnt);}
inline __m128i ShiftLeft(USVecType16 vec, int bit_cnt)
{ return _mm_slli_epi16(vec, bit_cnt);}
inline __m128i ShiftRight(SVecType16 vec, int bit_cnt)
{ return _mm_srai_epi16(vec, bit_cnt);}
inline __m128i ShiftRight(USVecType16 vec, int bit_cnt)
{ return _mm_srli_epi16(vec, bit_cnt);}

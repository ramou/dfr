/* $Header:  $ */
/******************************************************************
 
                        IsConstBase.C

*******************************************************************/

template<int Inc_Val>
Elem_T *FindSebBndry_imp(Elem_T *src_p, Elem_T *src_e, Elem_T seb_tmv) {
  if (TopMask(src_e - Inc_Val) == seb_tmv) return src_e;  // All the same.

  int nearNdx = 0;
  int farNdx = (src_e - src_p);
  /**
   * src_p[nearNdx] is in the SeBlock and src_p[farNdx] is not. We
   * will move nearNdx towards farNdx, but make sure it always indexes
   * an item in the SeBlock. Similarly, we will move farNdx towards
   * nearNdx, but make sure it always indexes an element outside the
   * SeBlock. When they are adjacent, we have found the boundry of the
   * SeBlock.
   *
   * In both while loops we have to protect against an infinite loop
   * when nearNdx and farNdx are adjacent, since if we aren't careful
   * the "midpoint" of the two adjacent values can round to a value
   * that continues the loop. We deal with this by adding/subtracting
   * Inc_Val so that rounding is in the direction which will terminate
   * the loop. 
   **/

  int curNdx = farNdx / 2;
  while (abs(farNdx - nearNdx) > 20) {
    while (TopMask(src_p + curNdx) == seb_tmv) {
      nearNdx = curNdx;
      // Round in far direction to force loop exit when adjacent.
      //if (curNdx==((nearNdx + farNdx + Inc_Val) / 2)) ErrTag();
      curNdx = (nearNdx + farNdx + Inc_Val) / 2;
    }
    if (__builtin_expect(curNdx == farNdx, false)) {
      // Only when nearNdx and farNdx are adjacent.
#ifdef OCD_TESTING
      if ((nearNdx + Inc_Val) != farNdx) ErrTag("FindSebBndry1");
      if (TopMask(src_p + nearNdx) != seb_tmv) ErrTag("FindSebBndry11");
      if (TopMask(src_p + farNdx) == seb_tmv) ErrTag("FindSebBndry12");
#endif
      return src_p + curNdx;
    }
    farNdx = curNdx;
    curNdx = (nearNdx + farNdx - Inc_Val) / 2;
    while (TopMask(src_p + curNdx) != seb_tmv) {
      farNdx = curNdx;
      // Round in near direction to force loop exit when adjacent.
      //if (curNdx==((nearNdx + farNdx - Inc_Val) / 2)) ErrTag();
      curNdx = (nearNdx + farNdx - Inc_Val) / 2;
    }
    nearNdx = curNdx;
    curNdx = (nearNdx + farNdx + Inc_Val) / 2;
  }
  nearNdx += Inc_Val;
  while (seb_tmv == TopMask(src_p[nearNdx])) nearNdx += Inc_Val;
#ifdef OCD_TESTING
  if (seb_tmv == TopMask(src_p[nearNdx])) ErrTag("FindSebBndry2a");
  if (seb_tmv != TopMask(src_p[nearNdx - Inc_Val])) ErrTag("FindSebBndry2b");
#endif
  return src_p + nearNdx;
}

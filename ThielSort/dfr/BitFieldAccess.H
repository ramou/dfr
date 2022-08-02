/* $Header:  $ */ 

#ifndef BitFieldAccessH
#define BitFieldAccessH

/********************************************************************
 *
 *                            BitFieldAccess.H
 *
 * This file contains 2 class definitions:
 *
 * 1 - BFDefinition : A class to specify a set of contiguous bit
 *     fields. It is templated with a set of field widths. It contains
 *     static functions to return fields widths and offsets (both in
 *     bits) for a specific field and asserts to make sure that the
 *     requested field id is valid. 
 *
 * 2 - BitFieldAccess : A class to access information about a specific
 *     bit field from a BFDefinition class, and to access that field
 *     from a particular data item. It is templated with a
 *     BFDefinition class and a field id.
 *
 * Fields are numbered from 1 to NFlds. Field number 0 is defined in
 * BFDefinition to specify low order bits to be ignored, but field 0
 * cannot be accessed via BitFieldAccess.
 ********************************************************************/

// BFDefinition                                                 BFDefinition

/**
 * Define the field structure. Template parameters are Fld*_W, where *
 * is a Fld_Id, and the value is the width of field Fld*_Id. Fld0_W is
 * the bit offset of the low order bit of the low order field, thus
 * the width of the field before field 1.
 **/
template<int ...w>class BFDefinition {
 public:
  // Number of defined fields (not counting field 0).
  static constexpr int NFields = sizeof...(w) -1;
  // Range check ndx and return fields width.
  template<int ndx>
    static int constexpr Width() { return Widths[ndx]; }
  // Return the bit offset of field fld_id.
  template<int Fld_Id> static constexpr int Offset() {
    return (Fld_Id <= 1) ? Widths[0] :
      Offset<(Fld_Id <= 1) ? 1 : Fld_Id - 1>() + Widths[Fld_Id - 1];
  }

private:
  static const constexpr int Widths[sizeof...(w)] = {w...};
};
  
/*
 * BitFieldAccess<> contains static functions to return contiguous bit
 * field values for any defined field. Fields are assumed to
 * be within a long value, so they fit in a register and can be
 * accessed via shift and mask. 
 * NOTE: Valid fields numbered 1:Fld_Defs::NFields.
 ********************************************************************/

// BitFieldAccess                                    BitFieldAccess
template<int Fld_Id, typename Fld_Defs>
class BitFieldAccess : private Fld_Defs {
  static_assert((Fld_Id > 0) && (Fld_Id <= Fld_Defs::NFields),
		"BitFieldAccess::Fld_ID out of range");
public:
  static constexpr const int FldId = Fld_Id;
  // Field width
  static constexpr const int FW = Fld_Defs:: template Width<Fld_Id>();
  // Field offset 
  static constexpr const int FldOffset =
    Fld_Defs:: template Offset<Fld_Id>();
  
  // Return a mask for the number of low order bits in field Fld_Id.
  static constexpr int Mask() {
    return (1 << FW) - 1;
  }
  // Return a mask for the unshifted field.
  static unsigned long InPlaceMask() {
    return ((1L << FW) -1) << FldOffset;
  }
  static constexpr int NumBuckets() { return 1 << FW; }
  // Return the bits in field Fld_Id of v, shifted to the low order
  // position. 
  template<typename T> static int GetField(T v) {
    return (v >> FldOffset) & Mask();
  }
  // Load target byte (if feasible) or load *p and call the above.
  template<typename T> static int GetField(T *p) {
    if ((FW == 8) && ((FldOffset & 7) == 0)) {
      // return 8 bit byte
      return reinterpret_cast<const unsigned char*>(p)[FldOffset >> 3];
    } else {
      return GetField(*p);
    }
  }
  // Duplicate original calls, even though they refer to bytes, not
  // arbitrary fields. Easier than changing the old code.
  template<typename T>static int GetDealB(T v) { return GetField(v); }
  template<typename T>static int GetDealB(T *p) { return GetField(p); }
  // Not the way it should be done for flexibility, but usable.
  // Return the bits in field Fld_Id+1 so counting can be done for the
  // next pass.
  template<typename T> static int GetCountB(T v) {
    // if a next field is available, FIdInc == 1, else 0.
    constexpr const int FIdInc = (Fld_Id < Fld_Defs::NFields) ? 1 : 0;
    if (FIdInc) {
      typedef BitFieldAccess<Fld_Id+FIdInc, Fld_Defs> CountDigitT;
      return CountDigitT::GetField(v);
    } else {
      ErrTag("No count digit past last digit", true);
      return 0;  // Not real, since ErrTag() will die.
    }
  }
};  // BitFieldAccess

#endif

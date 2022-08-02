/* $Header:  $ */ 

#ifndef DealByteAccessH
#define DealByteAccessH

/********************************************************************
 *
 *                            DealByteAccess.H
 *
 * DealByteAccess<> contains static functions to efficiently access
 * specific bytes from an ELEM. These functions should be moved to a
 * smarter ELEM class.
 ********************************************************************/

// DealByteAccess                                            DealByteAccess
template<typename ELEM, int Deal_Byte, int Count_Byte = Deal_Byte+1>
class DealByteAccess {

 public:

  static const int DealByte = Deal_Byte;
  static const int CountByte = Count_Byte;

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  //DealByteAccess();
  //~DealByteAccess();

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  static int GetDealB(const ELEM *p) {
    return reinterpret_cast<const unsigned char*>(p)[Deal_Byte];
  }
  static int GetDealB(ELEM v) {
    static_assert(sizeof(ELEM) <= sizeof(long),
                  "ELEM won't fit in a register");
    return (unsigned char)
      (reinterpret_cast<unsigned long>(v) >> (Deal_Byte * 8));
  }
  static int GetCountB(const ELEM *p) {
    return reinterpret_cast<const unsigned char*>(p)[Count_Byte];
  }
  static int GetCountB(ELEM v) {
    static_assert(sizeof(ELEM) <= sizeof(long),
                  "ELEM won't fit in a register");
    if (CountByte >= sizeof(long)) {
      ErrTag("Calling for DBA::Count_Byte overflow");
      return 0;
    } else {
      return (unsigned char)
        (reinterpret_cast<unsigned long>(v) >>
         ((CountByte >= sizeof(long))) ? 0 : (Count_Byte * 8));
    }
  }

};  // DealByteAccess

#endif

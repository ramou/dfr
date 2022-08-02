/* $Header:  $ */ 

#ifndef BcSortDescH
#define BcSortDescH

/********************************************************************
 *
 *                            BcSortDesc.H
 *
 * Base class for sort descriptor (Sd) classes. Includes
 * infrastructure stuff common to various Sd classes, but not
 * necessarily all of them. No template parameters so members are
 * visible with no irritating additions. Currently, just diverting
 * stuff.
 *
 ********************************************************************/

// BcSortDesc                                                    BcSortDesc
class BcSortDesc {
protected:

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.
  
  ////////////////////////////////////////////////////////////////////////
  // Diversion support

  /**
   * Find the number of deals to use by looking through
   * DiversionThresholds[] and return it.
   **/
  static int FindNDeals(CountType n) {
    int nDeals = 1;
    for (int i = 0; i < DIVERSION_ARRAY_SIZE; i++) {
      if (n < DiversionThresholds[i]) return nDeals;
      if (DiversionThresholds[i] == 0) return nDeals;
      nDeals++;
    }
    MyAbort("Invalid DiversionThresholds vector");  // noreturn!
    return 0;  // so never executed.
  }
  static void PrintDiversionThresholds() {
    if (TimingOption <= 1) return;
    std::cout << "DiversionThresholds[] = {";
    for (int i = 0; i < DIVERSION_ARRAY_SIZE - 1; i++) {
      std::cout << DiversionThresholds[i] << ", ";
    }
    std::cout << DiversionThresholds[DIVERSION_ARRAY_SIZE - 1] << "}\n";
  }

};  // BcSortDesc

#endif

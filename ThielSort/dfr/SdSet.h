/* $Header:  $ */ 

#ifndef SdSeth
#define SdSeth

/********************************************************************
 *
 *                            SdSet.h
 *
 * Collect the available sort descriptor .H files here and select 1
 * according to the SELECT_SORT_DESCRIPTOR define constant.
 ********************************************************************/

// Generic sort descriptor class template. Doesn't exist.
template<typename SRec_T, int Selector_Id>class SdClass;

// Various specializations of SdClass:

#include "SdDivertingHybridDp.H"

#if 0
#include "SdDivertingDp.H"
#include "SdNarrowDp.H"
#include "SdNoLdlIlBkts.H"
#include "SdStdSort.H"
#include "SdSimpleRadix.H"
#include "SdDealerExecRadix.H"
#include "SdSplitting.H"
#include "SdSegmentedDFR.H"
#include "SdDivertAndConquor.H"
#include "SdMsdLadleless.H"
#endif

template<typename SRec_T>
using SdTemplate = SdClass<SRec_T, SELECT_SORT_DESCRIPTOR>;

#endif

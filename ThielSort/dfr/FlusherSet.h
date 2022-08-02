/* $Header:  $ */ 

#ifndef FlusherSeth
#define FlusherSeth

/********************************************************************
 *
 *                            FlusherSet.h
 *
 * Collect the available flusher .H files here and select 1 according
 * to the SELECT_FLUSHER_EXEC define constant.
 ********************************************************************/

// Generic flusher class template. Doesn't exist.
template<typename SRec_T, int Selector_Id>class FlusherClass;

// Various specializations of FlusherClass:

#if SELECT_FLUSHER_EXEC == DK_FLUSHER
#include "DirectKissFlusher.H"
#elif SELECT_FLUSHER_EXEC == DBIN_FLUSHER
#include "DirectBinFlusher.H"
#elif SELECT_FLUSHER_EXEC == DSPLIT_FLUSHER
#include "DirectSplitFlusher.H"
#elif SELECT_FLUSHER_EXEC == DLADLELESS_FLUSHER
#include "DirectLadlelessFlusher.H"
#elif SELECT_FLUSHER_EXEC == DHYBRID_FLUSHER
#include "HybridFlusher.H"
#elif SELECT_FLUSHER_EXEC == DK2_FLUSHER
#include "DK2Flusher.H"
#endif

template<typename SRec_T>
using FlusherTemplate = FlusherClass<SRec_T, SELECT_FLUSHER_EXEC>;

#endif

/* $Header:  $ */ 

#ifndef LegoDefsH
#define LegoDefsH

/********************************************************************
 *
 *                            LegoDefs.h
 *
 * There are a variety of experimental sort implementations which can
 * be created by combining modules which implement required actions
 * for the various steps (e.g., dealers, flushers, and mover modules
 * (e.g., simd or not simd)). This file briefly describes the various
 * actions required in (almost) all experiments, and the choices
 * available for each action. It defines module identifer macros
 * (MIDs) to identify each module and identifies (but doesn't define)
 * the SELECT macro which must be defined (with definitions on the
 * make command or in the various pattern definitions) to be one of
 * those MIDs.
 *
 * The MIDs are chosen so that they are unique across the set of
 * actions.  Not all combinations are legal, but we don't go into that
 * here. All the modules are templates specialized by the MID. The
 * general template is not defined.
 *
 ********************************************************************/

/**
 * DealerCore classes handle the deal of a single item. They will also
 * do bucket counting, the livebits analysis, and ladle prefetches if
 * needed. A DealerCore module must be assigned to macro
 * SELECT_DEALER_CORE.
 **/
// Selects direct store ladles, no prefetch.
#define DIRECT_STORE_DEALER_CORE 101
// Selects direct store ladles with ladle PF.
#define DIRECT_STORE_DEALER_CORE_WITH_PF 102
// Selects direct store bypassing ladles, storing directly to buckets.
#define DIRECT_STORE_DEALER_NO_LADLE 103
// Selects direct store with slicing dealer core
#define SLICING_DEALER_CORE 104
// Selects direct store with pf for current store
#define DCPF_REG_QUEUE_DEALER_CORE 105

/**
 * Flusher macros describe the flusher executive. One must be assigned
 * to SELECT_FLUSHER_EXEC. (The 'D' prefix is Direct.)
 **/
// Select DirectKissFlusher
#define DK_FLUSHER 201
// Select DirectBinFlusher
#define DBIN_FLUSHER 202
#if 0
// Select Leapfrog flusher
#define DLEAPFROG_FLUSHER 203
#endif
#define DSPLIT_FLUSHER 204
// Select hybrid flusher, where most dealing bypasses the ladle.
#define DHYBRID_FLUSHER 205
// Select really simple flusher modeled on hybrid flusher.
#define DK2_FLUSHER 206
// Next version
#define DK3_FLUSHER 207

/**
 * Movers move data from the ladle to the buckets under control of
 * FMoverT. One must be assigned to SELECT_FLUSH_MOVER.
 **/
// Select Flush1By1.
#define FLUSH_1BY1_MOVER 301
// Select FlushParaCL.
#define FLUSH_PARA_CL_MOVER 302
// Select FlushSimd.
#define FLUSH_SIMD_MOVER 303
// Select 8 byte NT store
#define FLUSH_NT8_MOVER 304
// Select FlushBinPair (2 ladles in a bin in parallel)
#define FLUSH_BIN_PAIR 305
// Select FastSimd mover.
#define FAST_SIMD_MOVER 306

/**
 * DataPool classes. Two versions of leapfrog, neither of which is as
 * fast as I expected. May be a way, but ...
 * One must be assigned to SELECT_DATAPOOL.
 *
 * NOTE: Leapfrog not supported at this time.
 **/
// Select DataPool.
#define CHUNK_DATAPOOL 401
#if 0
// Select version 1 leapfrog
#define LEADFROG_VN_1_DATAPOOL 402
// Select version 2 leapfrog
#define LEADFROG_VN_2_DATAPOOL 403
#endif

/**
 * Sort descriptor classes. The somewhat missnamed SdNarrowDp.H does
 * all 8 passes, but the wide / narrow choice is made separately.
 * One must be assigned to SELECT_SORT_DESCRIPTOR. Many are out of date!
 **/
// Select the somewhat missnamed SdNarrowDp.H: all 8 passes with
// ladles, but the wide / narrow choice is made separately.
#define FR_WITH_LADLES_SD 501
// Select diverting fast radix with ladles.
#define FR_DIVERTING_WITH_LADLES_SD 502
// Select sort descriptor for no diverting and no ladles.
#define FR_DIRECT_TO_BUCKETS 503
// Select sort descriptor for std::sort()
#define STD_SORT 504
// Select sort descriptor for simple imp of simple radix sort
#define SIMPLE_RADIX_SIMPLE_IMP 505
// Select sort descriptor for simple radix with a counting
// pass. Differs from above by using DealerExec to deal. No ladles.
#define SIMPLE_RADIX_DEALEREXEC_IMP 506
// Select sort descriptor to deal MSD pass to split data before using
// DFR on rest.
#define MSD_TO_SPLIT_THEN_DFR 507
// Select sort descriptor to sort a large vector in sections, then
// merge the sorted sections.
#define SEGMENTED_LARGE_VECTOR 508
// Selects the Sd that will segment large vectors using an initial MSD
// deal pass.
#define FR_DIVERTING_WITH_MSD_SEGMENTS_SD 509
// Selects the Sd using MSD segmentation with the ladleless deal.
#define FR_MSD_LADLELESS_SD 510
// Selects the Sd that uses the hybrid flusher. Initial changes from
// FR_DIVERTING_WITH_LADLES_SD look to be small, but a distinct SD
// seems cleaner.
#define FR_DIVERTING_WITH_HYBRID_FLUSHER_SD 510

#endif

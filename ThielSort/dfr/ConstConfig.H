/* $Header:  $ */ 

#ifndef ConstConfigH
#define ConstConfigH

#include "LegoDefs.h"

/********************************************************************
 *
 *                            ConstConfig.H
 *
 * This file declares namespace ConstConfig and populates it with
 * static const members to control the behaviour of the sort. Config.H
 * will add a sort descriptor (SD) class in the namespace based on the
 * CONFIG_PATTERN value. CONFIG_PATTERN must be set to one of the
 * available configuration pattern ids (SD_*) defined below via a
 * Extras="-DSD_*" option on the make command that compiles the sort.
 *
 * Note that the exploration patterns that were never intended to be
 * actual sorts have been removed from here and from Config.H. See
 * ../CurDfr for them. I have now removed the leapfrog options from
 * this version.
 *
 * See Config.H for more information.
 *
 * Structure of this file.
 *
 * 0 - Macro definitions used to identify each SD, and a note on the
 *     algorithms and modules that SD will use.  Only one SD of each
 *     type will be visible in a compile file.
 *
 * 1 - A "#if" nest to select custom defaults based on the
 *     CONFIG_PATTERN. In some cases defaults will be forced, but
 *     mostly they can be overridden by makefile defines (via
 *     Extras="..."). Values without custom defaults are ignored here.
 *
 * 2 - An unconditional block to supply defaults for all values which
 *     haven't been given defaults yet. This section lists all options
 *     and provides some documentation for them. Some options are
 *     irrelevant in some SDs, and some get specific values in some
 *     SDs.
 *
 * 3 - A set of static const definitions which declare all the
 *     constants. These members are global to the sort, ensuring that
 *     all parts are using consistent values even though a wide
 *     variety of customizations are available.
 *
 * 4 - A PrintConfig() function that is called to print the effective
 *     values of all the config parameters. It is a no-op unless
 *     DUMP_CONFIG is true.
 *
 * NOTE: Most of these defines are not used in the actual code. They
 * are used to define static const constants (later in this include
 * file). The constant's names are formed from the define names by
 * keeping the first letter and the first letter after an '_' as
 * capitals, changing all other letters to lower case, and removing
 * all '_' chars. Thus THIS_NAME is used to initialize the static
 * const constant ThisName. This doesn't apply to name defined in
 * LegoDefs.h (which are used in typedef and using statements) or to
 * the SD names defined immediately below (which are used to specify
 * lego names and other define values (in this file), but are not used
 * directly).
 ********************************************************************/

////////////// Configuration patterns : Choose 1 for CONFIG_PATTERN
////////////// Default = none (Compile fails)
#define SD_DIRECT_1BY1_PATTERN 1
#define SD_DIRECT_NTSIMD_PATTERN 2
#define SD_DIRECT_PARACL_PATTERN 5
// Use dealer core that can prefetch deal target after each store. The
// prefetch request is a template parameter to Deal1(), so we don't
// need to pf the same way each time.
#define SD_DIRECT_WITH_PF_1BY1_PATTERN 6
// As above, with simd flush.
#define SD_DIRECT_WITH_PF_NTSIMD_PATTERN 7
// As SD_DIRECT_WITH_PF_1BY1_PATTERN, but deal directly to bucket.
#define SD_DIRECT_TO_BUCKET_WITH_PF_PATTERN 8
// As SD_DIRECT_WITH_PF_1BY1_PATTERN, but use FlushParaCL
#define SD_DIRECT_WITH_PF_PARACL_PATTERN 9

// Leapfrog options deleted from this version

#define SD_DIRECT_DIVERTING_WITH_PF_1BY1_PATTERN 16
////#define SD_DRCT_NOLDL_FULL_SORT_PATTERN 17
#define SD_STD_SORT_PATTERN 18
#define SD_SIMPLE_SIMPLE_RADIX_PATTERN 19
#define SD_SIMPLE_DEALEREXEC_RADIX_PATTERN 20
#define SD_DIRECT_DIVERTING_WITH_PF_NTSIMD_PATTERN 21
#define SD_MSD_TO_SPLIT_THEN_DFR 22
#define SD_SEGMENTED_LARGE_VECTOR 23

// Start of ladle bin versions. All are DIRect DIVerting LadleBin
// (SD_DIR_DIV_LB ...)

#define SD_DIR_DIV_DBIN_WITH_PF_1BY1_PATTERN 24
#define SD_DIR_DIV_DBIN_WITH_PF_NTSIMD_PATTERN 25

// Start of DirectSplitFlusher versions.

#define SD_DIR_DIV_SPLITF_WITH_PF_1BY1_PATTERN 26
#define SD_DIR_DIV_SPLITF_WITH_PF_NTSIMD_PATTERN 27

// versions with NT8 mover
#define SD_DIR_DIV_SPLITF_WITH_PF_NT8_PATTERN 28
#define SD_DIR_DIV_LB_WITH_PF_NT8_PATTERN 29
#define SD_DIR_DIV_DBIN_WITH_PF_NT8_PATTERN 30

// Version with FlushBinPair mover
#define SD_DIR_DIV_DBIN_WITH_PF_BIN_PAIR_PATTERN 31

// Versions with SLICING_DEALER_CORE
#define SD_DIRECT_DIVERTING_WITH_SLICING_NTSIMD_PATTERN 32

// Versions with DCPF_REG_QUEUE_DEALER_CORE
#define SD_DIRECT_DIVERTING_WITH_DCPFCUR_NTSIMD_PATTERN 33

// versions with FR_DIVERTING_WITH_LADLES_SD
#define SD_DIRECT_WITH_SEGMENTATION_ON_MSD_PATTERN 34

//versions with MSD segmentation and DirectLadlelessFlusher.
#define SD_DIRECT_MSD_SEG_LADLELESS_PATTERN 35

//versions with hybrid flusher (deal direct to target bkts until there
//isn't room to hold worst case, then switch to ladle for that bkt).
#define SD_DIR_DIV_HYBRID_WITH_PF_1BY1_PATTERN 36
#define SD_DIR_DIV_HYBRID_WITH_PF_NTSIMD_PATTERN 37

#ifndef DUMP_CONFIG
#define DUMP_CONFIG true
#endif

// SystemPageSize controls memory page prefetch for bucket space.
// == -1 means suppress attempts to prefetch page, else system page size.
// Haven't made this work well for me. For now, default is -1.
#ifndef SYSTEM_PAGE_SIZE
#define SYSTEM_PAGE_SIZE -1
#endif

////////////// Configurations

/**
 * Note that the dealer cost can be influenced by the choice of
 * flusher! The combination is what must be tuned, not the individual
 * parts, at least when using ntsimd or non ntsimd flushers, due the
 * the cache pollution introduced if the buckets are being brought
 * into cache.
 **/

#if CONFIG_PATTERN == SD_DIRECT_1BY1_PATTERN
/////////////////////////////// Override defaults as required
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_NTSIMD_PATTERN
/////////////////////////////// Override defaults as required
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_WITH_LADLES_SD
// Don't prefetch target of NT store!
#define BKT_N_TO_PRIME -1
// Don't (for now, at least) try to load SIMD ops from IL ladle
#define IL_VALUE 1

#elif CONFIG_PATTERN == SD_DIRECT_PARACL_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_PARA_CL_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_WITH_PF_1BY1_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_WITH_PF_NTSIMD_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_WITH_LADLES_SD
// Don't prefetch target of NT store!
#define BKT_N_TO_PRIME -1
// Don't (for now, at least) try to load SIMD ops from IL ladle
#define IL_VALUE 1
// LDL_N_TO_PRIME applies to flush. Leave for now.
//#define LDL_N_TO_PRIME -1

#elif CONFIG_PATTERN == SD_DIRECT_WITH_PF_PARACL_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_PARA_CL_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_WITH_LADLES_SD

// Leapfrog options deleted from this version

#elif CONFIG_PATTERN ==SD_DIRECT_DIVERTING_WITH_PF_1BY1_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_DIVERTING_WITH_PF_NTSIMD_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#ifndef SELECT_FLUSHER_EXEC
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#endif
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_DIVERTING_WITH_SLICING_NTSIMD_PATTERN
#define SELECT_DEALER_CORE SLICING_DEALER_CORE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_DIVERTING_WITH_DCPFCUR_NTSIMD_PATTERN
#define SELECT_DEALER_CORE DCPF_REG_QUEUE_DEALER_CORE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIRECT_WITH_SEGMENTATION_ON_MSD_PATTERN
#define SELECT_DEALER_CORE DCPF_REG_QUEUE_DEALER_CORE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
// NOTE: Can't use Simd mover, since segment vectors aren't aligned:-((
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_MSD_SEGMENTS_SD

#elif CONFIG_PATTERN == SD_DIRECT_MSD_SEG_LADLELESS_PATTERN
#define SELECT_DEALER_CORE DCPF_REG_QUEUE_DEALER_CORE
#define SELECT_FLUSHER_EXEC DLADLELESS_FLUSHER
// NOTE: no mover used!
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_MSD_LADLELESS_SD

#elif CONFIG_PATTERN ==SD_STD_SORT_PATTERN
#define SELECT_SORT_DESCRIPTOR STD_SORT
// Rest because the system expects selectors to be defined
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_NO_LADLE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#elif CONFIG_PATTERN ==SD_SIMPLE_SIMPLE_RADIX_PATTERN
#define SELECT_SORT_DESCRIPTOR SIMPLE_RADIX_SIMPLE_IMP
// Rest because the system expects selectors to be defined
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_NO_LADLE
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL

#elif CONFIG_PATTERN == SD_MSD_TO_SPLIT_THEN_DFR
#define SELECT_SORT_DESCRIPTOR MSD_TO_SPLIT_THEN_DFR
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL

#elif CONFIG_PATTERN == SD_SEGMENTED_LARGE_VECTOR
#define SELECT_SORT_DESCRIPTOR SEGMENTED_LARGE_VECTOR
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL

#elif CONFIG_PATTERN == SD_DIR_DIV_DBIN_WITH_PF_1BY1_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DBIN_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_DBIN_WITH_PF_NTSIMD_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DBIN_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_SPLITF_WITH_PF_1BY1_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DSPLIT_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_SPLITF_WITH_PF_NT8_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DSPLIT_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_NT8_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_LB_WITH_PF_NT8_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DK_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_NT8_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_DBIN_WITH_PF_NT8_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DBIN_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_DBIN_WITH_PF_BIN_PAIR_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DBIN_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_BIN_PAIR
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_LADLES_SD

#elif CONFIG_PATTERN ==SD_DIR_DIV_HYBRID_WITH_PF_1BY1_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#define SELECT_FLUSHER_EXEC DHYBRID_FLUSHER
#define SELECT_FLUSH_MOVER FLUSH_1BY1_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_HYBRID_FLUSHER_SD

#elif CONFIG_PATTERN == SD_DIR_DIV_HYBRID_WITH_PF_NTSIMD_PATTERN
#define SELECT_DEALER_CORE DIRECT_STORE_DEALER_CORE_WITH_PF
#ifndef SELECT_FLUSHER_EXEC
#define SELECT_FLUSHER_EXEC DHYBRID_FLUSHER
#undef WIDE_INDICES
#define WIDE_INDICES true
#endif
#define SELECT_FLUSH_MOVER FLUSH_SIMD_MOVER
#define SELECT_DATAPOOL CHUNK_DATAPOOL
#define SELECT_SORT_DESCRIPTOR FR_DIVERTING_WITH_HYBRID_FLUSHER_SD
#undef IL_VALUE
#define IL_VALUE 1

#else
// Cause abort at preprocesser level to avoid tons of meaningless errors.
#include "WeNeedAConfigPattern.H"

#endif

// OptionList                                                   OptionList
// Apply global defaults as required
#ifndef COUNT_TYPE
// Must use long long if sort N doesn't fit in int.
#define COUNT_TYPE int
#endif
#ifndef BATCH_SIZE
#define BATCH_SIZE 15000
#endif
#ifndef USE_EFFECTIVE_BATCH_SIZE
#define USE_EFFECTIVE_BATCH_SIZE true
#endif
#ifndef BUCKET_EXPANSION_MULT
#define BUCKET_EXPANSION_MULT 1.0
#endif
/**
 * A ladle / bin is flushable if TRIGGER_VALUE < the current ladle /
 * bin len. Since flushing is done in terms of MUs, any leftovers
 * won't really be flushed, so the given value will be rounded down to
 * a multiple of MUs.
 **/
#ifndef TRIGGER_VALUE
#define TRIGGER_VALUE 48
#endif
#ifndef IL_VALUE
#define IL_VALUE 1
#endif
#ifndef Max_N_BUCKETS
#define Max_N_BUCKETS 256
#endif
// If 1, pads structs to force alignment of fields in classes. Not simd!
#ifndef FORCE_WORD_ALIGNMENT
#define FORCE_WORD_ALIGNMENT 0
#endif
// Control version of load to simd reg from IL source: 1 - 6.
// All bad!! None of these seem viable.
#ifndef LOAD_SIMD_REG
#define LOAD_SIMD_REG 1
#endif
// Select code variants within a flush mover.
#ifndef FLUSH_MOVER_VARIANT
#define FLUSH_MOVER_VARIANT 1
#endif
// Control version of FlushBinPair.H flush mover.
// 1 - Unwind / latency is 2 items deep.
// 2 - Unwind / latency is 4 items deep.
// No longer relevant.
#ifndef FLUSH_BIN_PAIR_VN
#define FLUSH_BIN_PAIR_VN 1
#endif
/**
 * Control prefetch of source during the deal. We prefetch (locality
 * 0) a CL this many items in advance of dealing them.
 **/
#ifndef READ_PF_ADVANCE
// Following was near optimal on home desktop.
//#define READ_PF_ADVANCE 152
#define READ_PF_ADVANCE 240
#endif

/**
 * Force flush movers to move multiples of their native RtPerMu. 
 * Default is let native stand.
 **/
#ifndef EXTENDED_RT_PER_MU
#define EXTENDED_RT_PER_MU -1
#endif

#ifndef WIDE_INDICES
#define WIDE_INDICES false
#endif

/**
 * Using the ladle index before incing it sometimes seems to force an
 * extra mov instr in the inner deal loop. If DC_POST_INC == 1, we use
 * it then inc it. If == 0, we inc first. Used in some DealerCores.
 **/
#ifndef DC_POST_INC
#define DC_POST_INC 1
#endif

// Insertion sort parameters
/**
 * If using diverting FR and the number of passes is odd, we can
 * combine the copy of the DataPool to the final result with the
 * insertion sort (USE_COPY_INSERTION_SORT == 1) or do the copy and
 * then do the insertion sort as a separate pass
 * (USE_COPY_INSERTION_SORT == 0). Default is 1.
 **/
#ifndef USE_COPY_INSERTION_SORT
#define USE_COPY_INSERTION_SORT 0
#endif
/**
 * Select size of "short run" to just sift current item down and
 * assume no other sifting is needed in same SeBlock. At 2, the
 * assumption is always true. If > 2, we may need 1 (or more) sifts in
 * same SeBlock, but still may pay by saving when next few are in
 * order even if SeBlock isn't ended.
 **/
#ifndef ISNS_SHORT_RUN
#define ISNS_SHORT_RUN 2
#endif
/**
 * Doing the IS in a bunch of little pieces adds a bit of overhead,
 * but doing a copy instead of the copyIS is more expensive, I
 * think. If nPasses is even, we might save a bit delaying the IS til
 * the entire segment is done. If it's odd, doing the copy back to the
 * src vector probably is costly. Default to doing it on each segment.
 **/
#ifndef DO_IS_ON_DIP_SEGMENTS
#define DO_IS_ON_DIP_SEGMENTS 1
#endif
#if 0
// Not useful as tested. Maybe in some other context???
/**
 * We use different functions for the insertion sort depending on
 * whether the Average Count Per Bucket (or per SeBlock) (Apcb) is
 * expected to be small or large. If Apcb < ACPB_THRESHOLD, we use
 * DoBracketedSort_imp(). Otherwise LRBracketedSort_imp() (Long Run).
 * Initial value (0.85) is a blind guess! Set to 0.0 or 100.0 to
 * always use the LR or not LR versions, resp.
 **/
#ifndef ACPB_THRESHOLD
#define ACPB_THRESHOLD 0
#endif

/**
 * What size of scan batch should we use with LRBracketedSort()
 * (Insertion sort step)?
 **/
#ifndef LR_SCAN_SIZE
#define LR_SCAN_SIZE 2048
#endif
#endif

/**
 * Control prefetching of target bucket space during flush. 
 *
 * If the flush is using NT instructions, target bucket prefetching
 * must be turned off by setting BKT_N_TO_PRIME to -1, since it has no
 * value but still causes the cache to be polluted and uses up memory
 * bandwidth.
 *
 * We need enough prefetch lead time to allow the target CLs to be
 * fetched before we store to them. We address this by finding the
 * list of ladles to be flushed and prefetching the targets for the
 * first BktNToPrime flushable ladles. The main loop prefetches the
 * bucket space for a ladle and flushes a ladle. (A prefetch precedes
 * a flush in the loop to permit BktNToPrime = 0 to do the prefetch
 * immediatly before the flush.) This doesn't give an exact lead time
 * since the size of the flushable ladles will vary, but it is a
 * reasonable approximation, and much simpler. Two values control this
 * mechanism:
 *
 * 1 - BKT_N_TO_PRIME : the number of buckets prefetched in
 * advance. If BKT_N_TO_PRIME is 0, the bucket is prefetched
 * immediately before it is filled. If 1, 1 bucket is prefetched after
 * bucket b is prefetched and before bucket b is filled, etc. If
 * BKT_N_TO_PRIME is < 0, bucket prefetching is suppressed completely.
 *
 * 2 - BKT_LOCALITY : the locality hint for the prefetch. Zero tells
 * the cpu that we will use the space only briefly, and it should
 * avoid polluting the lower level cache with it. For buckets, this is
 * probably right, but since we expect to move several data items into
 * each CL, we don't want them removed immediately after the first
 * access, so perhaps 1 is better. Also, if BKT_N_TO_PRIME > 0, we
 * need the CL to hang around until we get there.
 **/
#ifndef BKT_N_TO_PRIME
#define BKT_N_TO_PRIME -1
//#define BKT_N_TO_PRIME 10
#endif
#ifndef BKT_LOCALITY
#define BKT_LOCALITY 0
#endif
/**
 * Control prefetching of ladle source during flush.
 *
 * Similar to above, except that ladle space is expected to be in some
 * level of cache most of the time, since it is reused each batch,
 * while bucket space (at the item level, at least) is used only once
 * per pass, so must normally be prefetched from main memory. Thus we
 * use locality = 3 to try to keep the active parts of the ladle
 * somewhere in cache, and this isn't configurable. Also, since on
 * average the prefetch comes from a "closer" location, the ladle
 * fetch lead time (LDL_N_TO_PRIME) is expected to be <=
 * BKT_N_TO_PRIME. The current code requires (and enforces) this.
 **/
#ifndef LDL_N_TO_PRIME
#define LDL_N_TO_PRIME -1
//#define LDL_N_TO_PRIME 4
#endif

/**
 * Control prefetching of ladle source during deal.
 *
 * 1 - control how often we do the ladle prefetch. If we do it too
 *     often we waste time prefetching CLs multiple times (may not
 *     have much cost). If we don't do it often enough, we may miss
 *     doing it soon enough to be useful (significant cost, although
 *     most ladle CLs will probably be somewhere in the cache
 *     system). The IL_VALUE will affect this tradeoff.
 *
 * 2 - control how far ahead of current store point we prefetch.
 *
 **/
/**
 * Do prefetch every nth deal. Currently 0-2, where 0 suppresses the
 * prefetch.
 **/
#ifndef DEALER_PF_FREQ
#define DEALER_PF_FREQ 1
#endif
/**
 * Specify prefetch lead when pf is done in DealerExec. Number of
 * elements in same bucket to jump forward for prefetch
 * point. Multiplied by IL_VALUE at time of pf. >= 0, where 0
 * suppresses the prefetch.
 *
 * NOTE: At some point this was set to default to IL_VALUE, which was
 * wrong ehen IL_VALUE was used to adjust it later. Some timings will
 * reflect this bug. Fixed 2021-11-18.
 **/
#ifndef DEALER_PF_LEAD
#define DEALER_PF_LEAD 1
#endif
/**
 * If dealer does a "far" pf from time to time, this is > 0. Normally
 * should target 1 or 2 cache lines ahead of the current store
 * position is my guess.
 **/
#ifndef DEAL_PF_FAR
#define DEAL_PF_FAR -1
#endif
/**
 * If NT_DEAL is true the DIRECT_STORE_DEALER_CORE_WITH_PF dealer will
 * use an NT store op if WIDE_INDICES. Default false.  
 *
 * ***** Confirmed my earlier try. The 64 bit NT store sucks. It has a
 * very large latency, but that shouldn't matter unless it locks the
 * source register that long?? At E7 its 2.5 times as costly. At E8
 * it's 5 times (only a single trial, but...).
 **/
#ifndef NT_DEAL
#define NT_DEAL false
#endif

#ifndef PREFER_NEW_BUFFER
#define PREFER_NEW_BUFFER false
#endif
/**
 *                    Prefetch ladle CLs from time to time
 *
 * The ladle should be in cache all the time, but will often not be in
 * L1 cache, especially when BatchSize is not small. So it may be
 * useful to prefetch ladle CLs ahead during the deal. The dealer
 * pauses every so often to prefetch another CL for each ladle. The
 * key question is how many deal operations are required before we can
 * expect the prefetched CL to usually be in
 * place. LDL_N_IN_DEAL is set to that value. The prefetch algorithm
 * takes the NBuckets value into consideration when deciding how far
 * ahead to prefetch within each ladle. If the number of ladles is
 * large, we may only need to prefetch the next CL for each ladle. If it
 * is smaller, we need to prefetch multiple CLs to give each ladle
 * room for LDL_N_IN_DEAL / NBuckets items between prefetches.
 *
 * Note that when LDL_N_IN_DEAL / NBuckets is small, using it as the
 * lead amount per bucket means that we are usually not prefetching a
 * new CL at all. When larger, we usually are. Timings show the larger
 * LDL_N_IN_DEAL works better (up to about 9 per bucket). But this
 * doesn't really address whether we need to prefetch further ahead
 * than the next LDL_N_IN_DEAL. EXTRA_N_PER_LDL is the amount added to
 * (LDL_N_IN_DEAL / NBuckets) to support looking at this. (This
 * question arises only if LDL_N_IN_DEAL isn't a multiple of
 * NBuckets. It probably should be, so I think EXTRA_N_PER_LDL isn't
 * very useful.)
 *
 * Since LDL_N_IN_DEAL / NBuckets will be truncated to an int, it is
 * probably best if it is x * 2**n, where n is at least 5 or 6. I use
 * n = 8 for initial testing.  In principle, the prefetch should
 * reflect the size of the sort record, but I am testing assuming an
 * 8-byte record.
 *
 * The IL_VALUE also affects this. In principle, we only need to
 * prefetch the largest ladle in a bin. But that gets a bit messy,
 * since the mapping of ladles to bins isn't known in the ladle
 * manager where I do this. And of course the LDL_N_IN_DEAL is
 * multiplied by IlN.
 *
 * At the moment, this is only used at the start of a flush. The
 * prefetches during a deal use DEALER_PF_LEAD and ignore this.
 **/

// If LDL_N_IN_DEAL < 0, suppress this prefetch.
#ifndef LDL_N_IN_DEAL
#define LDL_N_IN_DEAL -1
// (5 * 256) looked useful at one point.
#endif
#ifndef EXTRA_N_PER_LDL
#define EXTRA_N_PER_LDL 0
#endif

/**
 * Select competing versions of starting ladle prefetches done by
 * DoneThisFlush(). Choices:
 * - 0: No ladle prefetch
 * - 1: Use PfLadleSpace()
 * - 2: Use PfCurStorePoint()
 * - 3: Use both.
 **/
#ifndef DTF_PREFETCH_OPTION
#define DTF_PREFETCH_OPTION 2
#endif

#ifndef USE_SIMD_PREP
#define USE_SIMD_PREP false
#endif

/**
 * When diverting, we select the number of passes to use before
 * diverting based on the number of items to be sorted. Best results
 * seem to come when N divided by the total number of buckets across
 * all passes (i.e., the product of the bucket sizes of the passes
 * used) is a smallish, but not too small, value. Stuart's figures
 * indicate that we should go from 3 passes to 4 at about
 * N=1.5E8. Since 256 cubed comes to 16M, the smallish number in this
 * case would be about 9.
 *
 * An array holds the minimum N to use for p passes, stored from
 * lowest to highest. We interogate these values in ascending order:
 * If the array was { 2000, 600000, 150000000, 0 }, we would use 1
 * pass for n up to 2000, 2 passes for 2000 <= n < 600000, etc. When
 * we find a 0 in the array, the number of passes it the count
 * associated with the last non zero n. In the above example, any n >=
 * 150000000 will use 4 passes.
 *
 * Do to the above, providing configuration constants for each limit
 * seems to be an accident waiting to happen since changing a value
 * without considering its neighbors in the array may lead to
 * unexpected pass counts. Thus the define for initializing this array
 * it a comma separated list of values. If provided, at least 1 value
 * must be given. At this time, we expect no more than 8, but fewer
 * than 8 is permitted.
 **/
/**
 * My desktop (Intel(R) Core(TM) i5-7600 CPU @ 3.50GHz):
 *
 * 1/2 pass boundry: Not much difference between 1000, 2000 and
 * 3000. Insertion sort seems to be only a small part of the cost.
 *
 * 2/3 pass boundry: Not much difference around 60000. Not researched
 * in detail, but 60000 seems pretty good.
 **/

#ifndef DIVERSION_ARRAY
#define DIVERSION_ARRAY 2000, 60000, 16000000, 0
#endif
// Following configurable in principle, but why?
#define DIVERSION_ARRAY_SIZE 8

// Default is no timing. Values : 0 (none), 1 (total only), >1 (pass
// and total times)
#ifndef TIMINGS
#define TIMINGS 0
#endif

#if 0
// Skeleton for setting defaults like above
#ifndef
#define
#endif
#endif

// ConstConfig                                                  ConstConfig

/**
 * Namespace ConstConfig collects the static const values that may be
 * specified at compile time, and a few values and types that are
 * dependent only on them. An SdConfig class variant (selected by
 * CONFIG_PATTERN) will be added to ConstConfig in Config.H.
 **/
namespace ConstConfig {
  
  // Declare values as static consts
  static constexpr int IlN = IL_VALUE;
  static constexpr int IlNShift = FindPower(IlN);
  
  static constexpr int Trigger = TRIGGER_VALUE;
  static constexpr int BatchSize = BATCH_SIZE;
  // Used variable NBuckets in a job to sample behaviour. That
  // requires Max_N_BUCKETS to allocate vectors, but precludes fixing
  // NBuckets at this level. Since we may use variable NBuckets later,
  // I'm leaving it that way.
  static constexpr int MaxNBuckets = Max_N_BUCKETS;
  static constexpr int ReadPfAdvance = READ_PF_ADVANCE;
  static constexpr int LdlNToPrime = LDL_N_TO_PRIME;
  static constexpr int BktNToPrime = BKT_N_TO_PRIME;
  static constexpr int BktLocality = BKT_LOCALITY;
  static constexpr int LdlNInDeal = LDL_N_IN_DEAL;
  static constexpr int ExtraNPerLdl = EXTRA_N_PER_LDL;
  static constexpr int DtfPrefetchOption = DTF_PREFETCH_OPTION;
  static constexpr int DealerPfLead = DEALER_PF_LEAD;
  static constexpr int DealerPfFreq = DEALER_PF_FREQ;
  static constexpr bool NtDeal = NT_DEAL;
  static constexpr int ExtendedRtPerMu = EXTENDED_RT_PER_MU;
  // With IlN > 1 and simd flush, which method to load the simd reg? {1,2,3}
  static constexpr int LoadSimdReg = LOAD_SIMD_REG;
    
  // Titan seems to prefer true here. Windows much prefers false!
  static constexpr bool PreferNewBuffer = PREFER_NEW_BUFFER;
  static constexpr bool WideIndices = WIDE_INDICES;
  static constexpr bool SystemPageSize = SYSTEM_PAGE_SIZE;
  static constexpr bool UseSimdPrep = USE_SIMD_PREP;
  static constexpr int TimingOption = TIMINGS;
  
  typedef COUNT_TYPE CountType;

  static bool DumpConfig = DUMP_CONFIG;  // May be overridden with -s

  namespace PrintConfigHelpers {  // Support for PrintConfig()
    void TestEndl(int n) {
      static constexpr int EndlTrigger = 60;
      static int LineW = EndlTrigger;
      if ((LineW + n + 4) >= EndlTrigger) {  
	std::cout << '\n';
	LineW = n + 4;
      } else {
	std::cout << ' ';
	LineW += n + 4;
      }
    }
    void PrintOne(const char *name, int value) {
      TestEndl(strlen(name) + 3);  // Mostly smallish ints
      std::cout << name << "= " << value;
    }
    void PrintOne(const char *name, double value) {
      TestEndl(strlen(name) + 6); 
      std::cout << name << "= " << value;
    }
    // Since builtin conversions preclude overloading 
    void PrintOneLong(const char *name, long value) {
      TestEndl(strlen(name) + 3);  // Mostly smallish ints
      std::cout << name << "= " << value;
    }
    void PrintOne(const char *name, bool value) {
      TestEndl(strlen(name) + 5);
      std::cout << name << "= ";
      if (value == 0) std::cout << "false";
      else if (value == 1) std::cout << "true";
      else std::cout << "???" << value << "???";
    }
    void PrintOne(const char *name, const char *value) {
      TestEndl(strlen(name) + strlen(value));
      std::cout << '{' << name << "= " << value << '}';
    }
  };

  void PrintConfig(long n_to_deal, const char *lm_name_p,
		   const char *fm_name_p, const char *dlr_core_name_p,
		   const char *f_mvr_name_p, const char *sd_name_p) {
    if (!DumpConfig) return;
    using namespace PrintConfigHelpers;
    
    PrintOne("Sort_Descriptor", sd_name_p);
    PrintOne("CONFIG_PATTERN", CONFIG_PATTERN);
    //PrintOneLong("n_to_deal", n_to_deal);  // Not a config!!
    PrintOne("LADLE_MANAGER", lm_name_p);
    PrintOne("FLUSHER_MANAGER", fm_name_p);
    //PrintOne("MOVER_MODULE", MOVER_MODULE);
    PrintOne("FLUSH_MOVER", f_mvr_name_p);
    PrintOne("DealerCore", dlr_core_name_p);
    PrintOne("BATCH_SIZE", BATCH_SIZE);
    PrintOne("USE_EFFECTIVE_BATCH_SIZE", USE_EFFECTIVE_BATCH_SIZE);
    PrintOne("TRIGGER_VALUE", TRIGGER_VALUE);
    PrintOne("IL_VALUE", IL_VALUE);
    PrintOne("EXTENDED_RT_PER_MU", EXTENDED_RT_PER_MU);
    PrintOne("LOAD_SIMD_REG", LOAD_SIMD_REG);    
    PrintOne("USE_COPY_INSERTION_SORT", USE_COPY_INSERTION_SORT);
    PrintOne("ISNS_SHORT_RUN", ISNS_SHORT_RUN);
    PrintOne("DO_IS_ON_DIP_SEGMENTS", DO_IS_ON_DIP_SEGMENTS);
    //PrintOne("FLUSH_BIN_PAIR_VN", FLUSH_BIN_PAIR_VN);
//PrintOne("ACPB_THRESHOLD", ACPB_THRESHOLD);
    PrintOne("PREFER_NEW_BUFFER", PREFER_NEW_BUFFER);
    PrintOne("READ_PF_ADVANCE", READ_PF_ADVANCE);
    PrintOne("LDL_N_IN_DEAL", LDL_N_IN_DEAL);
    PrintOne("EXTRA_N_PER_LDL", EXTRA_N_PER_LDL);
    PrintOne("DTF_PREFETCH_OPTION", DTF_PREFETCH_OPTION);
    PrintOne("LDL_N_TO_PRIME", LDL_N_TO_PRIME);
    PrintOne("BKT_N_TO_PRIME", BKT_N_TO_PRIME);
    PrintOne("BKT_LOCALITY", BKT_LOCALITY);
    PrintOne("DEALER_PF_LEAD", DEALER_PF_LEAD);
    //PrintOne("DEAL_PF_FAR", DEAL_PF_FAR);
    PrintOne("DEALER_PF_FREQ", DEALER_PF_FREQ);
    //PrintOne("NT_DEAL", NT_DEAL);
    PrintOne("Max_N_BUCKETS", Max_N_BUCKETS);
    PrintOne("WIDE_INDICES", WIDE_INDICES);
    PrintOne("SYSTEM_PAGE_SIZE", SYSTEM_PAGE_SIZE);
    PrintOne("DC_POST_INC", DC_POST_INC);
    PrintOne("USE_SIMD_PREP", USE_SIMD_PREP);
    PrintOne("BUCKET_EXPANSION_MULT", BUCKET_EXPANSION_MULT);
    TestEndl(999);  // Force end of line
    //PrintOne("", );

  }

};

#endif

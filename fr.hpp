#ifndef FASTRADIX_HPP
#define FASTRADIX_HPP

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>

#ifdef DEBUG
#include <bitset>
#include <iomanip>
#endif

#ifdef TIMER
#include <chrono>
#endif


#ifdef _WINDOWS
typedef unsigned __int64 uint64_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int8 uint8_t;
#else
#include <stdint.h>
#endif

//This should be ~14 lowering while I dev

#ifdef TEST_THRESHOLD
const int DIVERSION_THRESHOLD=TEST_THRESHOLD;
#else
const int DIVERSION_THRESHOLD = 12;
#endif

#ifdef TEST_STDSORT_THRESHOLD
const int STDSORT_DIVERSION_THRESHOLD=TEST_STD_SORTTHRESHOLD;
#else
const int STDSORT_DIVERSION_THRESHOLD = 250;
#endif

const int DISTRIBUTION_SENSITIVE_THRESHOLD = 4096; 	/* 	If a length of data to be processed is smaller than this and
						      		we have no other data, we don't want to fart around figuring
								out the actual distribution, just assume uniform and deal with 
								the consequences. This can come up on on the first pass of the
								top bytes or possibly the first pass of the bottom bytes.
							*/
const int LADLE_SIZE = 256;

#ifdef TEST_LADLE_THRESHOLD
const int LADLE_THRESHOLD = TEST_LADLE_THRESHOLD;
#else
const unsigned LADLE_THRESHOLD = 4000000000; //roughlt LADLE_SIZE*256 for a 1-to-1
#endif

const int PREFETCH_LOOKAHEAD = 1;

template <typename INT, typename ELEM>
void dfr(ELEM *source, auto length) {

	/**
		1) Diverting on small lengths
	*/

	auto insertionSort = [](register ELEM *source, const register auto &length) {	
        	register ELEM buf;
        	register INT val;
        	register unsigned i = 1;
        	register unsigned cursor;

        	for(; i < ((length <= DIVERSION_THRESHOLD)?length:DIVERSION_THRESHOLD); i++) {
        		cursor=i;
			buf = source[cursor];
			val = *(reinterpret_cast<INT*>(source + (cursor)));
			while(cursor > 0 &&
			val < *(reinterpret_cast<INT*>(source + (cursor-1)))
			) {
				source[cursor]=source[cursor-1];
				cursor--;
			}
			source[cursor]=buf;
        	}

        	for(; i < length; i++) {
			cursor = i;
                	buf = source[cursor];
                	val = *(reinterpret_cast<INT*>(source + (cursor)));
                	while(
                      	val < *(reinterpret_cast<INT*>(source + (cursor-1)))
                     	) {
                        	source[cursor]=source[cursor-1];
                        	cursor--;
                	}
                	source[cursor]=buf;
        	}
	};

	if(length <= DIVERSION_THRESHOLD) {
		insertionSort(source, length);
		return;
	}


	//The number of bytes in the data type used
	const unsigned int bytecount = sizeof(INT);

	//The buffer for doing radix passes
	ELEM *destination = new ELEM[length];

	//These will start as pairs of positions pointing to the start of estimated
	//buckets. As values are placed, the second entry will increment. After all
	//data is dealt, this will be a set of start/end points for each bucket
	//I suggest all start values, then all end values as a cache-efficient format
	//When estimating uniform buckets, we can use math and the first 256 values
	//To make this fancier/faster.
	//
	//We need two sets of these, one for each buffer, and we'll swap them
	ELEM **sourceBuckets = new ELEM*[512];
	ELEM **destinationBuckets = new ELEM*[512];

	//Initially these will be the starting positions into the overflow buffer
	//Once dumped into the overflow buffer, they will then be the ending positions
	//when dealing out of the overflow buffer
	ELEM **overflowBuckets = new ELEM*[256];

	ELEM **ladleBuckets = new ELEM*[256];

	//If we ever actually count elements upon placing (e.g. a last pass), we should
	//use this before converting it to *buckets*
	unsigned* bucketCounts = new unsigned[256];
        std::memset(bucketCounts, 0, 256*sizeof(unsigned)); //This is needed!

	//We actually count for overflow, so we must use this. We'll just convert to 
	//overflowBuckets when done because that's easier to process.
	unsigned* overflowCounts = new unsigned[256];
        std::memset(overflowCounts, 0, 256*sizeof(unsigned)); //This is needed!

	//We don't initialize this till we need it. We should check overflowMaxSize 
	//and expend it as needed, re-initializing the buffer. There may be benefit
	//to being fancy about how we choose to grow this buffer, bur for now I don't 
	//care. We may also do something like initializing it for 2% of data, which
	//should be a realistic expectation given the newer sampling approach.
	int overflowMaxSize = 0;
	ELEM *overflowBuffer = NULL;

	ELEM *ladleBuffer = NULL;

	const auto resetLadleBuckets = [&ladleBuffer, &ladleBuckets] {
		ladleBuffer = new ELEM[LADLE_SIZE*256];
		ladleBuckets[0] = ladleBuffer;
		for(unsigned i = 1; i < 256; i++) ladleBuckets[i] = ladleBuckets[i-1]+LADLE_SIZE;
	};

	if(length > LADLE_THRESHOLD) {
		resetLadleBuckets();
	}

	/**
		3) Estimating relevant top (high-order) bytes
	*/

	//For scanning:
	//checking for live bits
	INT bitmask = source[0];
	INT livebits = 0;

        // Random pieces
        std::random_device rd;
        std::default_random_engine gen(rd());
        std::uniform_int_distribution<> dist(1, length);


	//After just 4 samples (3+initial bitmask), we know pretty reasonably
	//which bits have even odds of being there. only a 1/8 chance that
	//a even chance bit is not shown. On average, that means 8 bits are
	//going to be mis-represented, but some of them will be in the low
	//order, where it has no impact. If a bit is live, but there's a very
	//low chance that  it would be detected, then we actually want to under-
	//represent it because it's not going to spread the sizes as evenly, so
	//it is to our advantage to under-represent that case with the liveness
	//check.

	#ifdef NOTRANDOM
        livebits |= bitmask^ *(reinterpret_cast<INT*>(source + 1));
        livebits |= bitmask^ *(reinterpret_cast<INT*>(source + 2));
        livebits |= bitmask^ *(reinterpret_cast<INT*>(source + 3));
	#else
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	#endif

	auto neededBits = ceil(std::log2((float)(length)/(float)(DIVERSION_THRESHOLD))); //ceil(log(1792/14)/log(2)),
	//auto neededBits = ceil(std::log2((float)(length)/(float)(STDSORT_DIVERSION_THRESHOLD-1)));
	//The above causese a bug. I don't understand why...

	auto foundLiveBits = 0;
	auto neededBytes = sizeof(INT);
	const auto allBits = 8*sizeof(INT);

	//STUART:: look into __builtin_popcount or __builtin_popcountll 
	for(auto i = 1; i < allBits; i++) {
		auto val = ((livebits >> (int)(allBits - i)) & 1);
		foundLiveBits += val;
		if(foundLiveBits==neededBits) {
			neededBytes = ceil((i+1)/8.0);
			break;
		}
	}

	int countedByte;
	const auto topBytes = new uint8_t[sizeof(INT)];
	std::memset(topBytes, 0, sizeof(INT)*sizeof(uint8_t));
	auto topBytesSize = 0;
	const auto bottomBytes = new uint8_t[sizeof(INT)];
	std::memset(bottomBytes, 0, sizeof(INT)*sizeof(uint8_t));
        auto bottomBytesSize = 0;

	const INT myi = 1;
	const INT highBitMask = (((((myi << (allBits-1))-1)<<1)+1)<<((sizeof(INT)-neededBytes)*8) );
        auto overflow = source;

	/**
		4) GETTING A SAMPLE DISTRIBUTION
	*/

	/*
		We divide the length roughly (using bitshift) by 256 and
		allocate the space in destination evenly.

		The extra space after the rough divide is distributed by
		adding one spot to each of the first few buckets (the
		remainder) and then just allocating the rougher estimate
		for the latter buckets, thus the full length is allocated.
	*/
	auto estimateUniform = [] (const auto &destination, const auto &length, const auto &destinationBuckets) {
		const auto bucketGuess = length>>8;
		const auto remainder = length-(bucketGuess<<8);
		auto d = destination;
		auto buildBucketEstimates = [&destinationBuckets, &d] (const auto &bucketSize, const auto &start, const auto &end, const auto &destinationBuckets){
			for(auto i = start; i < end; i++) {
				destinationBuckets[i << 1]       = d;//start
				d += bucketSize;		//these buckets are one bigger than the latter buckets
				destinationBuckets[(i << 1) + 1] = d;//end, which will be start of next bucket
			}
		};

		buildBucketEstimates(bucketGuess+1, 0, remainder, destinationBuckets);
		buildBucketEstimates(bucketGuess, remainder, 256, destinationBuckets);
	};

	/*
		This will eventually need to capture more data
	*/
	auto doEstimates = [&estimateUniform, &destination, &length, &destinationBuckets]() {
		if(length < DISTRIBUTION_SENSITIVE_THRESHOLD) {
			estimateUniform(destination, length, destinationBuckets);
		} else {
			//THIS IS NOT DONE YET! 4b!
			estimateUniform(destination, length, destinationBuckets);
		}
	};

	auto swap = [&source, &destination, &sourceBuckets, &destinationBuckets, &length] {
		std::swap(source, destination);
		std::swap(sourceBuckets, destinationBuckets);
	};

	const auto processLadle = 
		[&ladleBuffer, &ladleBuckets, &destinationBuckets, &overflow, &overflowBuffer, &overflowMaxSize, &source, &destination, &overflowCounts, &length]
		(const auto &currentByte) {

#ifdef DEBUG_LADLE
std::cout << 
"Processing ladle for byte " << +currentByte
<< std::endl;
#endif
		for(auto i = 0; i < 256; i++) {
			__builtin_prefetch((const void*)(destinationBuckets[i]),0,1);
			ELEM *start = ladleBuffer+i*LADLE_SIZE;
			__builtin_prefetch((const void*)(start),0,0);
			ELEM *end = ladleBuckets[i];
			unsigned len = end-start;
			if(len > 0) {
#ifdef DEBUG_LADLE
std::cout << "ladle bucket " << i << ": ";
for(int i = 0; i < len; i++) std::cout << *(start+i) << " ";
std::cout << std::endl;
#endif

				const uint8_t *target = reinterpret_cast<const uint8_t*>(start) + currentByte;
				ELEM **currentDestination = destinationBuckets + ((*target) << 1);
				unsigned lengthLeft = *(currentDestination+1)-*currentDestination;
				if(len <= lengthLeft) { //we fit
#ifdef DEBUG_LADLE
std::cout << "The data chunk of " << len << " elements fits in it's regular bucket located at " << +(*currentDestination)-destination << std::endl;
#endif
					std::memcpy(*currentDestination, start, len*sizeof(ELEM));
					*currentDestination += len;
#ifdef DEBUG_LADLE
std::cout << "\tIts regular bucket is now located at " << +(*currentDestination)-destination << std::endl;
#endif

				} else {//we don't fit. paste as much as we can, then use overflow
#ifdef DEBUG_LADLE
std::cout << "The data chunk of " << len << " elements did not fit in its regular bucket located at " << +(*currentDestination)-destination << " so we could only deal in " << lengthLeft << " elements" << std::endl;
#endif
					const unsigned overflowDeal = (len-lengthLeft);
					std::memcpy(*currentDestination, start, lengthLeft*sizeof(ELEM));
					*currentDestination += lengthLeft;
#ifdef DEBUG_LADLE
std::cout << "\tIts regular bucket is now located at " << +(*currentDestination)-destination << std::endl;
#endif

					if((overflowBuffer+overflowMaxSize <= overflow) && (overflow <= overflowBuffer+(2 * overflowMaxSize))) {
						unsigned overflowLeft = overflowBuffer+(2 * overflowMaxSize)-overflow;
						if((len-lengthLeft) <= overflowLeft) { //fits in initial overflow
#ifdef DEBUG_LADLE
std::cout << "\tIt did fit into the regular overflow at " << overflow - (overflowBuffer+overflowMaxSize) << " where we dealt " << overflowDeal << " elements" << std::endl;
#endif
							std::memcpy(overflow, start+lengthLeft, overflowDeal*sizeof(ELEM));
							overflow+= overflowDeal;
						} else {
#ifdef DEBUG_LADLE
std::cout << "\tIt did not fit into the regular overflow at " << overflow - (overflowBuffer+overflowMaxSize) << " so we could only deal  " << overflowLeft << " elements" << std::endl;
#endif
							std::memcpy(overflow, start+lengthLeft, (overflowLeft)*sizeof(ELEM));

#ifdef DEBUG_LADLE
std::cout << "\tWe then dealt the remaining " << (len-lengthLeft-overflowLeft) << " elements to the beginning of source" << std::endl;
#endif

							std::memcpy(source, start+lengthLeft+overflowLeft, (overflowDeal-overflowLeft)*sizeof(ELEM));
							overflow=source+(overflowDeal-overflowLeft);
						}
					} else { // we have to put it in overflow overflow, which always fits.
#ifdef DEBUG_LADLE
std::cout << "\tWe then dealt the remaining " << overflowDeal <<  " elements into the overflow overflow (the beginning of source) " << overflow-source << std::endl;
#endif
						std::memcpy(overflow, start+lengthLeft, overflowDeal*sizeof(ELEM));
						overflow += overflowDeal;
					}
					//Either way the overflow count is adjusted for everything that didn't get dealt into the bucket.
					overflowCounts[*target]+=overflowDeal;
				}
				ladleBuckets[i] = start;
			}
		}
#ifdef DEBUG_LADLE
std::cout << "destination: ";
for(int i = 0; i < length; i++) std::cout << *(destination+i) << " ";
std::cout << std::endl;
std::cout << "overfow(source): ";
for(int i = 0; i < overflow-source; i++) std::cout << *(source+i) << " ";
std::cout << std::endl;
#endif

	};


        const auto passOverInputDealWithOverflowAndGatherLiveBits = 
		[&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize, &livebits, &bitmask]
		(ELEM* start, ELEM *end, const register uint8_t &currentByte, const register auto &targetBuckets) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const register unsigned len = (end-start);
		register ELEM **currentDestination;

                for(register unsigned i = 0; i < len; ++i, livebits |= bitmask ^ reinterpret_cast<INT>(*(element++)), target+=sizeof(ELEM)) {
			currentDestination = targetBuckets + ((*target) << 1);
			if(*currentDestination < *(currentDestination+1)) {
				*((*currentDestination)++) = *element;
 			} else {
				overflowCounts[*target]++;
				*overflow = *element;
				overflow++;
				if(overflow == overflowBuffer+(2 * overflowMaxSize)) {
					overflow = source;
				}
			}
                }
        };

        const auto passOverInputDealWithOverflowAndGatherLiveBitsLadle =
                [&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize, &livebits, &bitmask, &processLadle, &ladleBuckets, &length, &ladleBuffer]
                (ELEM* start, ELEM *end, const register uint8_t &currentByte, const register auto &targetBuckets) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const register unsigned len = (end-start);
		register ELEM **currentDestination;

		if(length > LADLE_THRESHOLD) {
			for(register auto i = 0; i < len; ++i, livebits |= bitmask ^ reinterpret_cast<INT>(*(element++)), target+=sizeof(ELEM)) {
				if(*(ladleBuckets + (*target)) == (ladleBuffer + (*target + 1)*LADLE_SIZE )) {
					processLadle(currentByte);
					*((*(ladleBuckets + (*target)))++) = *element;
				} else {
					*((*(ladleBuckets + (*target)))++) = *element;
				}
			}
			processLadle(currentByte);
		} else {
                	for(register auto i = 0; i < len; ++i, livebits |= bitmask ^ reinterpret_cast<INT>(*(element++)), target+=sizeof(ELEM)) {
                        	currentDestination = targetBuckets + ((*target) << 1);
                        	if(*currentDestination < *(currentDestination+1)) {
                                	*((*currentDestination)++) = *element;
                        	} else {
                                	overflowCounts[*target]++;
                                	*overflow = *element;
                                	overflow++;
                                	if(overflow == overflowBuffer+(2 * overflowMaxSize)) {
                                        	overflow = source;
                                	}
                        	}
                	}
		}
        };


        const auto passOverInputDealWithOverflow =
                [&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize]
                (ELEM* start, ELEM *end, const register uint8_t &currentByte, const register auto &targetBuckets) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const register unsigned len = (end-start);
		register ELEM **currentDestination;

                for(register auto i = 0; i < len; ++i, element++, target+=sizeof(ELEM)) {
                        currentDestination = targetBuckets + ((*target) << 1);
                        if(*currentDestination < *(currentDestination+1)) {
                                *((*currentDestination)++) = *element;
                        } else {
                                overflowCounts[*target]++;
                                *overflow = *element;
                                overflow++;
                                if(overflow == overflowBuffer+(2 * overflowMaxSize)) {
                                        overflow = source;
                                }
                        }
                }
        };

        const auto passOverInputsDealWithOverflow = 
		[&sourceBuckets, &overflowBuckets, &passOverInputDealWithOverflow]
		(ELEM *thisSource, ELEM *thisBuffer, const register auto &currentByte, const register auto &targetBuckets) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(auto i = 0; i < 256; i++) {
                                endSourceBucket = sourceBuckets[i<<1];
                                endOverflowBucket= overflowBuckets[i];

                                if(endSourceBucket-startSourceBucket) {
					passOverInputDealWithOverflow(startSourceBucket, endSourceBucket, currentByte, targetBuckets);
                                }

                                if(endOverflowBucket-startOverflowBucket) {
                                        passOverInputDealWithOverflow(startOverflowBucket, endOverflowBucket, currentByte, targetBuckets);
				}

				startSourceBucket=sourceBuckets[(i<<1)+1];
				startOverflowBucket = endOverflowBucket;
                        }

        };

        const auto passOverInputDealWithOverflowAndCounting =
                [&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize, &bucketCounts, &destination]
                (ELEM* start, ELEM *end, const register uint8_t &currentByte, const register auto &targetBuckets, const uint8_t &countByte) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
		const register uint8_t *countTarget = reinterpret_cast<const uint8_t*>(element) + countByte;
                const register unsigned len = (end-start);
		register ELEM **currentDestination;

                for(register auto i = 0; i < len; ++i, element++, target+=sizeof(ELEM), bucketCounts[*countTarget]++, countTarget+=sizeof(ELEM)) {
         		currentDestination = targetBuckets + ((*target) << 1);
                        if(*currentDestination < *(currentDestination+1)) {
                                *((*currentDestination)++) = *element;
                        } else {
                                overflowCounts[*target]++;
                                *overflow = *element;
                                overflow++;
                                if(overflow == overflowBuffer+(2 * overflowMaxSize)) {
                                        overflow = source;
                                }
                        }
                }
        };

        const auto passOverInputsDealWithOverflowAndCounting =
                [&sourceBuckets, &overflowBuckets, &passOverInputDealWithOverflowAndCounting]
                (ELEM *thisSource, ELEM *thisBuffer, const auto &currentByte, const register auto &targetBuckets, const register uint8_t &countByte) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(auto i = 0; i < 256; i++) {
				endSourceBucket = sourceBuckets[i<<1];
				endOverflowBucket= overflowBuckets[i];

				if(endSourceBucket-startSourceBucket) {
                                        passOverInputDealWithOverflowAndCounting(startSourceBucket, endSourceBucket, currentByte, targetBuckets, countByte);
				}
				if(endOverflowBucket-startOverflowBucket) {
                                        passOverInputDealWithOverflowAndCounting(startOverflowBucket, endOverflowBucket, currentByte, targetBuckets, countByte);
				}
                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };

        const auto passOverInputDealExact =
                []
                (ELEM* start, ELEM *end, const register uint8_t &currentByte, const register auto &targetBuckets) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const register unsigned len = (end-start);
		register ELEM **currentDestination;

                for(register auto i = 0; i < len; ++i, ++element, target+=sizeof(ELEM)) {
			currentDestination = targetBuckets + (*target);
			*((*currentDestination)++) = *element;
                }
        };

        const auto passOverInputsDealExact =
                [&sourceBuckets, &overflowBuckets, &passOverInputDealExact]
                (ELEM *thisSource, ELEM *thisBuffer, const register auto &currentByte, const register auto &targetBuckets) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(auto i = 0; i < 256; i++) {
                                        passOverInputDealExact(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], currentByte, targetBuckets);
                                        passOverInputDealExact(startOverflowBucket, endOverflowBucket= overflowBuckets[i], currentByte, targetBuckets);

                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };

        const auto passOverInputCounting =
                [&bucketCounts]
                (ELEM* start, ELEM *end, const register uint8_t &countByte) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *countTarget = reinterpret_cast<const uint8_t*>(element) + countByte;

                const register unsigned len = (end-start);

                for(register unsigned i = 0; i < len; ++i, element++, bucketCounts[*countTarget]++, countTarget+=sizeof(ELEM));
        };

        const auto passOverInputsCounting =
                [&sourceBuckets, &overflowBuckets, &passOverInputCounting]
                (ELEM *thisSource, ELEM *thisBuffer, const uint8_t &countByte) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(unsigned i = 0; i < 256; i++) {
                                        passOverInputCounting(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], countByte);
                                        passOverInputCounting(startOverflowBucket, endOverflowBucket= overflowBuckets[i], countByte);

                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };


        const auto passOverInputDealExactAndGatherLiveBits =
                [&livebits, &bitmask]
                (ELEM* start, ELEM *end, const register uint8_t &currentByte, const register auto &targetBuckets) {
                if(start == end) return;
                const register ELEM* element = start;
                const register uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const register unsigned len = (end-start);
		register ELEM **currentDestination;

                for(register auto i = 0; i < len; ++i, livebits |= bitmask ^ reinterpret_cast<INT>(*(element++)), target+=sizeof(ELEM)) {
			currentDestination = targetBuckets + (*target);
			*((*currentDestination)++) = *element;
                }
        };

        const auto passOverInputsDealExactAndGatherLiveBits =
                [&sourceBuckets, &overflowBuckets, &passOverInputDealExactAndGatherLiveBits]
                (ELEM *thisSource, ELEM *thisBuffer, const auto &currentByte, const auto &targetBuckets) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(auto i = 0; i < 256; i++) {
                                        passOverInputDealExactAndGatherLiveBits(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], currentByte, targetBuckets);
                                        passOverInputDealExactAndGatherLiveBits(startOverflowBucket, endOverflowBucket= overflowBuckets[i], currentByte, targetBuckets);

                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };

	//STUART:: I noted that extracting the darget buckets pointer to a register yielded advantage elsewhere... is this true here too??
        const auto passOverInputDealSimple =
                []
                (ELEM* start, ELEM *end, const register auto &targetBuckets) {
                if(start == end) return;
                const register ELEM* element = start;
                const register unsigned len = (end-start);

                for(register auto i = 0; i < len; ++i, *((*targetBuckets)++)=*(element++));
        };

        const auto passOverInputsDealSimple =
                [&sourceBuckets, &overflowBuckets, &passOverInputDealSimple]
                (ELEM *thisSource, ELEM *thisBuffer, const auto &targetBuckets) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(auto i = 0; i < 256; i++) {
                                        passOverInputDealSimple(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], targetBuckets);
                                        passOverInputDealSimple(startOverflowBucket, endOverflowBucket= overflowBuckets[i], targetBuckets);

                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };

	const auto convertCountsToContiguousBuckets = [&destination](const auto &counts, const auto &buckets, const auto &buffer) {
		*buckets=buffer;
		register auto currentBucket = buckets+1;
		register auto previousBucket = buckets+0;

		std::for_each(counts, counts+255, [&currentBucket, &previousBucket](register auto &count) {
			*currentBucket = *previousBucket + count;
			currentBucket++;
			previousBucket++;
			count=0;
		});
	};

	const auto processOverflow = [&](const auto &currentByte){
		//Check if overflow is in source or if it still fits within existing buffer.
		if((source < overflow) && ( overflow <= (source+length))) {
			//We used up the old overflow so we need a bigger overflow.

			//Make a new buffer
			int newsize = overflowMaxSize + overflow-source;
			ELEM *oldOverflowBuffer = overflowBuffer;

			overflowBuffer = new ELEM[2*newsize];

			//build overflow buckets
			convertCountsToContiguousBuckets(overflowCounts, overflowBuckets, overflowBuffer);

			//pass left-to right from the overflow buffer, then from the front of the source also used as overflow buffer.
			passOverInputDealExact(oldOverflowBuffer+overflowMaxSize, oldOverflowBuffer+(2*overflowMaxSize), currentByte, overflowBuckets);
			if(overflow-source)passOverInputDealExact(source, overflow, currentByte, overflowBuckets);

			//Assign the new buffer size.
			overflowMaxSize = newsize;
			//point overflow at the new space
			overflow = overflowBuffer+newsize;

			//Kill old buffer
			if(oldOverflowBuffer) {
				delete [] oldOverflowBuffer;
			}
		} else {
			if(overflowMaxSize != 0){ //If we processed overflow without ever using overflow we skip this stuff
				convertCountsToContiguousBuckets(overflowCounts, overflowBuckets, overflowBuffer);
				passOverInputDealExact(overflowBuffer+overflowMaxSize,(source==overflow)?(overflowBuffer+(2*overflowMaxSize)):overflow, currentByte, overflowBuckets);
				overflow = overflowBuffer+overflowMaxSize;
			}
		}
	};

        // Allocates bytes, in order, to either the list of Top Bytes
        // or Bottom Bytes, based on the Threshold Byte which doesn't
        // get added to either because presumably that byte was used
        // directly or indirectly to scan for live bits.
        //
        // This is so cheap let's not do anything fancy
	auto buildByteLists = [&](const auto &threshold) {
		for(auto i = 0; i < bytecount; i++) {
			if((livebits >> (i*8)) & 255) {
				if(i < threshold) {
					bottomBytes[bottomBytesSize] = i;
					bottomBytesSize++;
				} else if(i >= threshold) {
					topBytes[topBytesSize] = i;
					topBytesSize++;
				}
			}
		}
	};

	bool hasByteLists = false;
 	auto fr = [&](const auto &bytes, const auto &numBytes) {
		auto currentByte = hasByteLists?bytes[0]:sizeof(INT)-neededBytes;
		if(hasByteLists) neededBytes = numBytes;
		bool setByteListHere = false;

	//Should we assign overflow's starting position here? Is the shrinking of source the cause of this pain or is it always pointing to a suitable buffer?

#ifdef TIMINGS
        std::chrono::high_resolution_clock::time_point end;
        std::chrono::high_resolution_clock::time_point start;
#endif


	        if(neededBytes == 1) { // Just do a single count pass first and in a very un-fast-radix-way deal exactly into destination
			countedByte=currentByte;
			passOverInputCounting(source, source+length, currentByte);
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

			if(hasByteLists) {
				passOverInputDealExact(source, source+length, currentByte, destinationBuckets);
			} else {
				passOverInputDealExactAndGatherLiveBits(source, source+length, currentByte, destinationBuckets);
				buildByteLists(currentByte);
				hasByteLists=true;
#ifdef DEBUG
std::cout << "After 1 pass Top Bytes: " << topBytesSize << " Bottom Bytes: " << bottomBytesSize << std::endl;
#endif

			}
			swap();
	        } else if(neededBytes == 2) {   // We need a pass that does count capturing, then deals back in place.
						// We ignore when highest byte is dead; too annoying to check for

#ifdef DEBUG
#ifdef TIMINGS
start = std::chrono::high_resolution_clock::now();
#endif
#endif

			countedByte = hasByteLists?bytes[1]:currentByte+1;
			doEstimates();
			passOverInputDealWithOverflowAndCounting(source, source+length, currentByte, destinationBuckets, countedByte);
			processOverflow(currentByte);
			swap();
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "\tPass 1 of 2 Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
start = std::chrono::high_resolution_clock::now();
#endif
#endif


			if(hasByteLists) {
				passOverInputsDealExact(source, overflowBuffer, countedByte, destinationBuckets);
			} else {
				passOverInputsDealExactAndGatherLiveBits(source, overflowBuffer, countedByte, destinationBuckets);
				buildByteLists(currentByte);
				hasByteLists=true;
			}
			swap();




#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "\tPass 2 of 2 Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif
#endif

#ifdef DEBUG
std::cout << "After 2 passes Top Bytes: " << topBytesSize << " Bottom Bytes: " << bottomBytesSize << std::endl;
#endif



	        } else {
			//Do a first pass
#ifdef TIMINGS
        start = std::chrono::high_resolution_clock::now();
#endif
			doEstimates();
			if(hasByteLists) {
					passOverInputDealWithOverflow(source, source+length, currentByte, destinationBuckets);
			} else {
					passOverInputDealWithOverflowAndGatherLiveBitsLadle(source, source+length, currentByte, destinationBuckets);
					buildByteLists(currentByte);
					hasByteLists=true;
					setByteListHere=true;
#ifdef DEBUG
std::cout << "Top Bytes: " << topBytesSize << " Bottom Bytes: " << bottomBytesSize << std::endl;
#endif

			}

			processOverflow(currentByte);
			swap();
#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "First Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif
#endif



		//We have to account for the newly generated numBytes being either 0, 1, 2 or 3 (we did one pass, but numbytes still has it... 
		//our one pass could theoretically have been over dead bits), which means we found enough dead bits to eliminate
		//passes that would have got us caught up above in the numBytes == 1 or numBytes == 2 section. Sure it's annoying, but actually
		//doing the deadbit passes would be more annoying, so rejoice!

		if(setByteListHere && currentByte != bytes[0]) { //We did a pass that wasn't needed. Reverse and recurse
			destinationBuckets[0] = destination; //simpleDeal needs this
			//passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
			passOverInputsDealSimple(source, overflowBuffer, destinationBuckets);
			swap();
			return true; //Tell the caller to re-call fr now that we're "fixed" things.
		}

		if (numBytes == 1) {	//We did a pass with overflow but we were only supposed to do one pass anyway
						// Since we didn't count we need to fix that
				//Ok, we did a pass, it was needed, but we should have made it an
				// exact pass. We have to fix that here.
			destinationBuckets[0] = destination; //simpleDeal needs this
			//passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
			passOverInputsDealSimple(source, overflowBuffer, destinationBuckets);
			swap();
			std::memcpy(destination, source, length*sizeof(ELEM));
			swap();

		} else if (numBytes == 2) {	//We did one of two passes. 
			//We did a real pass, it just would have helped to have counted first.
			countedByte = bytes[1];
			//passOverInputs(source, overflowBuffer, countedByte, noDeal, countForByte);
			passOverInputsCounting(source, overflowBuffer, countedByte);
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);
			//passOverInputs(source, overflowBuffer, countedByte, dealExact, noStats);
			passOverInputsDealExact(source, overflowBuffer, countedByte, destinationBuckets);
			swap();
		} else {

				for(int i = 1; i < (numBytes - 2); i++) {
#ifdef DEBUG
#ifdef TIMINGS
        start = std::chrono::high_resolution_clock::now();
#endif
#endif

						doEstimates();
						passOverInputsDealWithOverflow(source, overflowBuffer, bytes[i], destinationBuckets);
						processOverflow(bytes[i]);
						swap();
#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "Middle Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif
#endif

				}


				//Do last two passes
#ifdef DEBUG
#ifdef TIMINGS
        start = std::chrono::high_resolution_clock::now();
#endif
#endif

				countedByte=bytes[numBytes-1];
				doEstimates();
				passOverInputsDealWithOverflowAndCounting(source, overflowBuffer, bytes[numBytes-2], destinationBuckets, bytes[numBytes-1]);
				processOverflow(bytes[numBytes-2]);
				swap();
				convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "Next-To-Last Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
start = std::chrono::high_resolution_clock::now();
#endif
#endif


				//Deal last pass in top bytes
				passOverInputsDealExact(source, overflowBuffer, bytes[numBytes-1], destinationBuckets);
				swap();
#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "Last Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif
#endif

			}
		}

		return false;
	};


	//
	// Process  Top Bytes
	//

#ifdef DEBUG
#ifdef TIMINGS
        std::chrono::high_resolution_clock::time_point end;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
#endif
#endif

	if(fr(topBytes, topBytesSize)) fr(topBytes, topBytesSize);
#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "Time Top: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
start = std::chrono::high_resolution_clock::now();
#endif
#endif

	if(bottomBytesSize > 0) {


		//At this point, source has our data in disjoint buckets and destination is ready.
		//buckets contains the start and end positions of each disjoint bucket in source.
		//overflowBuffer is the contiguous secondary source, with end buckets defined by
		//overflowBuckets, and is twice overflowMaxSize so we can deal into it to before
		//we make use of the input as a buffer source to prevent writing over a live bucket.

		// ((reinterpret_cast<INT>(element))>

		
		const register INT r_highBitMask = highBitMask;
		const auto getHighBits = [&r_highBitMask](const register ELEM &element) {
			return (reinterpret_cast<INT>(element)) & r_highBitMask;
		};

		register INT currentBits = getHighBits(*source);

		register bool potentialLongRun = false;
		register bool definiteLongRun = false;
		ELEM* startSmallRuns;
		ELEM* endSmallRuns;
		ELEM* startDefiniteLongRun;
		ELEM* endDefiniteLongRun;
		const register unsigned step = (DIVERSION_THRESHOLD>>1);
		startSmallRuns = source;
		register ELEM* element = source+(DIVERSION_THRESHOLD>>1)-1;
		register INT newBits;
		for(register unsigned i = step-1; i < length-step; i += step, element += step) {
			newBits = getHighBits(*element);
			if(currentBits ^ newBits) { //Things changed
				if(definiteLongRun) { // A long run has ended
					//walk back by one to find the exact end of the long run
					endDefiniteLongRun=element;
					while(endDefiniteLongRun > source && (getHighBits(*(endDefiniteLongRun-1)) ^ currentBits)) endDefiniteLongRun--;

					//process the long run.
                                        if((endDefiniteLongRun-startDefiniteLongRun) < STDSORT_DIVERSION_THRESHOLD) {
                                                std::sort(startDefiniteLongRun, endDefiniteLongRun);
                                        } else {
						const auto oldSource = source;
						const auto oldDestination = destination;
						const auto oldLength = length;
						source = startDefiniteLongRun;
						destination = destination + (startDefiniteLongRun-source);
						length = endDefiniteLongRun - startDefiniteLongRun;

						if(overflowMaxSize>0) overflow = overflowBuffer + overflowMaxSize;
						else overflow = source;

						fr(bottomBytes, bottomBytesSize);
						if(bottomBytesSize%2) {
							std::memcpy(destination, source, length*sizeof(ELEM));
						}

						source = oldSource;
						destination = oldDestination;
						length = oldLength;
					}
					//Start a new short run
					startSmallRuns=endDefiniteLongRun;

					definiteLongRun=false;
				}
				currentBits = newBits;
				potentialLongRun = false;
			} else { // Duplicate hit, figure out if we're a long run or just beginning to be one
				if(!definiteLongRun) {
					if(potentialLongRun) {
						//Mark the potential start of a long run by walking back to find a change
						//Also run the previous short runs
						startDefiniteLongRun=element;
						while((startDefiniteLongRun > source) && !(getHighBits(*(startDefiniteLongRun-1)) ^ currentBits)) startDefiniteLongRun--;
						//Insertion Sort Long Runs
						if(startSmallRuns != startDefiniteLongRun) {
							endSmallRuns = startDefiniteLongRun;
							insertionSort(startSmallRuns, endSmallRuns-startSmallRuns);
						}
						definiteLongRun=true;
						potentialLongRun=false;
					} else {
						potentialLongRun = true;
					}
				} 
			}
		}

		// We're at the end of the list, figure out if the last run
		// is long or not then process.
		if(definiteLongRun) {
			endDefiniteLongRun = source+length;
                	//process the long run.
                	const auto oldSource = source;
                	const auto oldDestination = destination;
                	const auto oldLength = length;
                	source = startDefiniteLongRun;
                	destination = destination + (startDefiniteLongRun-source);
                	length = endDefiniteLongRun - startDefiniteLongRun;
			if(overflowMaxSize>0) overflow = overflowBuffer + overflowMaxSize;
			else overflow = source;

                	fr(bottomBytes, bottomBytesSize);
                	if((bottomBytesSize)%2) {
				std::memcpy(destination, source, length*sizeof(ELEM));
			}
                	source = oldSource;
                	destination = oldDestination;
                	length = oldLength;
		} else {
			endSmallRuns = source+length;
			insertionSort(startSmallRuns, endSmallRuns-startSmallRuns);
		}
	}

	if(topBytesSize%2) {
		std::memcpy(destination, source, sizeof(ELEM)*length);
		swap();
	}

#ifdef DEBUG
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cout << "Time Bottom: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif
#endif


	//We're done! Delete all the stuff we made
	delete [] destination;
	delete [] sourceBuckets;
	delete [] destinationBuckets;
	delete [] overflowBuckets;
	delete [] overflowBuffer;
	delete [] topBytes;
	delete [] bottomBytes;
	delete [] bucketCounts;
	delete [] overflowCounts;
	delete [] ladleBuffer;
	delete [] ladleBuckets;
}
#endif

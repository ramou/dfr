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
const int DIVERSION_THRESHOLD = 10;
#endif

const int SMALL_SAMPLE = 26;
const int LARGE_SAMPLE = 1024;	 // 1024 because I like large multiple of 256
const int DISTRIBUTION_SENSITIVE_THRESHOLD = 4096; 	/* 	If a length of data to be processed is smaller than this and
						      		we have no other data, we don't want to fart around figuring
								out the actual distribution, just assume uniform and deal with 
								the consequences. This can come up on on the first pass of the
								top bytes or possibly the first pass of the bottom bytes.
							*/

template <typename INT, typename ELEM>
void dfr(ELEM *source, auto length) {

	/**
		1) Diverting on small lengths
	*/

	auto insertionSort = [](ELEM *source, const auto &length) {
        	ELEM buf;
        	INT val;
        	for(unsigned i = 1; i < length; i++) {
			unsigned cursor = i;
                	buf = source[cursor];
                	val = *(reinterpret_cast<INT*>(source+(cursor=i)));
                	while(cursor > 0 &&
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
	ELEM **overflowBuckets = new ELEM*[256]; //Made bigger to avoid valgrind error :\

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

	auto foundLiveBits = 0;
	auto neededBytes = sizeof(INT);
	const auto allBits = 8*sizeof(INT);
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

	const auto passOverInput = [](ELEM* start, ELEM *end, const uint8_t &currentByte, const auto &deal, const auto &gatherStats) {
		if(start == end) return;
		const ELEM* element = start;
		const uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte; 
		const unsigned len = (end-start);

		for(unsigned i = 0; i < len; ++i, ++element, target+=sizeof(ELEM)) {
			deal(*target, *element);
			gatherStats(*target, *element);
		}
	};

        const auto passOverInputDealWithOverflowAndGatherLiveBits = 
		[&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize, &livebits, &bitmask]
		(ELEM* start, ELEM *end, const uint8_t &currentByte, const auto &destinationBuckets) {
                if(start == end) return;
                const ELEM* element = start;
                const uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const unsigned len = (end-start);

                for(unsigned i = 0; i < len; ++i, livebits |= bitmask ^ reinterpret_cast<INT>(*(element++)), target+=sizeof(ELEM)) {
			ELEM **currentDestination = destinationBuckets + ((*target) << 1);
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

        const auto passOverInputDealWithOverflow =
                [&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize]
                (ELEM* start, ELEM *end, const uint8_t &currentByte, const auto &destinationBuckets) {
                if(start == end) return;
                const ELEM* element = start;
                const uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
                const unsigned len = (end-start);

                for(unsigned i = 0; i < len; ++i, ++element, target+=sizeof(ELEM)) {
                        ELEM **currentDestination = destinationBuckets + ((*target) << 1);
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
		(ELEM *thisSource, ELEM *thisBuffer, const auto &currentByte, const auto &destinationBuckets) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(unsigned i = 0; i < 256; i++) {
                                        passOverInputDealWithOverflow(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], currentByte, destinationBuckets);
                                        passOverInputDealWithOverflow(startOverflowBucket, endOverflowBucket= overflowBuckets[i], currentByte, destinationBuckets);

                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };

        const auto passOverInputDealWithOverflowAndCounting =
                [&source, &overflow, &overflowCounts, &overflowBuffer, &overflowMaxSize, &bucketCounts]
                (ELEM* start, ELEM *end, const uint8_t &currentByte, const auto &destinationBuckets, const uint8_t &countByte) {
                if(start == end) return;
                const ELEM* element = start;
                const uint8_t *target = reinterpret_cast<const uint8_t*>(element) + currentByte;
		const uint8_t *countTarget = reinterpret_cast<const uint8_t*>(element) + countByte;

                const unsigned len = (end-start);

                for(unsigned i = 0; i < len; ++i, element++, target+=sizeof(ELEM), bucketCounts[*countTarget]++, countTarget+=sizeof(ELEM)) {
                        ELEM **currentDestination = destinationBuckets + ((*target) << 1);
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
                (ELEM *thisSource, ELEM *thisBuffer, const auto &currentByte, const auto &destinationBuckets, const uint8_t &countByte) {
                        auto startSourceBucket = thisSource;
                        ELEM* endSourceBucket;
                        auto startOverflowBucket = thisBuffer;
                        ELEM* endOverflowBucket;

                        for(unsigned i = 0; i < 256; i++) {
                                        passOverInputDealWithOverflowAndCounting(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], currentByte, destinationBuckets, countByte);
                                        passOverInputDealWithOverflowAndCounting(startOverflowBucket, endOverflowBucket= overflowBuckets[i], currentByte, destinationBuckets, countByte);

                                        startSourceBucket=sourceBuckets[(i<<1)+1];
                                        startOverflowBucket = endOverflowBucket;
                        }
        };


	const auto passOverInputs = [&](ELEM *thisSource, ELEM *thisBuffer, const auto &currentByte, const auto &deal, const auto &gatherStats) {
			auto startSourceBucket = thisSource;
			ELEM* endSourceBucket;
			auto startOverflowBucket = thisBuffer;
			ELEM* endOverflowBucket;

			for(unsigned i = 0; i < 256; i++) {
					passOverInput(startSourceBucket, endSourceBucket = sourceBuckets[i<<1], currentByte, deal, gatherStats);
					passOverInput(startOverflowBucket, endOverflowBucket= overflowBuckets[i], currentByte, deal, gatherStats);

					startSourceBucket=sourceBuckets[(i<<1)+1];
					startOverflowBucket = endOverflowBucket;
			}
	};

	const auto dealWithOverflow = [&](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[targetBucketIndex << 1];
		if(currentDestination < destinationBuckets[(targetBucketIndex << 1) + 1]) {
			*currentDestination = element;
			destinationBuckets[targetBucketIndex << 1]++;
		} else {
			overflowCounts[targetBucketIndex]++;
			*overflow = element;
			overflow++;
			if(overflow == overflowBuffer+(2 * overflowMaxSize)) {
				overflow = source;
			}
		}
	};

	#ifdef DEBUG
	const auto outputDeal = [&destinationBuckets](const auto &targetBucketIndex, const auto &element) {
		std::cout << (reinterpret_cast<INT>(element)) << " ";
        };
	#endif

	const auto simpleDeal = [&destinationBuckets](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[0];
                *currentDestination=element;
                destinationBuckets[0]++;
	};

	const auto dealExact = [&destinationBuckets, &destination](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[targetBucketIndex];
		*currentDestination=element;
		destinationBuckets[targetBucketIndex]++;
	};

	const auto dealToOverflow = [&overflowBuckets](auto &targetBucketIndex, const auto &element) {
		auto currentDestination = overflowBuckets[targetBucketIndex];
		*currentDestination=element;
		overflowBuckets[targetBucketIndex]++;
	};

	const auto noDeal = [](const auto &targetBucketIndex, const auto &element){};
	const auto noStats = [](const auto &targetBucketIndex, const auto &element){};

	const auto gatherLiveBits = [&livebits, &bitmask](const auto &targetBucketIndex, const auto &element){
		livebits |= bitmask ^ reinterpret_cast<INT>(element);
	};

	const auto countForByte = [&] (const auto &targetBucketIndex, const auto &element) {
		auto countedBucketIndex = ((reinterpret_cast<INT>(element))>>(countedByte*8)) & 255;
		bucketCounts[countedBucketIndex]++;
	};

	const auto convertCountsToContiguousBuckets = [&destination](const auto &counts, const auto &buckets, const auto &buffer) {
		*buckets=buffer;
		auto currentBucket = buckets+1;
		auto previousBucket = buckets+0;

		std::for_each(counts, counts+255, [&currentBucket, &previousBucket](auto &count) {
			*currentBucket = *previousBucket + count;
			currentBucket++;
			previousBucket++;
			count=0;
		});
	};

	const auto processOverflow = [&](const auto &currentByte){
		//Check if overflow is in source or if it still fits within existing buffer.
		if((source <= overflow) && ( overflow < (source+length))) {

			//We used up the old overflow so we need a bigger overflow.

			//Make a new buffer
			int newsize = overflowMaxSize + overflow-source;
			ELEM *oldOverflowBuffer = overflowBuffer;

			overflowBuffer = new ELEM[2*newsize];

			//build overflow buckets
			convertCountsToContiguousBuckets(overflowCounts, overflowBuckets, overflowBuffer);

			//pass left-to right from the overflow buffer, then from the front of the source also used as overflow buffer.
			passOverInput(oldOverflowBuffer+overflowMaxSize, oldOverflowBuffer+(2*overflowMaxSize), currentByte, dealToOverflow, noStats);
			passOverInput(source, overflow, currentByte, dealToOverflow, noStats);

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
				passOverInput(overflowBuffer+overflowMaxSize,overflow, currentByte, dealToOverflow, noStats);
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

	        if(neededBytes == 1) { // Just do a single count pass first and in a very un-fast-radix-way deal exactly into destination
			countedByte=currentByte;
			passOverInput(source, source+length, currentByte, noDeal, countForByte);
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

			if(hasByteLists) {
				passOverInput(source, source+length, currentByte, dealExact, noStats);
			} else {
				passOverInput(source, source+length, currentByte, dealExact, gatherLiveBits);
				buildByteLists(currentByte);
				hasByteLists=true;
			}
			swap();
	        } else if(neededBytes == 2) {   // We need a pass that does count capturing, then deals back in place.
						// We ignore when highest byte is dead; too annoying to check for
			countedByte=currentByte+1;
			doEstimates();
			passOverInput(source, source+length, currentByte, dealWithOverflow, countForByte);
			processOverflow(currentByte);
			swap();
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

			if(hasByteLists) {
				passOverInputs(source, overflowBuffer, countedByte, dealExact, noStats);
			} else {
				passOverInputs(source, overflowBuffer, countedByte, dealExact, gatherLiveBits);
				buildByteLists(currentByte);
				hasByteLists=true;
			}
			swap();
	        } else {
			//Do a first pass
#ifdef TIMINGS
        std::chrono::high_resolution_clock::time_point end;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
#endif

			doEstimates();
			if(hasByteLists) {
					passOverInput(source, source+length, currentByte, dealWithOverflow, noStats);
			} else {
					passOverInputDealWithOverflowAndGatherLiveBits(source, source+length, currentByte, destinationBuckets);
					buildByteLists(currentByte);
					hasByteLists=true;
					setByteListHere=true;
#ifdef DEBUG
std::cout << "Top Bytes: " << topBytesSize << " Bottom Bytes: " << bottomBytesSize << std::endl;
#endif

			}
			processOverflow(currentByte);
			swap();
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cerr << "First Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif



		//We have to account for the newly generated numBytes being either 0, 1, 2 or 3 (we did one pass, but numbytes still has it... 
		//our one pass could theoretically have been over dead bits), which means we found enough dead bits to eliminate
		//passes that would have got us caught up above in the numBytes == 1 or numBytes == 2 section. Sure it's annoying, but actually
		//doing the deadbit passes would be more annoying, so rejoice!

		if(setByteListHere && currentByte != bytes[0]) { //We did a pass that wasn't needed. Reverse and recurse
			destinationBuckets[0] = destination; //simpleDeal needs this
			passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
			swap();
			return true; //Tell the caller to re-call fr now that we're "fixed" things.
		}

		if (numBytes == 1) {	//We did a pass with overflow but we were only supposed to do one pass anyway
						// Since we didn't count we need to fix that
				//Ok, we did a pass, it was needed, but we should have made it an
				// exact pass. We have to fix that here.
			destinationBuckets[0] = destination; //simpleDeal needs this
			passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
			swap();
			std::memcpy(destination, source, length*sizeof(ELEM));
			swap();

		} else if (numBytes == 2) {	//We did one of two passes. 
			//We did a real pass, it just would have helped to have counted first.
			countedByte = bytes[1];
			passOverInputs(source, overflowBuffer, countedByte, noDeal, countForByte);
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);
			passOverInputs(source, overflowBuffer, countedByte, dealExact, noStats);
			swap();
		} else {

				for(int i = 1; i < (numBytes - 2); i++) {
#ifdef TIMINGS
        start = std::chrono::high_resolution_clock::now();
#endif

						doEstimates();
						passOverInputsDealWithOverflow(source, overflowBuffer, bytes[i], destinationBuckets);
						processOverflow(bytes[i]);
						swap();
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cerr << "Middle Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif

				}


				//Do last two passes
#ifdef TIMINGS
        start = std::chrono::high_resolution_clock::now();
#endif

				countedByte=bytes[numBytes-1];
				doEstimates();
				//passOverInputs(source, overflowBuffer, bytes[numBytes-2], dealWithOverflow, countForByte);
				passOverInputsDealWithOverflowAndCounting(source, overflowBuffer, bytes[numBytes-2], destinationBuckets, bytes[numBytes-1]);
				processOverflow(bytes[numBytes-2]);
				swap();
				convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cerr << "Next-To-Last Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
start = std::chrono::high_resolution_clock::now();
#endif


				//Deal last pass in top bytes
				passOverInputs(source, overflowBuffer, countedByte, dealExact, noStats);
				swap();
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cerr << "Last Pass Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
#endif

			}
		}

		return false;
	};


	//
	// Process  Top Bytes
	//

#ifdef TIMINGS
        std::chrono::high_resolution_clock::time_point end;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
#endif


	if(fr(topBytes, topBytesSize)) fr(topBytes, topBytesSize);
#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cerr << "Time Top: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
start = std::chrono::high_resolution_clock::now();
#endif

	if(bottomBytesSize > 0) {


		//At this point, source has our data in disjoint buckets and destination is ready.
		//buckets contains the start and end positions of each disjoint bucket in source.
		//overflowBuffer is the contiguous secondary source, with end buckets defined by
		//overflowBuckets, and is twice overflowMaxSize so we can deal into it to before
		//we make use of the input as a buffer source to prevent writing over a live bucket.

		// ((reinterpret_cast<INT>(element))>

		auto getHighBits = [&highBitMask](const ELEM &element) {
			return (reinterpret_cast<INT>(element)) & highBitMask;
		};

		INT currentBits = getHighBits(*source);

		bool potentialLongRun = false;
		bool definiteLongRun = false;
		ELEM* startSmallRuns;
		ELEM* endSmallRuns;
		ELEM* startDefiniteLongRun;
		ELEM* endDefiniteLongRun;

		startSmallRuns = source;

		for(ELEM* element = source+1; element < source+length; element += (DIVERSION_THRESHOLD>>1)) {
			const INT newBits = getHighBits(*element);
			if(currentBits ^ newBits) { //Things changed
				if(definiteLongRun) { // A long run has ended
					//walk back by one to find the exact end of the long run
					endDefiniteLongRun=element;
					while(!(getHighBits(*(endDefiniteLongRun-1)) ^ currentBits)) endDefiniteLongRun--;

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
					if(bottomBytesSize%2) {
						std::memcpy(destination, source, length*sizeof(ELEM));
					}

					source = oldSource;
					destination = oldDestination;
					length = oldLength;

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

#ifdef TIMINGS
end = std::chrono::high_resolution_clock::now();
std::cerr << "Time Bottom: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";
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
}
#endif

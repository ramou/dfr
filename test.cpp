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

/*
	Useful to make randomn test data
*/
template<typename UINT>
void make_random(UINT *data, unsigned N) {
        UINT max = std::numeric_limits<UINT>::max();
	std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<UINT> distribution(0,max);
        for(unsigned i = 0; i<N;i++)data[i]=distribution(generator);
}


auto estimated_pass = [] (
	auto &source, 
	auto &dest, 
	const auto length, 
	auto &currentByteCounter, 
	auto &overflowByteCounter, 
	auto bucketGuess, 
	const auto target_byte, 
	const auto next_byte) {

	unsigned overflowOffset=0;
	std::for_each(source, source+length, [&](auto &n) {
		//Get the bytes that we want to work with 
		const uint8_t target = *(reinterpret_cast<const uint8_t*>(&n)+target_byte);
		const uint8_t next = *(reinterpret_cast<const uint8_t*>(&n)+next_byte);

		if(currentByteCounter[target] < (target+1)*bucketGuess) {
			std::cout << "underflow" << std::endl;
			unsigned *countPtr = currentByteCounter + target;
			dest[*countPtr] = n;
			++(*countPtr);
		} else {
			std::cout << "overflow" << std::endl;
			source[overflowOffset++]=n;
			++overflowByteCounter[target];
		}
	});

	return overflowOffset;
};

template <typename INT, typename ELEM>
void insertionSort(ELEM *source, const auto &length) {
	ELEM buf;
	auto cursor=0;
	INT val;
	for(auto i = 1; i < length; i++) {
                buf = source[cursor = i];
		val = *(reinterpret_cast<INT*>(&buf));
                while(cursor > 0 &&
		      val < *(reinterpret_cast<INT*>(source + (cursor-1)))
		     ) {
                        source[cursor]=source[cursor-1];
                        cursor--;
                }
                source[cursor]=buf;
        }
}

//This should be ~14 lowering while I dev
const int DIVERSION_THRESHOLD = 4;

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

	if(length <= DIVERSION_THRESHOLD) {
		insertionSort<INT, ELEM>(source, length);
		return;
	}


	/**
		2) Allocate Memory
	*/

	//The number of bytes
	const unsigned int bytecount = sizeof(INT);

	//The buffer for doing radix passes
	ELEM *destination = new ELEM[length];
	//std::memset(destination, 0, length*sizeof(ELEM)); //Not sure I need this

	//These will start as pairs of positions pointing to the start of estimated
	//buckets. As values are placed, the second entry will increment. After all
	//data is dealt, this will be a set of start/end points for each bucket
	//I suggest all start values, then all end values as a cache-efficient format
	//When estimating uniform buckets, we can use math and the first 256 values
	//To make this fancier/faster.
	//
	//We need two sets of these, one for each buffer, and we'll swap them
	ELEM **sourceBuckets = new ELEM*[512];
	//std::memset(sourceBuckets, 0, 512*sizeof(ELEM*)); //Not sure I need this

	ELEM **destinationBuckets = new ELEM*[512];
	//std::memset(destinationBuckets, 0, 512*sizeof(ELEM*)); //Not sure I need this


	//Initially these will be the starting positions into the overflow buffer
	//Once dumped into the overflow buffer, they will then be the ending positions
	//when dealing out of the overflow buffer
	ELEM **overflowBuckets = new ELEM*[256];
	//std::memset(overflowBuckets, 0, 256*sizeof(ELEM*));//Not sure I need this

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

	/* Correctly random. Switch back to this when bugs resolved!
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	*/

	// Fixed check to force repeatable bugs
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + 1));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + 2));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + 3));

	auto neededBits = ceil(std::log2((float)(length)/(float)(DIVERSION_THRESHOLD))); //ceil(log(1792/14)/log(2)),

	//std::cout << "We need " << neededBits << " bits for diverting " << length << " with diversion threshold " << DIVERSION_THRESHOLD << std::endl;

	auto foundLiveBits = 0;
	auto neededBytes = sizeof(INT);
		#ifdef DEBUG
		std::cout << "When sampling we found the following live bits: " << std::endl << std::bitset<sizeof(INT)*8>(livebits) << std::endl;
		#endif

	const auto allBits = 8*sizeof(INT);

                #ifdef DEBUG
                std::cout << "We are working with a max of " << allBits << " bits." << std::endl;
                #endif


	for(auto i = 1; i < allBits; i++) {
		auto val = ((livebits >> (int)(allBits - i)) & 1);
		foundLiveBits += val;
		if(foundLiveBits==neededBits) {
				#ifdef DEBUG
				std::cout << "We processed " << i << " bits before finding enough so our top bytes are " << neededBytes << std::endl;
				#endif
			neededBytes = ceil((i+1)/8.0);
			break;
		}
	}
        //std::cout << "Given that sampled livebits are " << std::endl << std::bitset<sizeof(INT)*3>(livebits) << std::endl;
        //std::cout << "we will need " << neededBytes << " passes." << std::endl;

	int countedByte;
	const auto topBytes = new int[sizeof(INT)];
	std::memset(topBytes, 0, sizeof(INT)*sizeof(int));
	auto topBytesSize = 0;

	const auto bottomBytes = new int[sizeof(INT)];
	std::memset(bottomBytes, 0, sizeof(INT)*sizeof(int));
        auto bottomBytesSize = 0;

	const INT highBitMask = !((allBits-1) > neededBytes);
        auto currentByte = sizeof(INT)-neededBytes;
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

                #ifdef DEBUG
                std::cout << "HERE! Building estimates for destination " << destination  << std::endl;
		std::cout << "Bucket Guess is: " << bucketGuess << " and we have a remainder of " << remainder << std::endl;
                #endif


		auto buildBucketEstimates = [&destinationBuckets, &d] (const auto &bucketSize, const auto &start, const auto &end, const auto &destinationBuckets){
			for(auto i = start; i < end; i++) {
				destinationBuckets[i << 1]       = d;//start
				d += bucketSize;		//these buckets are one bigger than the latter buckets
				destinationBuckets[(i << 1) + 1] = d;//end, which will be start of next bucket
			}
		};

		buildBucketEstimates(bucketGuess+1, 0, remainder, destinationBuckets);
		buildBucketEstimates(bucketGuess, remainder, 256, destinationBuckets);

		#ifdef DEBUG
		for(int i = 0; i < 256; i++) {
			std::cout << " ob " << i << " " << (destinationBuckets[i << 1]-destination) << " " << (destinationBuckets[(i << 1) + 1]-destination) << std::endl;
		}

		#endif

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
	};

	auto swap = [&source, &destination, &sourceBuckets, &destinationBuckets, &length] {
                        #ifdef DEBUG
                        for(int i = 0; i < length; i++) {
                                std::cout << destination[i] << " ";
                        }
                        std::cout << std::endl;
                        #endif

		std::swap(source, destination);
		std::swap(sourceBuckets, destinationBuckets);
	};

	auto passOverInput = [](const auto &start, const auto &end, const auto &currentByte, const auto &deal, const auto &gatherStats) {
		std::for_each(start, end, [&](auto &element){
			auto targetBucketIndex = (((reinterpret_cast<INT>(element))>>(currentByte*8)) & 255);
			deal(targetBucketIndex, element);
			gatherStats(targetBucketIndex, element);
		});
	};

	auto passOverInputs = [&](const auto thisSource, const auto thisBuffer, const auto &currentByte, const auto &deal, const auto &gatherStats) {
			auto startSourceBucket = thisSource;
			ELEM* endSourceBucket;
			auto startOverflowBucket = thisBuffer;
			ELEM* endOverflowBucket;

                #ifdef DEBUG
                //std::cout << "HERE! Passing over inputs with source " << source << " and destination " << destination  << std::endl;
                #endif


			for(int i = 0; i < 256; i++) {
					endSourceBucket = sourceBuckets[i<<1];
					endOverflowBucket = overflowBuckets[i];


                #ifdef DEBUG
		if(startSourceBucket!=endSourceBucket)
 	        //       std::cout << "HERE! Passing over source bucket " << i << " from position " << (startSourceBucket-source) << " to " << (endSourceBucket-source)  << std::endl;
                #endif

					passOverInput(startSourceBucket, endSourceBucket, currentByte, deal, gatherStats);
                #ifdef DEBUG
		if(startOverflowBucket != endOverflowBucket)
                //	std::cout << "HERE! Passing over overflow bucket " << i << " from position " << (startOverflowBucket-overflowBuffer) << " to " << (endOverflowBucket-overflowBuffer)  << std::endl;
                #endif

					passOverInput(startOverflowBucket, endOverflowBucket, currentByte, deal, gatherStats);

					startSourceBucket=sourceBuckets[(i<<1)+1];
					startOverflowBucket = endOverflowBucket;
			}
	};

	auto dealWithOverflow = [&](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[targetBucketIndex << 1];
		#ifdef DEBUG
			std::cout << "Attempting to deal " << targetBucketIndex << "(" << (reinterpret_cast<INT>(element)) << ") into: " << (currentDestination-destination) << std::endl; 
		#endif
		if(currentDestination < destinationBuckets[(targetBucketIndex << 1) + 1]) {
			#ifdef DEBUG
			std::cout << "\tDeal Successful" << std::endl;
			#endif
			*currentDestination = element;
			destinationBuckets[targetBucketIndex << 1]++;
		} else {
			#ifdef DEBUG
			std::cout << "\tDeal Failed, dealing to overflow position: " << overflow << std::endl;
			#endif
			overflowCounts[targetBucketIndex]++;
			*overflow = element;
			overflow++;

			if(overflow == overflowBuffer+(2 * overflowMaxSize)) {
				overflow = source;
			}
		}
	};

	auto simpleDeal = [&destinationBuckets](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[0];
                *currentDestination=element;
                destinationBuckets[0]++;
	};

	#ifdef DEBUG
	auto outputDeal = [&destinationBuckets](const auto &targetBucketIndex, const auto &element) {
		std::cout << (reinterpret_cast<INT>(element)) << " ";
        };
	#endif


	auto dealExact = [&destinationBuckets, &destination](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[targetBucketIndex];
                #ifdef DEBUG
                        std::cout << "Attempting to deal " << targetBucketIndex << "(" << (reinterpret_cast<INT>(element)) << ") into: " << (currentDestination-destination) << std::endl;
                #endif
		*currentDestination=element;
		destinationBuckets[targetBucketIndex<<1]++;
	};

	auto dealToOverflow = [&overflowBuckets](auto &targetBucketIndex, const auto &element) {
		auto currentDestination = overflowBuckets[targetBucketIndex];
		*currentDestination=element;
		overflowBuckets[targetBucketIndex]++;
	};

	auto noDeal = [](const auto &targetBucketIndex, const auto &element){};
	auto noStats = [](const auto &targetBucketIndex, const auto &element){};

	auto gatherLiveBits = [&livebits, &bitmask](const auto &targetBucketIndex, const auto &element){
		livebits |= bitmask ^ reinterpret_cast<INT>(element);
	};

	auto countForByte = [&] (const auto &targetBucketIndex, const auto &element) {
		auto countedBucketIndex = ((reinterpret_cast<INT>(element))>>(countedByte*8)) & 255;
		bucketCounts[countedBucketIndex]++;
	};

	auto convertCountsToContiguousBuckets = [&destination](const auto &counts, const auto &buckets, const auto &buffer) {
		*buckets=buffer;
		auto currentBucket = buckets+1;
		auto previousBucket = buckets+0;

		int i = 0;

		std::for_each(counts, counts+256, [&currentBucket, &previousBucket, &i, &destination](auto &count) {
                #ifdef DEBUG
			std::cout << "Bucket " << i <<  " will have " << count << " values and will start at position " << (((*previousBucket))-destination) << std::endl;
			i++;
                #endif
			*currentBucket = *previousBucket + count;
			currentBucket++;
			previousBucket++;
			count=0;
		});
	};

	auto processOverflow = [&](const auto &currentByte){

                #ifdef DEBUG
                std::cout << "Processing Overflow At Byte " <<  currentByte  << std::endl;
                #endif


                #ifdef DEBUG
                std::cout << "\tSource: " <<  source  << std::endl;
		std::cout << "\tOverflowBuffer: " <<  overflowBuffer  << std::endl;
		std::cout << "\tOverflow: " <<  overflow  << std::endl;
                #endif


		//Check if overflow is in source or if it still fits within existing buffer.
		if((source <= overflow) && ( overflow < (source+length))) {
			//We used up the old overflow so we need a bigger overflow.

			//Make a new buffer
			int newsize = overflowMaxSize + overflow-source;
			ELEM *oldOverflowBuffer = overflowBuffer;

                #ifdef DEBUG
                std::cout << "\tAllocating overflow buffer of size " << 2*newsize << std::endl;
                #endif

			overflowBuffer = new ELEM[2*newsize];

		#ifdef DEBUG
                std::cout << "\tCreated overflow buffer at " << overflowBuffer << std::endl;
                #endif


                #ifdef DEBUG
                std::cout << "\tConverting overflowCounts to overflowBuckets" << std::endl;
                #endif


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
			#ifdef DEBUG
	                std::cout << "killing overflow at " << oldOverflowBuffer << std::endl;
        	        #endif
			if(oldOverflowBuffer) {
                #ifdef DEBUG
                std::cout << "\t We delete the overflowBuffer at " << oldOverflowBuffer << std::endl;
                #endif

				delete [] oldOverflowBuffer;
				//free(oldOverflowBuffer);
			} else {
                #ifdef DEBUG
                std::cout << "\t No need to delete overflow that wasn't allocated" << std::endl;
                #endif

			}
		} else {
                #ifdef DEBUG
                std::cout << "\tOverflow is pointing in overflow buffer: " <<  overflow  << std::endl;
                #endif
			if(overflowMaxSize == 0) return; //We have no overflow, overflow is probably null
			convertCountsToContiguousBuckets(overflowCounts, overflowBuckets, overflowBuffer);
			passOverInput(overflowBuffer+overflowMaxSize,overflow, currentByte, dealToOverflow, noStats);
		}
	};

        // Allocates bytes, in order, to either the list of Top Bytes
        // or Bottom Bytes, based on the Threshold Byte which doesn't
        // get added to either because presumably that byte was used
        // directly or indirectly to scan for live bits.
        //
        // This is so cheap let's not do anything fancy
	auto buildByteLists = [&](const auto &threshold) {
			#ifdef DEBUG
			std::cout << "building byte list around threshold " << threshold<<  std::endl;
                	std::cout << "\twith these bits: " << std::bitset<sizeof(INT)*8>(livebits) << std::endl;
                        std::cout << "\ttopBytesSize old: " << topBytesSize <<  std::endl;
                        std::cout << "\tbottomBytesSize old: " << bottomBytesSize <<  std::endl;
			#endif
		for(auto i = 0; i < bytecount; i++) {
			if((livebits >> (i*8)) & 255) {
				if(i < threshold) {
					bottomBytes[bottomBytesSize] = i;
					bottomBytesSize++;
				} else if(i >= threshold) {
			#ifdef DEBUG
                        std::cout << "found live byte " << ((livebits >> (i*8)) & 255) << " at " << i <<  std::endl;
                        #endif

					topBytes[topBytesSize] = i;
					topBytesSize++;
				}
			}
		}
			#ifdef DEBUG
                        std::cout << "\ttopBytesSize new: " << topBytesSize <<  std::endl;
                        std::cout << "\tbottomBytesSize new: " << bottomBytesSize <<  std::endl;
			#endif
	};

	bool hasByteLists = false;

 	auto fr = [&](const auto &bytes, const auto &numBytes) {

		const auto thresholdByte = sizeof(INT)-neededBytes;

		#ifdef DEBUG
                std::cout << "HERE! Running FR with neededBytes: " << neededBytes  << " and threshold byte " << thresholdByte << " and numBytes " << numBytes << std::endl;
		std::cout << "But numBytes might not be reliable" << std::endl;
		std::cout << "topBytes: " << topBytes << std::endl;
		std::cout << "bottomBytes: " << bottomBytes << std::endl;
		std::cout << "bytes: " << bytes << std::endl;
                #endif

	        if(neededBytes == 1) { // Just do a single count pass first and in a very un-fast-radix-way deal exactly into destination
			countedByte=currentByte;
			passOverInput(source, source+length, thresholdByte, noDeal, countForByte);
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

			if(hasByteLists) {
				passOverInput(source, source+length, thresholdByte, dealExact, noStats);
			} else {
				passOverInput(source, source+length, thresholdByte, dealExact, gatherLiveBits);
				buildByteLists(thresholdByte);
				hasByteLists=true;
			}
	        } else if(neededBytes == 2) {   // We need a pass that does count capturing, then deals back in place.
										// We ignore when highest byte is dead; too annoying to check for
			countedByte=currentByte+1;
			doEstimates();
			passOverInput(source, source+length, thresholdByte, dealWithOverflow, countForByte);
			processOverflow(thresholdByte);
			swap();
			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

			if(hasByteLists) {
				passOverInputs(source, overflowBuffer, thresholdByte-1, dealExact, noStats);
			} else {
				passOverInputs(source, overflowBuffer, thresholdByte-1, dealExact, gatherLiveBits);
				buildByteLists(thresholdByte);
				hasByteLists=true;
			}
			swap();
	        } else {


                #ifdef DEBUG
                std::cout << "HERE! Passing over first input with source " << source << " and destination " << destination  << std::endl;
                #endif


			//Do a first pass
			doEstimates();
			if(hasByteLists) {
					currentByte = thresholdByte;
					passOverInput(source, source+length, currentByte, dealWithOverflow, noStats);
			} else {
					currentByte = bytes[0];
					passOverInput(source, source+length, currentByte, dealWithOverflow, gatherLiveBits);
					buildByteLists(thresholdByte);
					hasByteLists=true;
			}
			processOverflow(currentByte);
			swap();

			#ifdef DEBUG
                        std::cout << "Input from both sources normalized: " << std::endl;
                        passOverInputs(source, overflowBuffer, countedByte, outputDeal, noStats);
                        std::cout << std::endl;
                        #endif


                #ifdef DEBUG
		if(hasByteLists){
                	std::cout << "HERE! After the first pass we may have built the bytelists and now numBytes is " << numBytes << std::endl;
			std::cout << "\ttopBytesSize: " << topBytesSize <<  std::endl;
			std::cout << "\tbottomBytesSize: " << bottomBytesSize <<  std::endl;
		}
                #endif

		//We have to account for the newly generated numBytes being either 0, 1, 2 or 3 (we did one pass, but numbytes still has it... 
		//our one pass could theoretically have been over dead bits), which means we found enough dead bits to eliminate
		//passes that would have got us caught up above in the numBytes == 1 or numBytes == 2 section. Sure it's annoying, but actually
		//doing the deadbit passes would be more annoying, so rejoice!

		if(topBytesSize == 0) { //We did a pass with overflow, but we shouldn't have done any passes.
			destinationBuckets[0] = destination;
			passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
			swap();
		} else if (topBytesSize == 1) {	//We did a pass with overflow but we were only supposed to do one pass anyway
							// Since we didn't count we need to fix that
							// Further, if that one pass was the lowest-order bit, we may have mistakenly
							// Done a dead pass and need to account for that

                #ifdef DEBUG
                if(hasByteLists){
                        std::cout << "topBytesSize == 1 so we need to fix our first pass. " << std::endl;
                }
                #endif


			if(currentByte != bytes[0]) {//This means currentByte==0, otherwise it wouldn't have been in the initial 
							// bytes/threshold byte. basically did a bad pass and just need to swap back
							// The good news is that we shouldn't have done that pass so we just need a
							// simple deal. Since this simple deal gets enough info to fix stuff, let's
							// also do the deal we should have done.

                #ifdef DEBUG
                if(hasByteLists){
                        std::cout << "\t our first pass was the lowest order bit, but it was a dead pass." << std::endl;
                }
                #endif


				destinationBuckets[0] = destination;
				countedByte = bytes[0];
	                        passOverInputs(source, overflowBuffer, currentByte, simpleDeal, countForByte);
        	                swap();
				convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);
				passOverInput(source, source+length, bytes[0], dealExact, noStats);
				swap();
			} else {			//Ok, we did a pass, it was needed, but we should have made it an
							// exact pass. We have to fix that here.

                #ifdef DEBUG
                if(hasByteLists){
                        std::cout << "\t our first pass was a live pass, so we just need to fix it a bit." << std::endl;
                }
                #endif

				destinationBuckets[0] = destination;
				passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
				swap();
				std::memcpy(destination, source, length*sizeof(ELEM));
				swap();
			}

		} else if (numBytes == 2) {	//We did one of two passes. If it was the lowest order pass it may
							// have been a dead pass.
			if(currentByte != bytes[0]) {//We did a useless pass and have to undo that.
				destinationBuckets[0] = destination;
				passOverInputs(source, overflowBuffer, currentByte, simpleDeal, noStats);
				swap();
				doEstimates();
				countedByte = bytes[1];
				passOverInput(source, source+length, thresholdByte, dealWithOverflow, countForByte);
	                        processOverflow(thresholdByte);
        	                swap();
                                convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);
				passOverInputs(source, overflowBuffer, bytes[1], dealExact, noStats);
				swap();
			} else {	//We did a real pass, it just would have helped to have counted first.
				countedByte = bytes[1];
				passOverInputs(source, overflowBuffer, currentByte, noDeal, countForByte);
				convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);
				passOverInputs(source, overflowBuffer, currentByte, dealExact, noStats);
				swap();
			}

		} else {


                	#ifdef DEBUG
                	std::cout << "HERE! Processing Middle Runs: " << std::endl;
                	#endif

				for(int i = 1; i < (numBytes - 2); i++) {
						doEstimates();
                	#ifdef DEBUG
                	std::cout << "HERE! Passing over Inputs in pass " << i << std::endl;
                	#endif

						passOverInputs(source, overflowBuffer, bytes[i], dealWithOverflow, noStats);
                	#ifdef DEBUG
                	std::cout << "HERE! Processing Overflows in pass " << i  << std::endl;
                	#endif
						processOverflow(bytes[i]);
						swap();
				}

                	#ifdef DEBUG
                	std::cout << "Processing Last Two Runs for bytes " << bytes[numBytes-1] << " and " << bytes[numBytes-2] << std::endl;
                	#endif


				//Do last two passes
				countedByte=bytes[numBytes-1];
				doEstimates();
				passOverInputs(source, overflowBuffer, bytes[numBytes-2], dealWithOverflow, countForByte);
				processOverflow(bytes[numBytes-2]);
				swap();

				convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);


			#ifdef DEBUG
			std::cout << "Input from both sources normalized: " << std::endl;
			passOverInputs(source, overflowBuffer, countedByte, outputDeal, noStats);
			std::cout << std::endl;
			#endif

				//Deal last pass in top bytes
				passOverInputs(source, overflowBuffer, countedByte, dealExact, noStats);
				swap();

                        #ifdef DEBUG
                        std::cout << "Input from source: " << std::endl;
                        passOverInput(source, source+length, countedByte, outputDeal, noStats);
                        std::cout << std::endl;
                        #endif

			}
		}
	};


	//
	// Process  Top Bytes
	//

	fr(topBytes, topBytesSize);

                        #ifdef DEBUG
                        std::cout << "HERE2! We processed our " << topBytesSize << " top Bytes. Now to process the " << bottomBytesSize << " bottom bytes." << std::endl;
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

		for(ELEM* element = source+1; element < source+length; element += (DIVERSION_THRESHOLD>1)) {
			const INT newBits = getHighBits(*element);
			if(currentBits ^ highBitMask) { //Things changed
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

                        	#ifdef DEBUG
                        	std::cout << "We found a run of length " << length << ". Running bottom bytes." << std::endl;
                        	#endif


					fr(bottomBytes, bottomBytesSize);
					if((bottomBytesSize+topBytesSize)%2) {
						const int offset = startDefiniteLongRun-source;
						std::memcpy(destination+offset, source+offset, (endDefiniteLongRun-startDefiniteLongRun)*sizeof(ELEM));
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
							insertionSort<INT, ELEM>(startSmallRuns, endSmallRuns-startSmallRuns);
							if(topBytesSize%2) {
		                                                const int offset = startSmallRuns-source;
                		                                std::memcpy(destination+offset, source+offset, (endSmallRuns-startSmallRuns)*sizeof(ELEM));
							}
						}

						definiteLongRun=true;
						potentialLongRun=false;
					} else potentialLongRun = true;
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

                	fr(bottomBytes, bottomBytesSize);
                	if((bottomBytesSize+topBytesSize)%2) {
				const int offset = startDefiniteLongRun-source;
				std::memcpy(destination+offset, source+offset, (endDefiniteLongRun-startDefiniteLongRun)*sizeof(ELEM));
			}
                	source = oldSource;
                	destination = oldDestination;
                	length = oldLength;
		} else {
			endSmallRuns = source+length;
			insertionSort<INT, ELEM>(startSmallRuns, endSmallRuns-startSmallRuns);
			if(topBytesSize%2) {
				const int offset = startSmallRuns-source;
				std::memcpy(destination+offset, source+offset, (endSmallRuns-startSmallRuns)*sizeof(ELEM));
			}
		}
	}

	if((bottomBytesSize+topBytesSize)%2) {
		#ifdef DEBUG
		std::cout << "Do an extra swap back to the original array." << std::endl;
		#endif
		std::memcpy(destination, source, sizeof(ELEM)*length);
		swap();
	}


                        #ifdef DEBUG
                        std::cout << "HERE3! We're done, just cleaning up now." << std::endl;
			std::cout << "source is: " << source << std::endl;
			std::cout << "destination is: " << destination << std::endl;
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


template <typename T>
struct elem {
        T key;
        T *val;
};


typedef uint64_t targetType;

int main(int argc, char *argv[]) {

	auto targetLength = 512;

	if(argc >= 2) {
		targetLength = atoi(argv[1]);
	}

	targetType* values = new targetType[targetLength];

	if(argc > 2) {
		std::cout << "Making fixed input of size " << targetLength << " with values: " << std::flush;
		for(int i = 0; i < targetLength; i++) {
		  auto val = atoi(argv[i+2]);
		  std::cout << val << " " << std::flush;
		  values[i] = val;
		}
		std::cout << std::endl;
	} else {
		std::cout << "Making random input of size " << targetLength << std::endl;
		make_random(values, targetLength);
	}

	dfr<targetType, targetType>(values, targetLength);

	for(int i = 1; i < targetLength; i++) {
		if(values[i] < values[i-1]) {
			std::cout << " value " << i << " is out of place." << std::endl;
			break;
		}
	}

                        #ifdef DEBUG
			for(int i = 0; i < targetLength; i++) {
				std::cout << values[i] << " ";
			}
			std::cout << std::endl;
                        #endif


	delete [] values;

        return 0;
}


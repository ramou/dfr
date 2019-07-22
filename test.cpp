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

#define BYTE_SOURCE reinterpret_cast<const uint8_t*>(Source)

//Let's assume we allow this much overflow
const unsigned _OVERFLOW = 1;


/*
auto pass = [] (const auto start, const auto length, const auto byte_offsets, const auto number_of_offsets) {
	std::for_each(start, start+length, [&byte_offsets, &number_of_offsets](auto &n){
		std::cout << "Processing Number " << n << std::endl;
		int counter = 0;
		std::for_each(byte_offsets, byte_offsets+number_of_offsets, [&counter,&n](auto &t){
			uint8_t val = *(reinterpret_cast<const uint8_t*>(&n)+t);
			std::cout << "\t \t at byte " << t << ":" ;
			std::cout << " (" << std::setfill('0') << std::setw(3) << std::dec << (+val);
                        std::cout << ") " << std::bitset<8>(val);
                        std::cout << std::endl;
		});
	});
};
*/

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

	/*
		std::cout << "\t \t at byte " << target_byte << ":" ;
		std::cout << " (" << std::setfill('0') << std::setw(3) << std::dec << (+target);
		std::cout << ") " << std::bitset<8>(target);

		std::cout << std::endl;

		std::cout << "\t \t at byte " << next_byte << ":" ;
                std::cout << " (" << std::setfill('0') << std::setw(3) << std::dec << (+next);
                std::cout << ") " << std::bitset<8>(next);

		std::cout << std::endl;
	*/


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

/*
auto estimated_pass_two_source = [] (const auto source1, const auto source2, const auto length, const auto dest, const auto target_byte, const auto next_byte) {


};


auto counting_estimated_pass = (const auto source, const auto length, const auto dest, const auto target_byte, const auto next_byte) {


};

auto pass = (const auto source, const auto length, const auto dest, const auto counts, const auto target_byte, const auto next_byte) {


};
*/

template <typename INT, typename ELEM>
void insertionSort(ELEM *source, auto length) {
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

/*
	Let's assume that we're going to get data in the form of :
	struct ELEM {
		INT key,
		.. anything else
	}

	This storts an array of ELEMs as defined above.

	We will do a very quick scan for potentially active bits and choose enough bytes $t$ such that from highest order
	to lowest order byte we've accumulated enough. The minimum number of bits we'll try to accumulate is:
	\ceil{\log_{2}\frac{legnth}{DIVERSION_THRESHOLD}}

	For example, if after our active bit scan we see:
	7        6        5        4        3        2        1        0
	01110011 11111111 11100111 11111111 11101111 11111011 11111011 11111011
	where the bytes are labeled highest-order to lowest (7-0), and 1s  represent a known active bit and 0 represents not knowning 
	(or a lowere chance of it being active), and we had a length of 1792 and we used DIVERSION_THRESHOLD = 14, we would calculate 
	ceil(log(1792/14)/log(2)), and thus we'd try to find 7 active bits. In the above example, the highest order byte (7) is 
	01110011, which only has 5 identified active bits. As such, we would take the top two bytes as the top bytes, 7 and 6 for the
	first passes.

	Setting aside the determination of bucket counts, after these passes are done, the data will be sorted by the two highest order
	bytes only. The algorithm subsequently scans through the data to find large runs, diverting for small runs. If we divert to 
	Insertion Sort, we can scan with a skip of DIVERSION_THRESHOLD/2, looking only at the sorted bits. If they change, then we have
	skipped past one or more buckets and we can keep moving. If the bits don't change twice in a row, then we can continue scanning 
	forward until they to, then scan one-at-a-time on either end to find the boundaries of our big bucket. All small buckets before 
	this can be sorted and moved into the original array, and we are done processing that data (because the small buckets are 
	relatively ordered, Insertion Sort is suitable to sort them all at once). Large buckets will be processed using Fast Radix on the
	remaining bytes... since the top bytes are identical, they can be ignored.

	Note that the aside from the per-element cost (reduction) of Fast Radix, that algorithm also removes some of the counting 
	calculations that would otherwise be repeated for each large bucket found, thus the actual cost of performing repeated Fast Radix
	subsorts is not so noticeable (we argue that its impact is so low that these two improvements actually dovetail nicely).

	Once the big bucket is sorted, we continue to skip through the data in the same way till we get to the end, at which point we
	Insertion Sort any pending small buckets. 

	Let's look at an example using only two bytes, with the high-order byte being the top byte that things are sorted by. For the
	sake of simplicity, let's say DIVERSION_THRESHOLD=6.

00 (00)s  01011100	00101001 -get are starting bits for comparison
01        01100000	10110000
02        01100010	10011010
03 (01)c! 01110100	01000010 -our first check shows bits changed, so there are small buckets to deal with later
04 (09)c! 01110101	10100011
05 (10)c! 01110110	11110011
06 (02)c! 10011010	00101111 -02) Our second check also shows a change in bits 
   (11)c=                        -11) We found our match, this must be the beginning!
07        10011010	10101010
08        10011010	01000011
09 (03)c= 10011010	00001000 -03)We found matching bits once, so maybe we have a big bucket
10        10011010	00010011
11        10011010	01110100
12 (04)c= 10011010	10101001 -04)We found another set of contiguous matching bits, so we definitely have a big bucket! 
   (08)c=                        -08) found the match again, this is the end of the long bucket, time to scan for the beginning
13 (07)c! 11000011	00011110
14 (06)c! 11011011	10110011
15 (05)c! 11100011	00001110 -These bits are different, so the bucket must have ended earlier, scan back
16        11101001	11011011
17        11101001	01000110
18 (12)c! 11110010	11110010
19        11110011	01110110

The numbers in ()s represent the step, starting with (00) where we look at the bits in the first entry as our starting point. With 20
values we actually only check bits 12 times, with half those checks coming from scanning for the true ends of the large bucket, this
being the worst case. However, as the ratio of large to small buckets reduces, this should tend towards 2/DIVERSION_THRESHOLD*length
checks.

Insertion Sort will be called for 00 to 05. Fast Radix will be called on the remaining byte for 06 to 12. Insertion Sort will be called
again for 13 to 19. Note that Insertion Sort on many buckets of size 1 completes after n-1 comparisons, and is thus very very fast. The 
first Insertion Sort yeilds the first 6 values sorted in 5 comparisons and the second Insertion Sort takes 7 comparisons (one insertion)

If instead of 1 bottom byte there were more, Insertion Sort would perform just as fast, but the Fast Radix component would have to act
on each byte in urn. However, it will nearly always take the full amount of time on only a small subset of the data if diversion takes 
place. In the event that our top passes led to the majority of the data being in a few large runs, then this would run in only a small 
constant amount of time longer than had we just run Fast Radix on all the bytes at once.

	1) For arrays smaller than DIVERSION_THRESHOLD, we just insertion sort
	2) Allocate neceesary memory
		2a) We will create a buffer of the same size
		2b) We will create deal_indices for holding where to place elements
		2c) We will create overflow_indices for holding where to place overflow elements
	3) We will sample a very small number of elements (1+3) to estimate active bits.
	   This will let us choose the correct number of top bytes to process before
	   diversion.
	4) If length is less than DISTRIBUTION_SENSITIVE_THRESHOLD
		4a) we will assume a uniform distribution.
		4b) Otherwise, we will sample SMALL_SAMPLE values
			4b.a) If the sample suggests uniform, we will prepare uniformly sized buckets
			4b.b) Otherwise, we will sample up to LARGE_SAMPLE values (LARGE_SAMPLE-SMALL_SAMPLE) more samples 
			     and scale the results to the length, making those the estimated bucket sizes
	5) During the first pass of the top bytes (the least significant of the most significant), we will scan difinitively for 
	   active bits.
	6) Each subsequent pass will be on the next most significant byte unless it was determined to be completely inactive by the
	   bit scan. The sample step from 4 will be repeated, save that uniform distributions will made based in valid bits for that
	   byte when length is less than DISTRIBUTION_SENSITIVE_THRESHOLD and a large sample will immediately be checked for in larger
           sized inputs of the bit scan showed any invalid bits... arguably we could check for uniformity across the remaining bits, but
	   I think that's not worth the minor cost.
	7) When the top bytes have been sorted, we scan the input DiversionThreshold/2 at a time, checking for changes in the already 
	   sorted bits until the end of the data.
		7a) If a big bucket is detected, find its exact boundary
			7a.a) run insertion sort from the last sorted data to the beginning of the big bucket
				7a.a1) If in the buffer, memcopy back to original array
			7a.b) run Fast Radix on the remaining bytes in big bucket, qualifying that passes run like in 4, but perhaps 
			      we will do a second bit scan and use that instead of the previous bit scan since it will more accurately 
			      reflect this bucket on the remaining bytes.
				7c.a) If in the buffer, memcopy back to source
		7b) If the end is detected, do 7a.a, but from the last sorted data until the end.

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

	//If we ever actually count elements upon placing (e.g. a last pass), we should
	//use this before converting it to *buckets*
	unsigned* bucketCounts = new unsigned[256];
        std::memset(bucketCounts, 0, sizeof(bucketCounts));

	//We actually count for overflow, so we must use this. We'll just convert to 
	//overflowBuckets when done because that's easier to process.
	unsigned* overflowCounts = new unsigned[256];
        std::memset(overflowCounts, 0, sizeof(overflowCounts));


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
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));
	livebits |= bitmask^ *(reinterpret_cast<INT*>(source + dist(gen)));


	auto neededBits = ceil(std::log2((float)(length)/(float)(DIVERSION_THRESHOLD))); //ceil(log(1792/14)/log(2)),

	//std::cout << "We need " << neededBits << " bits for diverting " << length << " with diversion threshold " << DIVERSION_THRESHOLD << std::endl;

	auto foundLiveBits = 0;
	auto neededBytes = sizeof(INT);
		#ifdef DEBUG
		std::cout << "When sampling we found the following live bits: " << std::endl << std::bitset<sizeof(INT)*8>(livebits) << std::endl;
		#endif

	const auto allBits = pow(2, sizeof(INT));
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
	auto topBytes = new int[sizeof(INT)];
	auto topBytesSize = 0;

	auto bottomBytes = new int[sizeof(INT)];
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

		/*
		//checking for SMALL_SAMPLE... this is more efficient than
		unsigned int smallCounts[sizeof(INT)][SMALL_SAMPLE];
		//storing samples
		unsigned int smallCountSamples[sizeof(INT)][256] ;

		unsigned int step = length/SMALL_SAMPLE;

		// This doesn't work for small data!
		// Pointers to sample values
		ELEM *smallSamples[SMALL_SAMPLE];

		for(int i = 1; i < SMALL_SAMPLE; i++) {
				smallSamples[i]=&source[i*step+dist(gen)];
		}
		*/
		}

		#ifdef DEBUG
		std::cout << "These are the bucket starts and ends:" << std::endl;
		for(auto i = 0; i < 256; i++) {
			std::cout << "\t" << i << "\t" << (destinationBuckets[i<<1]-destination) << " to " <<  (destinationBuckets[(i<<1) + 1]-destination) << std::endl;
		}
		#endif
	};

	auto copyBuffer = [&destination] (const auto start, const auto end) {
		std::memcpy(destination + (end-start), start, end-start);
	};

	auto swap = [&source, &destination, &sourceBuckets, &destinationBuckets] {
		std::swap(source, destination);
		std::swap(sourceBuckets, destinationBuckets);
	};

	auto passOverInput = [](const auto &start, const auto &end, const auto &currentByte, const auto &deal, const auto &gatherStats) {
		std::for_each(start, end, [&](auto &element){
			auto targetBucketIndex = (((reinterpret_cast<INT>(element))>>currentByte) & 255);
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
                std::cout << "HERE! Passing over inputs with source " << source << " and destination " << destination  << std::endl;
                #endif


			for(int i = 0; i < 256; i++) {
					endSourceBucket = sourceBuckets[i<<1];
					endOverflowBucket = overflowBuckets[i];


                #ifdef DEBUG
		if(startSourceBucket!=endSourceBucket)
 	               std::cout << "HERE! Passing over source bucket " << i << " from position " << (startSourceBucket-source) << " to " << (endSourceBucket-source)  << std::endl;
                #endif

					passOverInput(startSourceBucket, endSourceBucket, currentByte, deal, gatherStats);
                #ifdef DEBUG
		if(startOverflowBucket != endOverflowBucket)
                	std::cout << "HERE! Passing over overflow bucket " << i << " from position " << (startOverflowBucket-overflowBuffer) << " to " << (endOverflowBucket-overflowBuffer)  << std::endl;
                #endif

					passOverInput(startOverflowBucket, endOverflowBucket, currentByte, deal, gatherStats);

					startSourceBucket=sourceBuckets[(i<<1)+1];
					startOverflowBucket = endOverflowBucket;
			}
	};

	auto dealWithOverflow = [&](const auto &targetBucketIndex, const auto &element) {
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

	auto dealExact = [&destinationBuckets](const auto &targetBucketIndex, const auto &element) {
		auto currentDestination = destinationBuckets[targetBucketIndex<<1];
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
		auto countedBucketIndex = ((reinterpret_cast<INT>(element))>>countedByte) & 255;
		bucketCounts[countedBucketIndex]++;
	};

	auto convertCountsToContiguousBuckets = [](const auto counts, const auto buckets, const auto buffer) {
		buckets[0]=buffer;
		auto currentBucket = &buckets[1];
		auto previousBucket = &buckets[0];

		std::for_each(counts, counts+256, [&currentBucket, &previousBucket](auto &count) {
			currentBucket[0]=previousBucket[0]+count;
			currentBucket++;
			previousBucket++;
		});
		std::memset(counts, 0, (sizeof(counts)*256));

	};

	auto processOverflow = [&](const auto &currentByte){

                #ifdef DEBUG
                std::cout << "HERE! " <<  currentByte  << std::endl;
                #endif


		//Check if overflow is in source or if it still fits within existing buffer.
		if((source <= overflow) && ( overflow < (source+length))) {
			//We used up the old overflow so we need a bigger overflow.

			//Make a new buffer
			int newsize = overflowMaxSize + overflow-source;
			ELEM *newOverflowBuffer;
			newOverflowBuffer = (ELEM*)malloc(2*newsize*sizeof(ELEM));

			// Swap in the new overflow buffer, keeping a pointer to the old one for now.
			ELEM *oldOverflowBuffer = overflowBuffer;
			overflowBuffer=newOverflowBuffer;

			//build overflow buckets
			convertCountsToContiguousBuckets(overflowCounts, overflowBuckets, overflowBuffer);

			//pass left-to right from the overflow buffer, then from the front of the source also used as overflow buffer.
			passOverInput(oldOverflowBuffer+overflowMaxSize, oldOverflowBuffer+(2*overflowMaxSize), currentByte, dealToOverflow, noStats);
			passOverInput(source, overflow, currentByte, dealToOverflow, noStats);

			//Assign the new buffer size.
			overflowMaxSize = newsize;
			//Kill old buffer
			// "If ptr is a null pointer, the function does nothing."
			// --http://www.cplusplus.com/reference/cstdlib/free/
			free(oldOverflowBuffer);
		} else {
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
		for(int i = 0; i < bytecount; i++) {
			if((livebits >> (i*8)) & 255) {
				if(i < threshold) {
					bottomBytes[bottomBytesSize] = i;
					bottomBytesSize++;
				} else if(i >= threshold) {
					topBytes[bottomBytesSize] = i;
					topBytesSize++;
				}
			}
		}
	};

	bool hasByteLists = false;

 	auto fr = [&](const auto &bytes, const auto &numBytes) {

		const auto thresholdByte = sizeof(INT)-neededBytes;

		#ifdef DEBUG
                std::cout << "HERE! Running FR with neededBytes: " << neededBytes  << " and threshold byte " << thresholdByte << " and numBytes " << numBytes << std::endl;
		std::cout << "But numBytes might not be reliable" << std::endl;
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
					passOverInput(source, source+length, (thresholdByte), dealWithOverflow, noStats);
			} else {
					passOverInput(source, source+length, bytes[0], dealWithOverflow, gatherLiveBits);
					buildByteLists(thresholdByte);
					hasByteLists=true;
			}

                #ifdef DEBUG
		if(hasByteLists){
                	std::cout << "HERE! After the first pass we may have built the bytelists and now numBytes is " << numBytes << std::endl;
		}
                #endif

		//We have to account for the newly generated numBytes being either 1 or 2, which means we found enough dead bits to eliminate
		//passes that would have got us caught up above in the numBytes == 1 or numBytes == 2 section. Sure it's annoying, but actually
		//doing the deadbit passes would be more annoying, so rejoice!


                #ifdef DEBUG
                std::cout << "HERE! Processing Overflow" << std::endl;
                #endif

			processOverflow(bytes[0]);

                #ifdef DEBUG
                std::cout << "HERE! Swapping sources" << std::endl;
                #endif

			swap();

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
                std::cout << "HERE2!" << std::endl;
                #endif


			//Do last two passes
			countedByte=bytes[numBytes-1];
			doEstimates();
			passOverInputs(source, overflowBuffer, bytes[numBytes-2], dealWithOverflow, countForByte);
			processOverflow(bytes[numBytes-2]);
			std::swap(source, destination);

			convertCountsToContiguousBuckets(bucketCounts, destinationBuckets, destination);

			//Deal last pass in top bytes
			passOverInputs(source, overflowBuffer, countedByte, dealExact, noStats);
			swap();
		}
	};


	//
	// Process  Top Bytes
	//

	fr(topBytes, topBytesSize);

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

				fr(bottomBytes, bottomBytesSize);
				if((bottomBytesSize+topBytesSize)%2) copyBuffer(startDefiniteLongRun, endDefiniteLongRun);

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
						if(topBytesSize%2) copyBuffer(startSmallRuns, endSmallRuns);
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
                if((bottomBytesSize+topBytesSize)%2) copyBuffer(startDefiniteLongRun, endDefiniteLongRun);

                source = oldSource;
                destination = oldDestination;
                length = oldLength;
	} else {
		endSmallRuns = source+length;
		insertionSort<INT, ELEM>(startSmallRuns, endSmallRuns-startSmallRuns);
		if(topBytesSize%2) copyBuffer(startSmallRuns, endSmallRuns);
	}

	if((bottomBytesSize+topBytesSize)%2) swap(); //So we don't break the delete

	//We're done! Delete all the stuff we made
	delete [] destination;
	delete [] sourceBuckets;
	delete [] destinationBuckets;
	delete [] overflowBuckets;
	free(overflowBuffer);
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

        return 0;
}


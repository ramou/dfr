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

const unsigned BYTECOUNTER_SIZE = 257;
typedef unsigned BucketCount[BYTECOUNTER_SIZE];
typedef unsigned OverflowBucketCount[BYTECOUNTER_SIZE+1];

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
void fr(ELEM *source, auto length) {

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
	BucketCount bucketCounts;
        std::memset(bucketCounts, 0, sizeof(bucketCounts));

	//We actually count for overflow, so we must use this. We'll just convert to 
	//overflowBuckets when done because that's easier to process.
	OverflowBucketCount overflowCounts;
        std::memset(overflowCounts, 0, sizeof(overflowCounts));

	//We don't initialize this till we need it. We should check overflowMaxSize 
	//and expend it as needed, re-initializing the buffer. There may be benefit
	//to being fancy about how we choose to grow this buffer, bur for now I don't 
	//care. We may also do something like initializing it for 2% of data, which
	//should be a realistic expectation given the newer sampling approach.
	int overflowMaxSize = 0;
	ELEM *overflowBuffer;


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

	std::cout << "When sampling we found the following live bits: " << std::endl << std::bitset<sizeof(INT)*8>(livebits) << std::endl;

	const auto allBits = pow(2, sizeof(INT));
	for(auto i = 1; i < allBits; i++) {
		auto val = ((livebits >> (int)(allBits - i)) & 1);
		foundLiveBits += val;
		if(foundLiveBits==neededBits) {
			std::cout << "We processed " << i << " bits before finding enough so our top bytes are " << neededBytes << std::endl;
			neededBytes = ceil((i+1)/8.0);
			break;
		}
	}
        //std::cout << "Given that sampled livebits are " << std::endl << std::bitset<sizeof(INT)*3>(livebits) << std::endl;
        //std::cout << "we will need " << neededBytes << " passes." << std::endl;


	/**
		4) GETTING A SAMPLE DISTRIBUTION
	*/

	auto incrementDestination = [] (auto &d, auto amount) {
		d += amount;
	};

	/*
		We divide the length roughly (using bitshift) by 256 and
		allocate the space in destination evenly.

		The extra space after the rough divide is distributed by
		adding one spot to each of the first few buckets (the
		remainder) and then just allocating the rougher estimate
		for the latter buckets, thus the full length is allocated.
	*/
	auto estimateUniform = [&sourceBuckets] (auto destination, auto length, auto inc) {
		auto bucketGuess = length>>8;
                auto remainder = length-(bucketGuess<<8);

                auto d = destination;
                for(auto i = 0; i < remainder; i++) {
                        sourceBuckets[i << 1]       = d; //start
                        inc(d, bucketGuess+1);     //these buckets are one bigger than the latter buckets
                        sourceBuckets[(i << 1) + 1] = d; //end, which will be start of next bucket
                }
                for(auto i = remainder; i < 256; i++) {
                        sourceBuckets[i << 1]       = d; //start
			inc(d, bucketGuess);
                        sourceBuckets[(i << 1) + 1] = d; //end, which will be start of next bucket
                }
	};

	/*
		This will eventually need to capture more data
	*/
	auto doEstimates = [&estimateUniform, &destination, &length, &incrementDestination]() {
		if(length < DISTRIBUTION_SENSITIVE_THRESHOLD) {
			estimateUniform(destination, length, incrementDestination);
		} else {
		//THIS IS NOT DONE YET! 4b!

			estimateUniform(destination, length, incrementDestination);

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
	};

	doEstimates();
	auto currentByte = 8-neededBytes;
	auto o = source;

	std::cout << "We need " << neededBytes << " byte(s)." << std::endl;

	auto passOverInput = [&source, &length, &o, &livebits, &bitmask, &currentByte, &sourceBuckets, &overflowCounts, &destination](auto deal, auto &count) {
		std::for_each(source, source+length, [&livebits, &bitmask, &currentByte, &sourceBuckets, &overflowCounts, &o, &destination, &deal, &count](auto &e){
	                auto t = ((reinterpret_cast<INT>(e))>>currentByte) & 255;
                	auto d = sourceBuckets[t<<1];
        	        livebits |= bitmask ^ reinterpret_cast<INT>(e);
			deal(d,t,e,o);
			count(d,t,e);
	        });
	};

	auto dealToOverflow = [&sourceBuckets, &overflowCounts](auto &d, auto &t, auto &e, auto &o) {
		if(d < sourceBuckets[(t<<1)+1]) {
			*d=e;
			sourceBuckets[t<<1]++;
		} else {
			std::cout << "Placed value with byte " << t << " in overflow bucket" << std::endl;
                        overflowCounts[t]++;
                        *o=e;
                        o++;
		}
	};

	auto noCount = [](auto &d, auto &t, auto &e){};

	passOverInput(dealToOverflow,noCount);


	/**
		5) Process first pass of top bytes
	*/

	std::cout << "We overflowed " << (o-source) << " times." << std::endl;
	std::cout << "We now known the following are live bits: " << std::endl << std::bitset<sizeof(INT)*8>(livebits) << std::endl;

	//deal overflow into overflow buffer
	overflowMaxSize=(o-source);
	overflowBuffer=(ELEM*)malloc(2*overflowMaxSize*sizeof(ELEM));

	std::cout << "We made a buffer of length " << (2*overflowMaxSize) <<  std::endl;

		//convert overflow counts into target buckets
	overflowBuckets[0]=overflowBuffer;
	auto currentBucket = &overflowBuckets[1];
	auto previousBucket = &overflowBuckets[0];

	std::for_each(overflowCounts, overflowCounts+256, [&currentBucket, &previousBucket](auto &count) {
		currentBucket[0]=previousBucket[0]+count;
		currentBucket++;
		previousBucket++;
	});

		//Stuart:: so far this looks right

	std::cout << "We succeeded in building overflow buckets."  << std::endl;


		//deal into overflow buffer
	std::for_each(source, o, [&currentByte, &overflowBuckets](auto &e) {
		auto t = ((reinterpret_cast<INT>(e))>>currentByte) & 255;
		auto d = overflowBuckets[t];
		*d=e;
		overflowBuckets[t]++;
	});

	std::cout << "We succeeded in dealing into the overflow buckets."  << std::endl;

	std::for_each(overflowBuffer, overflowBuffer+(overflowMaxSize-1), [&currentByte](auto &e) {
		auto t = ((reinterpret_cast<INT>(e))>>currentByte) & 255;
		std::cout << "overflow value: " << t << std::endl;
	});

	std::cout << "Now we swap buffers."  << std::endl;

	//swap source and destination
	std::swap(source, destination);
	std::cout << "Now we swap buffer counts."  << std::endl;
	std::swap(sourceBuckets, destinationBuckets);


	//Do the rest of current bytes (we started based on neededbytes)
	currentByte++;

	//At this point, source has our data in disjoint buckets and destination is ready.
	//buckets contains the start and end positions of each disjoint bucket in source.
	//overflowBuffer is the contiguous secondary source, with end buckets defined by 
	//overflowBuckets, and is twice overflowMaxSize so we can deal into it to before 
	//we make use of the input as a buffer source to prevent writing over a live bucket.

	//redo estimates

	std::cout << "Now we redo our estimates."  << std::endl;

	doEstimates();


	//Between passes, use std::swap as necessarty for source and destination.

	//Now we should define our lambdas.


	//ERASE THIS! This just fixes our current state into source since we've 
	//done an odd number of passes and we're still writing code!!!
	std::swap(source, destination);

	//We're done! Delete all the stuff we made
	delete [] destination;
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

	fr<targetType, targetType>(values, targetLength);

	for(int i = 1; i < targetLength; i++) {
		if(values[i] < values[i-1]) {
			std::cout << " value " << i << " is out of place." << std::endl;
			break;
		}
	}

        return 0;
}


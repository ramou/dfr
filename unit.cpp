#define CATCH_CONFIG_MAIN

#include <bits/stdc++.h> 
#include "catch.hpp"
#include "fr.hpp"

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


template <typename T>
struct elem {
        T key;
        T *val;
};


typedef uint64_t targetType;

int foo(int argc, char *argv[]) {

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

bool theSame(const auto &start1, const auto &start2, const auto &length) {
	for(auto i = 0; i <  length; i++) {
		if(start1[i] != start2[i]) {
			return false;
		}
	}
	return true;
}

TEST_CASE( "[up8to3] guess 8 live bytes, actually 3, including lowest", "[up8to3]" ) {

	uint64_t safe[15] = {65536, 1, 1, 4, 5, 6, 7, 8, 512, 512, 11, 12, 13, 4, 4};
	uint64_t test[15] = {65536, 1, 1, 4, 5, 6, 7, 8, 512, 512, 11, 12, 13, 4, 4};
	int n = sizeof(safe)/sizeof(safe[0]); 
	std::sort(safe, safe+n);
        #ifdef DEBUG
	std::cout << "We want to get: ";
	for(int i = 0; i < 15; i++) std::cout << safe[i] << " ";
	std::cout << std::endl;
        #endif

	dfr<uint64_t, uint64_t>(test, 15);

        #ifdef DEBUG
        std::cout << "We got: ";
        for(int i = 0; i < 15; i++) std::cout << test[i] << " ";
        std::cout << std::endl;
        #endif



	REQUIRE( theSame(safe, test, 15));
}

TEST_CASE( "[up8to2] guess 8 live bytes, actually 2, including lowest", "[up8to2]" ) {

        uint64_t safe[15] = {513, 1, 1, 4, 5, 6, 7, 8, 512, 512, 11, 12, 13, 4, 4};
        uint64_t test[15] = {513, 1, 1, 4, 5, 6, 7, 8, 512, 512, 11, 12, 13, 4, 4};
        int n = sizeof(safe)/sizeof(safe[0]);

        #ifdef DEBUG
        std::cout << "We started with: " << std::endl;
        for(int i = 0; i < 15; i++) std::cout << safe[i] << " ";
        std::cout << std::endl;
        #endif

        std::sort(safe, safe+n);
        #ifdef DEBUG
        std::cout << "We want to get: " << std::endl;
        for(int i = 0; i < 15; i++) std::cout << safe[i] << " ";
        std::cout << std::endl;
        #endif

        dfr<uint64_t, uint64_t>(test, 15);
        #ifdef DEBUG
        std::cout << "We got: " << std::endl;
        for(int i = 0; i < 15; i++) std::cout << test[i] << " ";
        std::cout << std::endl;
        #endif


        REQUIRE( theSame(safe, test, 15));
}

TEST_CASE( "[up8to1]guess 8 live bytes, actually 1, including lowest", "[up8to1]" ) {

        uint64_t safe[15] = {254, 1, 1, 4, 5, 6, 7, 8, 253, 52, 11, 12, 13, 4, 4};
        uint64_t test[15] = {254, 1, 1, 4, 5, 6, 7, 8, 253, 52, 11, 12, 13, 4, 4};
        int n = sizeof(safe)/sizeof(safe[0]);
        std::sort(safe, safe+n);
        dfr<uint64_t, uint64_t>(test, 15);

        REQUIRE( theSame(safe, test, 15));
}



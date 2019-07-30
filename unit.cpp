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


TEST_CASE( "[up8to1_startat1]guess 8 live bytes, actually 1, not including lowest", "[up8to1_startat1]" ) {

        uint64_t safe[15] = {1<<8, 1<<8, 1<<8, 1<<8, 5<<8, 6<<8, 7<<8, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        uint64_t test[15] = {1<<8, 1<<8, 1<<8, 1<<8, 5<<8, 6<<8, 7<<8, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        int n = sizeof(safe)/sizeof(safe[0]);
        std::sort(safe, safe+n);
        dfr<uint64_t, uint64_t>(test, 15);

        REQUIRE( theSame(safe, test, 15));
}

TEST_CASE( "[up8to2_startat1]guess 8 live bytes, actually 2, not including lowest", "[up8to2_startat1]" ) {

        uint64_t safe[15] = {1<<8, 1<<8, 1<<8, 1<<8, 5<<8, 257<<8, 7<<8, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        uint64_t test[15] = {1<<8, 1<<8, 1<<8, 1<<8, 5<<8, 257<<8, 7<<8, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        int n = sizeof(safe)/sizeof(safe[0]);
        std::sort(safe, safe+n);
        dfr<uint64_t, uint64_t>(test, 15);

        REQUIRE( theSame(safe, test, 15));
}

TEST_CASE( "[up8to3_startat1]guess 8 live bytes, actually 3, not including lowest", "[up8to3_startat1]" ) {

        uint64_t safe[15] = {1<<8, 1<<8, 1<<8, 1<<8, 5<<8, 2<<16, 7<<24, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        uint64_t test[15] = {1<<8, 1<<8, 1<<8, 1<<8, 5<<8, 2<<16, 7<<24, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        int n = sizeof(safe)/sizeof(safe[0]);
        std::sort(safe, safe+n);
        dfr<uint64_t, uint64_t>(test, 15);

        REQUIRE( theSame(safe, test, 15));
}

TEST_CASE( "[up7to1_startat1]guess 7 live bytes, actually 1, including lowest", "[up7to1_startat1]" ) {

        uint64_t safe[15] = {254<<8, 1<<8, 1<<8, 4<<8, 5<<8, 6<<8, 7<<8, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        uint64_t test[15] = {254<<8, 1<<8, 1<<8, 4<<8, 5<<8, 6<<8, 7<<8, 8<<8, 253<<8, 52<<8, 11<<8, 12<<8, 13<<8, 4<<8, 4<<8};
        int n = sizeof(safe)/sizeof(safe[0]);
        std::sort(safe, safe+n);
        dfr<uint64_t, uint64_t>(test, 15);

        REQUIRE( theSame(safe, test, 15));
}


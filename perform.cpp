#include <bits/stdc++.h> 
#include "fr.hpp"

bool theSame(const auto &start1, const auto &start2, const auto &length);

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

int main(int argc, char *argv[]) {

        auto targetLength = 512;

        if(argc >= 2) {
                targetLength = atoi(argv[1]);
        }

        targetType* values = new targetType[targetLength];
	targetType* safe = new targetType[targetLength];

        if(argc > 2) {
		#ifdef DEBUG
                std::cout << "Making fixed input of size " << targetLength << " with values: " << std::flush;
		#endif
                for(int i = 0; i < targetLength; i++) {
                  auto val = atoi(argv[i+2]);
		#ifdef DEBUG
                  std::cout << val << " " << std::flush;
		#endif
                  values[i] = val;
                }
		#ifdef DEBUG
                std::cout << std::endl;
		#endif
        } else {
		#ifdef DEBUG
                std::cout << "Making random input of size " << targetLength << std::endl;
		#endif
                make_random(values, targetLength);
        }

	for(int i = 0; i < targetLength; i++) {
		safe[i] = values[i];
	}

        dfr<targetType, targetType>(values, targetLength);
	std::sort(safe, safe+targetLength);

	if(!theSame(safe, values, targetLength)) {
		 std::cout << "Values don't match expected sorted values" << std::endl;
	}

        delete [] values;
	delete [] safe;

        return 0;
}

bool theSame(const auto &start1, const auto &start2, const auto &length) {
	for(auto i = 0; i <  length; i++) {
		if(start1[i] != start2[i]) {
			std::cout << "Value " << i << " is out of place." << std::endl;
			return false;
		}
	}
	return true;
}


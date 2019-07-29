
default: buildtest ;

buildtest:
	c++ -o3 -std=c++14 unit.cpp -o test

tests: buildtest ;
	./test

debug:
	c++ -g -std=c++14 -DDEBUG unit.cpp -o test
	./test $(filter-out $@,$(MAKECMDGOALS))

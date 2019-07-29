MISSING_VALGRIND="No valgrind in $(PATH), consider doing apt-get install valgrind (or yum install valgrind)"

default: buildtest ;

valgrind:
	@if [ -n "$(shell which valgrind)" ]; then \
	echo "c++ -o3 -std=c++14 -DNOTRANDOM perform.cpp -o perform" ; \
	c++ -o3 -std=c++14 -DNOTRANDOM perform.cpp -o perform ; \
	echo "valgrind --error-exitcode=2 -q --leak-check=yes ./perform 15 1 1 65536 4 5 6 7 8 512 512 11 12 13 4 4" ; \
	if ! valgrind --error-exitcode=2 -q --leak-check=yes ./perform 15 1 1 65536 4 5 6 7 8 512 512 11 12 13 4 4 ; then \
	exit 2; \
	fi; \
	else \
	echo $(MISSING_VALGRIND); \
	exit 1; \
	fi

buildtest:
	c++ -o3 -std=c++14 -DNOTRANDOM unit.cpp -o test

tests: buildtest ;
	./test

debug:
	c++ -g -std=c++14 -DDEBUG -DNOTRANDOM unit.cpp -o test
	./test $(filter-out $@,$(MAKECMDGOALS))

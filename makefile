MISSING_VALGRIND="No valgrind in $(PATH), consider doing apt-get install valgrind (or yum install valgrind)"

default: ;
	c++ -o3 -std=c++14 perform.cpp -o perform

timing: ;
	c++ -o3 -DTIMINGS -std=c++14 perform.cpp -o perform
	./perform $(filter-out $@,$(MAKECMDGOALS))

valgrind: ;
	@if [ -n "$(shell which valgrind)" ]; then \
	echo "c++ -g -std=c++14 -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test" ; \
	c++ -g -std=c++14 -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test ; \
	echo "valgrind --error-exitcode=2 -q --leak-check=yes ./test $(filter-out $@,$(MAKECMDGOALS))" ; \
	if ! valgrind --error-exitcode=2 -q --leak-check=yes ./test $(filter-out $@,$(MAKECMDGOALS)) ; then \
	exit 2; \
	fi; \
	else \
	echo $(MISSING_VALGRIND); \
	exit 1; \
	fi

buildtest: ;
	c++ -g -std=c++14 -DDEBUG -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test
	c++ -g -std=c++14 -DDEBUG -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 perform.cpp -o perform

tests: buildtest ;
	./test

debug: ;
	c++ -g -std=c++14 -DDEBUG -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test
	./test -s $(filter-out $@,$(MAKECMDGOALS))

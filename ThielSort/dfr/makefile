MISSING_VALGRIND="No valgrind in $(PATH), consider doing apt-get install valgrind (or yum install valgrind)"

# Host / CPU details
##### Must export HOSTNAME for make to see it!!
##### NOTE: on titan, this is good for some versions, not others:-((

##### NOTE: PREFER_NewBuffer:=true fails when sharing buffers across
##### recursions because we are allocating 2 temp buffers in
##### preference to using the original space as either the source or
##### the dest each pass. Since each DataPool doesn't have its own
##### allocator any more, we try to reuse the same buffer for both.
ifeq ($(HOSTNAME),HEX-GPGPU)
PREFER_NewBuffer:=true
ARCH_EXTRAS:=
else
PREFER_NewBuffer:=false 
endif
CPUExtras=-D PREFER_NEW_BUFFER=$(PREFER_NewBuffer) 

### The TimerMode variable will override default timing options for
### all compile requests. TimerMode can be 0 (no timings), 1 (total
### time only) or 2 (all relevant timings are done).
ifdef TimerMode
# make all defaults TimerMode
DfltTmrMode0 := $(TimerMode)
DfltTmrMode1 := $(TimerMode)
DfltTmrMode2 := $(TimerMode)
else
# Default timing mode for various goals
DfltTmrMode0 := 0
DfltTmrMode1 := 1
DfltTmrMode2 := 2
endif

TimeStamp := $(shell date)

Open=@echo "{  $(TimeStamp)" >> $(LogFile)
Doll=@ls -l $(ExecFile) >> $(LogFile)
Close=@echo "}  $(TimeStamp)" >> $(LogFile)

# Compiler defaults.
### Simd parameters (AVX, AVX2).
ifndef SimdLevel
SimdLevel:=AVX
endif
ifeq ($(SimdLevel), AVX)
SIMD_PARAMS:= -D SIMD_STREAM_STORE_SIZE=5 -D SIMD_MOVE_SIZE=5 \
        -D SIMD_SHIFT_SIZE=4 -D SIMD_MASK_SIZE=5 \
        -D SIMD_STREAM_LOAD_SIZE=4 -mavx
endif
ifeq ($(SimdLevel), AVX2)
SIMD_PARAMS:= -D SIMD_STREAM_STORE_SIZE=5 -D SIMD_MOVE_SIZE=5 \
        -D SIMD_SHIFT_SIZE=5 -D SIMD_MASK_SIZE=5 \
        -D SIMD_STREAM_LOAD_SIZE=5 -mavx2
endif

### Add AVX2 here to try gather instr for loading interleaved.
### Above needs mod to use AVX2 in other contexts (unless gcc uses it).
GccExtras=-Wno-attributes -I ../SsV3 -I ../SimdSupport \
	-flax-vector-conversions -m64 $(SIMD_PARAMS)

##	-I ~/Devel/src/InfraStructure/isominc/include \
# -std=c++14      : specify C++ standard 
# -Wno-attributes : Block 'maybe not inlineable' warning.

# Configuration defaults:
Extras=
ShotgunExtras=-DIL_VALUE=1\
	-DCONFIG_PATTERN=SD_DIRECT_DIVERTING_WITH_PF_1BY1_PATTERN

SimdGunExtras=-DIL_VALUE=1\
	-DCONFIG_PATTERN=SD_DIRECT_DIVERTING_WITH_PF_NTSIMD_PATTERN\
	 -DLDL_N_TO_PRIME=2

June12Extras:=-DCONFIG_PATTERN=SD_DIRECT_DIVERTING_WITH_PF_NTSIMD_PATTERN \
	-DIL_VALUE=1 -DFLUSH_MOVER_VARIANT=1 -DLDL_N_TO_PRIME=4 \
	-DBKT_N_TO_PRIME=10 -DDC_POST_INC=0 -DDEALER_PF_LEAD=2 \
        -DBATCH_SIZE=24000 -DN_TPATHS=4 -DDIP_TIMING=0 -DDEAL_PF_FAR=6

July18SimdExtras:=-DCONFIG_PATTERN=SD_DIR_DIV_HYBRID_WITH_PF_NTSIMD_PATTERN \
	-DIL_VALUE=1 -DFLUSH_MOVER_VARIANT=1 -DLDL_N_TO_PRIME=4 \
	-DBKT_N_TO_PRIME=10 -DDC_POST_INC=0 -DDEALER_PF_LEAD=2 \
	-DSimdLevel=AVX2 -DBATCH_SIZE=24000 

AllExtras = $(CPUExtras) $(GccExtras) $($(CannedExtras)) $(Extras)

##### Specify the log file for most output (not cppout.txt).
LogFile:=Times.txt

default: ;
	c++ -g $(AllExtras) -DOCD_TESTING Main.c  2>&1 | tee  cppout.txt

preprocess: ;
	c++ -E $(AllExtras) -DOCD_TESTING Main.c  > preprocess.txt

passtiming: ;
	c++ -O3 -DTIMINGS=$(DfltTmrMode2) $(AllExtras) Main.c -o dfrOpt
	@echo Compiling : Extras = $(Extras) | tee -a $(LogFile)
	./dfrOpt $(filter-out $@,$(MAKECMDGOALS))

timing: ;
	c++ -O3 -DTIMINGS=$(DfltTmrMode1) $(AllExtras) Main.c -o dfrOpt
	@echo Compiling : Extras = $(Extras) | tee -a $(LogFile)
	./dfrOpt $(filter-out $@,$(MAKECMDGOALS))

shotgun: ;
	c++ -O3 -DTIMINGS=$(DfltTmrMode1) $(ShotgunExtras) $(AllExtras) \
	Main.c -o shotgun
	@echo Compiling : Extras = $(Extras) | tee -a $(LogFile)
	./shotgun $(filter-out $@,$(MAKECMDGOALS))

simdgun: ;
	c++ -O3 -DTIMINGS=$(DfltTmrMode1) $(SimdGunExtras) $(AllExtras) \
	Main.c -o simdgun
	@echo Compiling : Extras = $(Extras) | tee -a $(LogFile)
	./simdgun $(filter-out $@,$(MAKECMDGOALS))

disass: ;
	c++ -O3 -DTIMINGS=$(DfltTmrMode0) $(AllExtras) -DGDBHOOK Main.c \
	-o dfrOpt
	./dfrOpt $(filter-out $@,$(MAKECMDGOALS))
	gdb ./dfrOpt.exe

valgrind: ;
	@if [ -n "$(shell which valgrind)" ]; then \
	echo "c++ -g $(AllExtras) -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test" ; \
	c++ -g $(AllExtras) -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test ; \
	echo "valgrind --error-exitcode=2 -q --leak-check=yes ./test $(filter-out $@,$(MAKECMDGOALS))" ; \
	if ! valgrind --error-exitcode=2 -q --leak-check=yes ./test $(filter-out $@,$(MAKECMDGOALS)) ; then \
	exit 2; \
	fi; \
	else \
	echo $(MISSING_VALGRIND); \
	exit 1; \
	fi

### Set executable to use for multitime and E2ToE9 goals
ExecFile:=./dfrOpt

## By default, verify sort all the time.
RunOptions=-v
E2E9:=100 1000 10000 100000 1000000 10000000 100000000 1000000000
## E2 and E3 too small for me to care!
## E9 is long on desktop. Drop it for mass testing.
E4E8:=10000 100000 1000000 10000000 100000000
## But E9 is the traditional goal, so this version includes it
E4E9:=10000 100000 1000000 10000000 100000000 1000000000

MultiSet:=1000 5000 10000 25000 50000 100000 \
	250000 500000 1000000 2500000 5000000 10000000 25000000 50000000 \
	100000000 250000000 500000000 750000000 1000000000

NSet:
	$(Open)
	$(Doll)
	$(ExecFile) | tee -a $(LogFile)
	for n in $(CurSet) $($(NamedSet)) ; do \
	$(ExecFile) $(RunOptions) -n $$n | tee -a $(LogFile) ; done
	$(Close)

debug: ;
	c++ -g $(AllExtras) -DDEBUG -DNOTRANDOM -DTEST_THRESHOLD=4 -DTEST_LADLE_THRESHOLD=10 unit.cpp -o test
	./test -s $(filter-out $@,$(MAKECMDGOALS))

N:=50000000
logrun:
	@date >> $(LogFile)
	$(ExecFile) $(N) | tee -a $(LogFile)

log3:
	@echo "***************************************log3" >> $(LogFile)
	$(Open)
	$(Doll)
	##Large time for first run with large N is lack of -s!
	$(ExecFile) $(N) | tee -a $(LogFile)
	$(ExecFile) -s $(N) | tee -a $(LogFile)
	$(ExecFile) -s $(N) | tee -a $(LogFile)
	@echo "************************************end log3" >> $(LogFile)
	$(Close)
help:
	@echo "'make timing Extras=\"...\"' compiles with pass timing, -O3 $(GccExtras)."
	@echo "'Extras' selects code versions and specifies other options. "
	@echo "Default for Extras is $(Extras) ."
	@echo "'make passtiming' is same as above. Distinction between total "
	@echo "time and per pass times no longer exists"
	@echo "Current code options (blank delimited in Extras):"
	@echo " "
	@echo "'make disass' compiles like timing above, but enables calls to GdbHook()"
	@echo "so I can locate where I am and disass a function. It invokes gdb after the compile"
	@echo "and test run. It doesn't log anything in $(LogFile)."
	@echo " "
	@echo "'make' compiles without optimization or timing code and with -g"
	@echo "It echos Extras and launches gdb."
	@echo " "
	@echo "'make logrun N=nnn' runs a job to sort nnn values and appends a timestamp and"
	@echo "any output (e.g., timings) to local file $(LogFile)"
	@echo "The default for $(LogFile) is Times.txt. "
	@echo "Use LogFile=xxx to log to file xxx"
	@echo "Use N=nnn to specify n for sort. Default = $(N)."
	@echo " "
	@echo "'make multitime' to run timing tests for several Ns"
	@echo "with $(EXECFILE). Default for EXECFILE is ./dfrOpt
 

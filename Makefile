# make          builds constp and constp2
# make check    builds, then runs the test suite up to n = 11
# make clean    removes the binaries
#
# No dependencies beyond a C++17 compiler. Two things are probed rather than assumed:
#
#   -march=x86-64-v2  is the right portability floor on x86 (SSE4.2 and POPCNT, present
#                     on anything since about 2009) but is not a valid flag on ARM, so
#                     it is only used if the compiler accepts it.
#   -fopenmp          Apple's clang does not ship it. Without it the build is
#                     single-threaded, which is correct, just slower. For threads on
#                     macOS: brew install gcc, then  make CXX=g++-14

CXX ?= g++

ARCH := $(shell printf 'int main(){}' > /tmp/_arch.cpp; \
          $(CXX) -march=x86-64-v2 /tmp/_arch.cpp -o /tmp/_arch 2>/dev/null && echo -march=x86-64-v2)
OMP  := $(shell printf 'int main(){}' > /tmp/_omp.cpp; \
          $(CXX) -fopenmp /tmp/_omp.cpp -o /tmp/_omp 2>/dev/null && echo -fopenmp)

CXXFLAGS ?= -O3 $(ARCH) -std=c++17

all: constp constp2
	@[ -n "$(OMP)" ] || echo "note: no OpenMP in $(CXX); built single-threaded (correct, just slower)"

constp:  constp.cpp  ; $(CXX) $(CXXFLAGS) -o $@ $<
constp2: constp2.cpp ; $(CXX) $(CXXFLAGS) $(OMP) -o $@ $<

check: all
	python3 verify.py 11
test: check

clean: ; rm -f constp constp2

.PHONY: all check test clean

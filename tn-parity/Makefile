# make          builds wronskian, cert and dense
# make check    about a minute: dense reference, certificate and DP all agree
#
# See the top-level Makefile for why the arch and OpenMP flags are probed rather than
# assumed. For threads on macOS: brew install gcc, then  make CXX=g++-14

CXX ?= g++

ARCH := $(shell printf 'int main(){}' > /tmp/_arch2.cpp; \
          $(CXX) -march=native /tmp/_arch2.cpp -o /tmp/_arch2 2>/dev/null && echo -march=native)
OMP  := $(shell printf 'int main(){}' > /tmp/_omp2.cpp; \
          $(CXX) -fopenmp /tmp/_omp2.cpp -o /tmp/_omp2 2>/dev/null && echo -fopenmp)

CXXFLAGS ?= -O3 $(ARCH) -std=c++17

all: wronskian cert dense
	@[ -n "$(OMP)" ] || echo "note: no OpenMP in $(CXX); built single-threaded (correct, just slower)"

wronskian: wronskian.cpp ; $(CXX) $(CXXFLAGS) $(OMP) -o $@ $<
cert:      cert.cpp      ; $(CXX) $(CXXFLAGS) $(OMP) -o $@ $<
dense:     dense.cpp     ; $(CXX) $(CXXFLAGS) -o $@ $<

check: all
	@for k in 3 5 7 9 11 13; do ./dense $$k; done
	@for k in 3 5 7 9 11 13 15; do ./cert $$k 5 | tail -1; done

clean: ; rm -f wronskian cert dense

.PHONY: all check clean

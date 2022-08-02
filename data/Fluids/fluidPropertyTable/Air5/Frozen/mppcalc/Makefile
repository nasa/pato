
# Change this path to your Mutation++ install path
MPP_INSTALL = $(MPP_DIRECTORY)/install

CXX       = g++ --std=c++11
CXX_FLAGS = -I $(MPP_INSTALL)/include/mutation++ -I $(MPP_INSTALL)/include -O3 -g
# use .dylib for Mac, .so for Unix
CXX_LIBS  = $(MPP_INSTALL)/lib/libmutation++.so

all : mppcalc

mppcalc : mppcalc.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS) $<

clean :
	$(RM) *.o mppcalc


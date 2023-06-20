#COMPILER MODE C++11
CXX=g++ -std=c++11

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
#HTSLIB_INC=$(HOME)/Tools/htslib-1.11
#HTSLIB_LIB=$(HOME)/Tools/htslib-1.11/libhts.a
HTSLIB_INC=/opt/miniconda3/pkgs/htslib-1.12-hc38c3fb_1/include
HTSLIB_LIB=/opt/miniconda3/pkgs/htslib-1.12-hc38c3fb_1/lib/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
#BOOST_INC=/usr/include
BOOST_INC=/opt/miniconda3/pkgs/boost-cpp-1.74.0-hff03dee_4/include
#BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
BOOST_LIB_IO=/opt/miniconda3/pkgs/boost-cpp-1.74.0-hff03dee_4/lib/libboost_iostreams.dylib
#BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
BOOST_LIB_PO=/opt/miniconda3/pkgs/boost-cpp-1.74.0-hff03dee_4/lib/libboost_program_options.dylib



#COMPILER & LINKER FLAGS

#Best performance is achieved with this. Use it if running on the same plateform you're compiling.
#CXXFLAG=-O3 -march=native

#Good performance and portable on most CPUs
CXXFLAG=-O3 -mavx2 -mfma 

#Portable version without avx2 (much slower)
#CXXFLAG=-O3

#LDFLAG=-O3
LDFLAG="-I/opt/local/include -L/opt/local/lib"


#DYNAMIC LIBRARIES
#DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl -lssl -lcrypto
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl

#SHAPEIT SOURCES & BINARY
BFILE=bin/shapeit4.2
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#COMPILATION RULES
all: $(BFILE)

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)

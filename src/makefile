.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lbz2 -llzma -lcurl -lstdc++
SOURCE		= main.cpp local_mode.cpp intron.cpp variant.cpp tools.cpp hollywood.cpp svm_bp_finder.cpp RNAcofold.cpp misc.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o SpliceApp $(LIB)
			
clean:
		rm -f *.o *~

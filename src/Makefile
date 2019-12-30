CXX             = g++ 
#CXXFLAGS        = -g -DDEBUG
LDFLAGS        = -L. -lm -lRNA
all: RNAmountAlignScan RNAmountAlign

RNAmountAlignScan: RNAmountAlignScan.o pair_align.o arguments.o aux.o mountain_height.o evd.o stats.o karlin_altschul.o
	${CXX} ${CXXFLAGS} -o RNAmountAlignScan RNAmountAlignScan.o pair_align.o arguments.o aux.o mountain_height.o evd.o stats.o karlin_altschul.o ${LDFLAGS}

RNAmountAlign: RNAmountAlign.o pair_align.o arguments.o aux.o mountain_height.o stats.o karlin_altschul.o evd.o multi_align.o profile.o profile_aligner.o 
	${CXX} ${CXXFLAGS} -o RNAmountAlign RNAmountAlign.o pair_align.o arguments.o aux.o mountain_height.o stats.o karlin_altschul.o evd.o multi_align.o profile.o profile_aligner.o ${LDFLAGS}

RNAmountAlignScan.o:RNAmountAlignScan.cpp RNAmountAlignScan.h arguments.h pair_align.h mountain_height.h aux.h evd.h stats.h karlin_altschul.h
	${CXX} ${CXXFLAGS} -c RNAmountAlignScan.cpp

RNAmountAlign.o:RNAmountAlign.cpp arguments.h pair_align.h mountain_height.h aux.h
	${CXX} ${CXXFLAGS} -c RNAmountAlign.cpp

pair_align.o:pair_align.cpp aux.h aux.cpp
	${CXX} ${CXXFLAGS} -c pair_align.cpp
		
arguments.o:arguments.cpp aux.h
	${CXX} ${CXXFLAGS} -c arguments.cpp 
	
aux.o:aux.cpp pair_align.cpp
	${CXX} ${CXXFLAGS} -c aux.cpp

mountain_height.o:mountain_height.cpp
	${CXX} ${CXXFLAGS} -c mountain_height.cpp

stats.o:stats.cpp stats.h
	${CXX} ${CXXFLAGS} -c stats.cpp
	
evd.o:evd.cpp evd.h
	${CXX} ${CXXFLAGS} -c evd.cpp

karlin_altschul.o: karlin_altschul.cpp karlin_altschul.h
	${CXX} ${CXXFLAGS} -c karlin_altschul.cpp

multi_align.o:multi_align.cpp multi_align.h aux.h profile_aligner.h profile_aligner.h pair_align.h
	${CXX} ${CXXFLAGS} -c multi_align.cpp

profile_aligner.o:profile_aligner.cpp profile_aligner.h aux.h aux.cpp profile_aligner.h profile_aligner.cpp
	${CXX} ${CXXFLAGS} -c profile_aligner.cpp


profile.o:profile.cpp profile.h aux.h aux.cpp
	${CXX} ${CXXFLAGS} -c profile.cpp

clean:
	rm -fr *.o RNAmountAlignScan RNAmountAlign

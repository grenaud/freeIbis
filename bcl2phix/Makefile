

CXX = g++ #-g 
CXXFLAGS+= -O3 -Wall 
LDLIBS+= -lm

all:	bcl2phix combinePhixFiles

%.o:    %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

bcl2phix:	bcl2phix.o
	${CXX}  -o $@ $^ $(LDLIBS)

combinePhixFiles:	combinePhixFiles.o
	${CXX}  -o $@ $^ $(LDLIBS)


clean:
	rm -f bcl2phix combinePhixFiles *.o 

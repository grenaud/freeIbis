all:	estimaErrorLookupMax LinRegFixed.o  findBestk

estimaErrorLookupMax:
	g++ -lm -O3 -o estimaErrorLookupMax estimaErrorLookupMax.cpp 

LinRegFixed.o:
	g++ -lm -O3 -c LinRegFixed.cpp 

findBestk:	LinRegFixed.o
	g++ -lm -o findBestk findBestk.cpp LinRegFixed.o

clean:
	rm -f estimaErrorLookupMax findBestk LinRegFixed.o


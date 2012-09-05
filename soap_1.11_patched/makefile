CC=	g++

BIN = $(DESTDIR)/usr/bin
FLAGS= -static -DMAXGAP=3 -DMAXHITS=10000 -DTHREAD -O3 #-g -fno-inline #-DTHREAD -O3 -g
REF_CHR=  -DDB_CHR
REF_CONTIG=	-DDB_CONTIG
REF_SHORT=	-DDB_SHORT
REF_HUGE=	-DDB_HUGE
OLIGOLEN= -DREAD_84 #-DREAD_84 -DREAD_72 -DREAD_60 -DREAD_48 -DREAD_36

LIBS=  
THREAD=	-lpthread

SOURCE = align dbseq main pairs param reads utilities
OBJS1= $(patsubst %,%.chr.o,$(SOURCE))
OBJS2= $(patsubst %,%.contig.o,$(SOURCE))
OBJS3= $(patsubst %,%.short.o,$(SOURCE))
OBJS4= $(patsubst %,%.huge.o,$(SOURCE))

all: soap soap.contig soap.short soap.huge soap_dealign

%.chr.o:%.cpp
	$(CC) $(FLAGS) $(LIBS) $(REF_CHR) $(OLIGOLEN) -c $< -o $@
soap: $(OBJS1)
	$(CC) $(FLAGS) $(LIBS) $(REF_CHR) $(OLIGOLEN) $^ -o $@ $(THREAD)

%.contig.o:%.cpp
	$(CC) $(FLAGS) $(LIBS) $(REF_CONTIG) $(OLIGOLEN) -c $< -o $@
soap.contig: $(OBJS2)
	$(CC) $(FLAGS) $(LIBS) $(REF_CONTIG) $(OLIGOLEN) $^ -o $@ $(THREAD)

%.short.o:%.cpp
	$(CC) $(FLAGS) $(LIBS) $(REF_SHORT) $(OLIGOLEN) -c $< -o $@
soap.short: $(OBJS3)
	$(CC) $(FLAGS) $(LIBS) $(REF_SHORT) $(OLIGOLEN) $^ -o $@ $(THREAD)
	
%.huge.o:%.cpp
	$(CC) $(FLAGS) $(LIBS) $(REF_HUGE) $(OLIGOLEN) -c $< -o $@
soap.huge: $(OBJS4)
	$(CC) $(FLAGS) $(LIBS) $(REF_HUGE) $(OLIGOLEN) $^ -o $@ $(THREAD)

soap_dealign: soap_dealign.o dealign.o param.chr.o
	$(CC) $(FLAGS) $(LIBS) $^ -o $@
				
clean:
	rm -f *.o soap soap.contig soap.short soap.huge soap_dealign
install:
	install -d $(BIN)
	install ./soap ./soap.contig ./soap.short ./soap.huge ./soap_dealign $(BIN)
#	cp soap soap.contig soap.short soap.huge soap_dealign ../../ #../bin

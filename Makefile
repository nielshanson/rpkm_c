CC = g++ -Wall -ggdb
#CC = g++ -pg
CC = g++ -Wall
CCFLAGS=


PROG = rpkm
SOURCES= utilities.c++ rpkm.c++ helper.c++ fastareader.c++ matchoutputparser.c++
OBJECTS= $(SOURCES:.c++=.o)
HEADERS= $(SOURCES:.c++=.h)

%.o: %.c++   $(SOURCES) types.h
	$(CC)  $< -c -o $@  

all: $(PROG)

clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS) $(HEADERS) types.h
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)


test1:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta --r1 data/IIYH-se.sam  -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest
test2:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta  --r2 data/IIYH-pe.sam -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest
test3:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta --r1 data/IIYH-se.sam --r2 data/IIYH-pe.sam -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest

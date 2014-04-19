CC = g++ -Wall -ggdb
#CC = g++ -pg
CC = g++ -Wall
CCFLAGS=


PROG = rpkm
OBJECTS= utilities.o rpkm.o helper.o fastareader.o
SOURCES= utilities.c++ rpkm.c++ helper.c++ fastareader.c++

%.o: %.c++   $(SOURCES)
	$(CC)  $< -c -o $@  

all: $(PROG)

clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)


test:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta --r1 data/IIYH-se.sam --r2 data/IIYH-pe.sam -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest

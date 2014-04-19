CC = g++ -Wall -ggdb
#CC = g++ -pg
CC = g++
CCFLAGS=


PROG = rpkm
OBJECTS= utilities.o rpkm.o
SOURCES= utilities.c++ rpkm.c++

%.o: %.c++   $(SOURCES)
	$(CC)  $< -c -o $@  

all: $(PROG)

clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)


test:
	$(PROG) test.fasta

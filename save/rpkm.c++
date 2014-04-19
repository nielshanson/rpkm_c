#include "utilities.h"

int main( int argc, char **argv ){
    if( argc <= 1 ){
        std::cerr << "Usage: "<<argv[0]<<" [infile]" << std::endl;
        return -1;
    }
    get_fasta_sequence_info(argv[1]);
 
} 

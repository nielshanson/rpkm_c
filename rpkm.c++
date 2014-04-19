#include "utilities.h"
#include "helper.h"

#include <map>

int main( int argc, char **argv ){
    // parse options
    Options options;
    if (argc < 9) { options.print_usage(argv[0]); exit(1); }
    if ( options.SetOptions(argc, argv)==false) { options.print_usage(argv[0]); exit(1); }

    //options.print_options();

    //get_fasta_sequence_info(options.contigs);
    map<string, unsigned int> orfnames;
    if( options.pathways_table.size()!=0 )
       read_orf_names(options.pathways_table, orfnames);

    map<string, CONTIG> contigs_dictionary;
    unsigned int genome_length = create_contigs_dictionary(options.contigs_file,  contigs_dictionary);
    std::cout << " Total genome length " << genome_length << std::endl;
    
 
} 

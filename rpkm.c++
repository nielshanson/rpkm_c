#include "utilities.h"
#include "helper.h"
#include <map>


bool compare_triplets(const TRIPLET &a, const TRIPLET &b) {
   return a.start < b.start ? true : false; 
}
 

int main( int argc, char **argv ){
    // parse options
    Options options;
    if (argc < 9) { options.print_usage(argv[0]); exit(1); }
    if ( options.SetOptions(argc, argv)==false) { options.print_usage(argv[0]); exit(1); }

    //options.print_options();

    //get_fasta_sequence_info(options.contigs);
    map<string, float> orfnames;
    if( options.pathways_table.size()!=0 )
       read_orf_names(options.pathways_table, orfnames);

    map<string, CONTIG> contigs_dictionary;
    unsigned int genome_length = create_contigs_dictionary(options.contigs_file,  contigs_dictionary);
   // std::cout << " Total genome length " << genome_length << std::endl;
   
    vector<MATCH> all_reads;
    all_reads.reserve(1000000);
    unsigned int num_mappable_reads =0;
    // creating the read multiplicity counts if there is a single end read file 
    unsigned int _se_num_mappable_reads =0;
    unsigned int _se_num_multireads =0;
    map<std::string, unsigned int> multireads;
    if( options.se_reads_map_file.size() != 0 ) {
        multireads.clear();

        all_reads.clear();
        std::cout << "\n\n" << "Searching for multireads from single read file " << options.se_reads_map_file << std::endl;
        _se_num_multireads = detect_multireads_blastoutput(options.se_reads_map_file, options.reads_map_file_format, all_reads,  multireads);

        if(!options.multi_reads) multireads.clear();
        _se_num_mappable_reads = process_blastoutput(options.se_reads_map_file, contigs_dictionary, options.reads_map_file_format, all_reads, multireads);
        num_mappable_reads +=  _se_num_mappable_reads;
    }
    
    // creating the read multiplicity counts if there is a paired ends read file 
    unsigned int _pe_num_mappable_reads =0;
    unsigned int _pe_num_multireads =0;
    if( options.pe_reads_map_file.size() != 0 ) {
        multireads.clear();
        all_reads.clear();
        std::cout << "\n\n" << "Searching for multireads from paired end reads file " << options.pe_reads_map_file << std::endl;
        _pe_num_multireads = detect_multireads_blastoutput(options.pe_reads_map_file, options.reads_map_file_format, all_reads,  multireads, true);
        if(!options.multi_reads) multireads.clear();
        _pe_num_mappable_reads = process_blastoutput(options.pe_reads_map_file, contigs_dictionary, options.reads_map_file_format, all_reads,  multireads);
        num_mappable_reads +=  _pe_num_mappable_reads;
         
    }
    
    std::cout << "\n\nSorting  the read matches .....";
    for(map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
       std::sort(it->second.M.begin(), it->second.M.end(), compare_triplets);
    /*    for(TRIPLETS::iterator it1= it->second.M.begin(); it1 != it->second.M.end(); it1++) {
          std::cout << it1->start << std::endl;
        }
*/
    }
    std::cout << "done\n";

    unsigned int total_covered_length = 0;
    unsigned int total_contig_length = 0;
    COVERAGE coverage;
    for(map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
       substring_coverage(contigs_dictionary, it->first, 1, it->second.L, coverage);
       total_covered_length += coverage.coverage*contigs_dictionary[it->first].L;
       total_contig_length += contigs_dictionary[it->first].L;
       //std::cout << coverage.coverage << "  " << coverage.numreads << "  " << coverage.substring_length << "   " << coverage.uncovered_length << std::endl;
    }
    
    map<string, float> _all_orfnames;
    unsigned int orf_length = 0;
    unsigned int _num_orfs = ORFWise_coverage(contigs_dictionary, options.orf_file, _all_orfnames, genome_length,  orf_length, num_mappable_reads);
    std::cout << "done computing orfwise coverage " << std::endl;

/*
    for(std::map<string, float>::iterator it = orfnames.begin(); it != orfnames.end(); it++) {
       std::cout << it->first << "  " << it->second << std::endl;
    }

*/
    if( options.pathways_table.size() > 0 )
        add_RPKM_value_to_pathway_table(options.pathways_table, options.output_file, orfnames);
    

   // average rpkm
    float _rpkm_sample = 0 ;
    float _count = 0 ;
    for(map<string, float>::iterator it= _all_orfnames.begin(); it != _all_orfnames.end(); it++)  {
        _rpkm_sample += it->second;
        _count += 1;
    }
    _count = (_count > 0) ? _count : 1;
    float _avg_rpkm_sample  = _rpkm_sample/_count;


    // rpkm average for the table
    float _rpkm_table = 0 ;
    _count = 0 ;
    for(map<string, float>::iterator it= orfnames.begin(); it != orfnames.end(); it++)  {
        _rpkm_table +=  _all_orfnames[it->first];
        _count += 1;
    }
    _count = (_count > 0) ? _count : 1;
    float _avg_rpkm_table = _rpkm_table/_count;

    _count = (float)orfnames.size();

    unsigned int _num_single_reads  = num_mappable_reads - _se_num_multireads - _pe_num_multireads;

    char buf[100000];

    std::cout << std::endl;
    sprintf(buf, "Number of Contigs               : %d ", (int)contigs_dictionary.size() );
    std::cout << buf<< std::endl;
    sprintf(buf, "Number of ORFs in sample        : %d ",_num_orfs);
    std::cout << buf  << std::endl;
    sprintf(buf, "Number of ORFs in pathway table : %d ",(int)_count); 
    std::cout << buf  << std::endl;
    sprintf(buf,"Perentage contig coverage       : %-5.2f%%", (float)total_covered_length/(float)total_contig_length);
    std::cout << buf  << std::endl;
    sprintf(buf,"Total Genome Length             : %d ",genome_length);
    std::cout << buf  << std::endl;
    sprintf(buf,"Total ORF Length                : %d ",orf_length);
    std::cout << buf  << std::endl;
    sprintf(buf,"Total num of mappable reads     : %d ",num_mappable_reads);
    std::cout << buf  << std::endl;
    sprintf(buf,"Total num of se mappable reads  : %d ",_se_num_mappable_reads);
    std::cout << buf  << std::endl;
    sprintf(buf,"Total num of pe mappable reads  : %d ",_pe_num_mappable_reads);
    std::cout << buf  << std::endl;
    sprintf(buf,"Number of single reads          : %d ",_num_single_reads);
    std::cout << buf  << std::endl;
    sprintf(buf,"Number of se multireads         : %d ",_se_num_multireads);
    std::cout << buf  << std::endl;
    sprintf(buf,"Number of pe multireads         : %d ",_pe_num_multireads);
    std::cout << buf  << std::endl;
    _se_num_mappable_reads = _se_num_mappable_reads > 0 ? _se_num_mappable_reads : 1; 
    sprintf(buf,"Percentage se multireads        : %-5.2f%% ",(float)_se_num_multireads*100/(float)_se_num_mappable_reads);
    std::cout << buf  << std::endl;
    _pe_num_mappable_reads = _pe_num_mappable_reads > 0 ? _pe_num_mappable_reads : 1; 
    sprintf(buf,"Percentage pe multireads        : %-5.2f%% ",(float)_pe_num_multireads*100/(float)_pe_num_mappable_reads);
    std::cout << buf  << std::endl;
    sprintf(buf,"Avg rpkm across ORFs in sample  : %.2f ",_avg_rpkm_sample);
    std::cout << buf  << std::endl;
    sprintf(buf,"Avg rpkm across ORFs pwy table  : %.2f ",_avg_rpkm_table);
    std::cout << buf  << std::endl;
 
} 

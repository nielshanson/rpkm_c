#include <map>
#include "types.h"
#include "utilities.h"
#include "helper.h"


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
    unsigned long genome_length = create_contigs_dictionary(options.contigs_file,  contigs_dictionary);
/*
    std::cout << " Total genome length " << genome_length << std::endl;
    std::cout << "number of contigs " << contigs_dictionary.size() << std::endl;
*/
    bool print_stats_file = false; 

    std::ostream *output;
    std::ofstream fout;
    if( options.stats_file.size() > 0) {
        print_stats_file = true; 
        fout.open(options.stats_file.c_str(), std::ifstream::out);
        output = &fout;
    }   
 
    vector<MATCH> all_reads;
    all_reads.reserve(80000000);
    //unsigned long num_mappable_reads =0;
    // creating the read multiplicity counts if there is a single end read file 
    map<std::string, unsigned long> multireads;
    
   
    RUN_STATS stats;
    RUN_STATS _stats;
    for( vector<string>::iterator it = options.read_map_files.begin() ; it != options.read_map_files.end(); it++)  { 
        if( it->size() == 0 ) continue;
        all_reads.clear();
        if( options.show_status ) 
          std::cout << "\n\n" << "Searching for multireads from reads file " << *it << std::endl;
        _stats = detect_multireads_blastoutput(*it, options.reads_map_file_format, all_reads, multireads, options.show_status);
        
        if( options.show_status ) _stats.printStats(&std::cout);

        if(print_stats_file) {
            *output << "\nStats for file :  " << *it << std::endl;
            _stats.printStats(output);
        }

        stats = stats + _stats;
        process_blastoutput(*it, contigs_dictionary, options.reads_map_file_format, all_reads,  multireads, options.show_status);
//        std::cout << "Number of mappable reads " << all_reads.size() << std::endl;
    }
  
   
    if( options.show_status ) 
       std::cout << "Composite stats for all files " << std::endl;

    if( print_stats_file)
       *output << "\nComposite stats for all files " << std::endl;
 
    if( options.show_status ) stats.printStats(&std::cout);

    if(print_stats_file) stats.printStats(output);
 

    if( options.show_status ) std::cout << "\n\nSorting  the read matches .....";
    for(map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
       std::sort(it->second.M.begin(), it->second.M.end(), compare_triplets);
    /*    for(TRIPLETS::iterator it1= it->second.M.begin(); it1 != it->second.M.end(); it1++) {
          std::cout << it1->start << std::endl;
        }
*/
    }
    std::cout << "done\n";

    unsigned long total_covered_length = 0;
    unsigned long total_contig_length = 0;
    COVERAGE coverage;
    for(map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
       substring_coverage(contigs_dictionary, it->first, 1, it->second.L, coverage);
       
       total_covered_length += coverage.coverage*contigs_dictionary[it->first].L;
       total_contig_length += contigs_dictionary[it->first].L;
       //std::cout << coverage.coverage << "  " << coverage.numreads << "  " << coverage.substring_length << "   " << coverage.uncovered_length << std::endl;
    }
    
    map<string, float> _all_orfnames;
    float total_distinct_reads = static_cast<double>(stats.num_distinct_reads_mapped) + static_cast<double>(stats.num_distinct_reads_unmapped);
    unsigned long orf_length = 0;
    unsigned long _num_orfs = ORFWise_coverage(contigs_dictionary, options.orf_file, _all_orfnames,\
                               genome_length,  orf_length, total_distinct_reads);
 //   std::cout << "done computing orfwise coverage " << std::endl;

/*
    for(std::map<string, float>::iterator it = orfnames.begin(); it != orfnames.end(); it++) {
       std::cout << it->first << "  " << it->second << std::endl;
    }

*/

    if( options.pathways_table.size() > 0 && options.output_file.size() > 0 )
        add_RPKM_value_to_pathway_table(options.pathways_table, options.output_file, _all_orfnames);
    

   // average rpkm
    float _rpkm_sample = 0 ;
    float _count = 0 ;
    for(map<string, float>::iterator it= _all_orfnames.begin(); it != _all_orfnames.end(); it++)  {
        _rpkm_sample += it->second;
       // std::cout << _rpkm_sample << std::endl;
        _count += 1;
    }
    _count = (_count > 0) ? _count : 1;
    float _avg_rpkm_sample  = _rpkm_sample/_count;



     // Now print out the ORFwise rpkm value
    if( options.orf_rpkm_file.size() > 0 )
        writeOut_ORFwise_RPKM_values(options.orf_rpkm_file, _all_orfnames);



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


    char buf[100000];

    std::cout << std::endl;
    sprintf(buf, "Number of Contigs               : %ld ", (long int)contigs_dictionary.size() );
    std::cout << buf<< std::endl;

    sprintf(buf, "Number of ORFs in sample        : %ld ",_num_orfs);
    std::cout << buf  << std::endl;

    sprintf(buf, "Number of ORFs in pathway table : %ld ",(long int)_count); 
    std::cout << buf  << std::endl;

    sprintf(buf,"Total contig cover length       : %ld ", total_covered_length/100);
    std::cout << buf  << std::endl;

    sprintf(buf,"Total Contig Length             : %ld ",total_contig_length);
    std::cout << buf  << std::endl;

    sprintf(buf,"Total Genome Length             : %ld ",genome_length);
    std::cout << buf  << std::endl;

    sprintf(buf,"Perentage contig coverage       : %-5.2f%%", (float)total_covered_length/(float)total_contig_length);
    std::cout << buf  << std::endl;

    sprintf(buf,"Total ORF Length                : %ld ",orf_length);
    std::cout << buf  << std::endl;

    sprintf(buf,"Total num of mappable reads     : %d ",stats.num_distinct_reads_mapped);
    std::cout << buf  << std::endl;

    sprintf(buf,"Total num of unmappable reads   : %d ",stats.num_distinct_reads_unmapped) ;
    std::cout << buf  << std::endl;

    sprintf(buf,"Number of total reads           : %d ",stats.num_distinct_reads_mapped + stats.num_distinct_reads_unmapped);
    std::cout << buf  << std::endl;

    total_distinct_reads = total_distinct_reads ==0 ? 1 : total_distinct_reads; 

    sprintf(buf,"Percentage of mapped reads      : %-5.2f ",\
          100*(static_cast<double>(stats.num_distinct_reads_mapped)/total_distinct_reads )) ;
    std::cout << buf  << std::endl;

    sprintf(buf,"Number of multireads            : %d ", stats.num_multireads);
    std::cout << buf  << std::endl;

    sprintf(buf,"Percentage of multireads        : %-5.2f%% ",100*(float)stats.num_multireads/total_distinct_reads);
    std::cout << buf  << std::endl;

    sprintf(buf,"Avg rpkm across ORFs in sample  : %.2f ",_avg_rpkm_sample);
    std::cout << buf  << std::endl;

    sprintf(buf,"Avg rpkm across ORFs pwy table  : %.2f ",_avg_rpkm_table);
    std::cout << buf  << std::endl;
 

   if( print_stats_file ) {
        *output << std::endl;
        *output << "Mapping stats for the sample " << std::endl;
        sprintf(buf, "Number of Contigs               : %ld ", (long int)contigs_dictionary.size() );
        *output << buf<< std::endl;
    
        sprintf(buf, "Number of ORFs in sample        : %ld ",_num_orfs);
        *output << buf  << std::endl;
    
        sprintf(buf, "Number of ORFs in pathway table : %ld ",(long int)_count); 
        *output << buf  << std::endl;
    
        sprintf(buf,"Total contig cover length       : %ld ", total_covered_length/100);
        *output << buf  << std::endl;
    
        sprintf(buf,"Total Contig Length             : %ld ",total_contig_length);
        *output << buf  << std::endl;
    
        sprintf(buf,"Total Genome Length             : %ld ",genome_length);
        *output << buf  << std::endl;
    
        sprintf(buf,"Perentage contig coverage       : %-5.2f%%", (float)total_covered_length/(float)total_contig_length);
        *output << buf  << std::endl;
    
        sprintf(buf,"Total ORF Length                : %ld ",orf_length);
        *output << buf  << std::endl;
    
        sprintf(buf,"Total num of mappable reads     : %d ",stats.num_distinct_reads_mapped);
        *output << buf  << std::endl;
    
        sprintf(buf,"Total num of unmappable reads   : %d ",stats.num_distinct_reads_unmapped) ;
        *output << buf  << std::endl;
    
        sprintf(buf,"Number of total reads           : %d ",stats.num_distinct_reads_mapped + stats.num_distinct_reads_unmapped);
        *output << buf  << std::endl;
    
        total_distinct_reads = total_distinct_reads ==0 ? 1 : total_distinct_reads; 
    
        sprintf(buf,"Percentage of mapped reads      : %-5.2f ",\
              100*(static_cast<double>(stats.num_distinct_reads_mapped)/total_distinct_reads )) ;
        *output << buf  << std::endl;
    
        sprintf(buf,"Number of multireads            : %d ", stats.num_multireads);
        *output << buf  << std::endl;
    
        sprintf(buf,"Percentage of multireads        : %-5.2f%% ",100*(float)stats.num_multireads/total_distinct_reads);
        *output << buf  << std::endl;
    
        sprintf(buf,"Avg rpkm across ORFs in sample  : %.2f ",_avg_rpkm_sample);
        *output << buf  << std::endl;
    
        sprintf(buf,"Avg rpkm across ORFs pwy table  : %.2f ",_avg_rpkm_table);
        *output << buf  << std::endl;

        fout.close();
 }




} 

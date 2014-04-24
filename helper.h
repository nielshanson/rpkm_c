#ifndef _HELPER
#define _HELPER
#include <map>
#include <string>
#include <ostream>
#include <iterator>
#include "types.h"
#include "fastareader.h"
#include "matchoutputparser.h"

using namespace std;
void read_orf_names(std::string pathways_table, std::map<string, float> &orfnames) ;
unsigned int create_contigs_dictionary(std::string contigs_file,  std::map<std::string, CONTIG> &contigs_dictionary);

unsigned int detect_multireads_blastoutput(const std::string &blastoutput_file, const std::string &format,\
      vector<MATCH> &all_reads, map<std::string, unsigned int> &multireads,  bool paired_reads=false);

unsigned int process_blastoutput(const std::string & reads_map_file, std::map<string, CONTIG> &contigs_dictionary,\
     const std::string &reads_map_file_format, vector<MATCH> &all_reads, map<string, unsigned int> & multireads);

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,\
                         unsigned int start, unsigned int end, COVERAGE &coverage);

unsigned int ORFWise_coverage( map<string, CONTIG> &contigs_dictionary, const string &orf_file,\
                              map<string, float> &orfnames, unsigned int genome_length,\
                              unsigned int &orf_length,  unsigned int num_mappable_reads);


void add_RPKM_value_to_pathway_table(const string &pathways_table, const string &output_file, map<string, float> &orfnames);


#endif //_HELPER

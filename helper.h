#ifndef _HELPER
#define _HELPER
#include <map>
#include <string>
#include "types.h"
#include "fastareader.h"
#include <iterator>

using namespace std;
void read_orf_names(std::string pathways_table, std::map<string, unsigned int> &orfnames) ;
unsigned int create_contigs_dictionary(std::string contigs_file,  std::map<std::string, CONTIG> &contigs_dictionary);






#endif //_HELPER

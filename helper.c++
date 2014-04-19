#include "helper.h"
using namespace std;



void read_orf_names(string pathways_table_filename, map<string, unsigned int> &orfnames) {

    char buf[1000000];
    std::ifstream input(pathways_table_filename.c_str());
    if(!input.good()){
        std::cerr << "Error opening '"<<pathways_table_filename<<"'. Bailing out." << std::endl;
        return ;
    }   

    string stringCOMMENT("PWY_NAME");
    vector<char *> fields; 

    std::string line;
    while( std::getline( input, line ).good() ){
       if( matchString(line, stringCOMMENT, true) ) continue;
       fields.clear();
       split(line, fields, buf); 
       if( fields.size() <= 5) continue; 

       for(vector<char *>::iterator it=fields.begin()+5; it!=fields.end() ; it++) {
           orfnames[std::string(*it)];
       }
    }
    std::cout << "num orfnames " << orfnames.size() << std::endl ;

    input.close();

}





unsigned int create_contigs_dictionary(std::string contigs_file,  std::map<std::string, CONTIG> &contigs_dictionary) {
    
     FastaReader fastareader(contigs_file);
     map<string, unsigned int> contig_lengths;
     map<string, unsigned int>::iterator  it_contig_lens;
     
     fastareader.get_fasta_sequence_info(contig_lengths);
     unsigned int genome_length = 0;
     CONTIG contig;
     for(it_contig_lens = contig_lengths.begin(); it_contig_lens != contig_lengths.end(); it_contig_lens++ ) {
       // std::cout << it_contig_lens->first <<  "   " << it_contig_lens->second << std::endl;
        genome_length +=  it_contig_lens->second;
        contig.L = it_contig_lens->second;
        contigs_dictionary[it_contig_lens->first] = contig; 
     }

   //  map<string, CONTIG>::iterator  it_contig_dict;

//     for(it_contig_dict = contigs_dictionary.begin(); it_contig_dict != contigs_dictionary.end(); it_contig_dict++ ) {
 //         std::cout << it_contig_dict->first <<  "   " << it_contig_dict->second.L << std::endl;
  //   }

     return genome_length;
}

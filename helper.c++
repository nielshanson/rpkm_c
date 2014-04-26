#include "helper.h"
using namespace std;


#define _MAX 100000000

void read_orf_names(string pathways_table_filename, map<string, float> &orfnames) {

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
           orfnames[std::string(*it)] = 0;
       }
    }

    input.close();

}


unsigned long create_contigs_dictionary(std::string contigs_file,  std::map<std::string, CONTIG> &contigs_dictionary) {
    
     FastaReader fastareader(contigs_file);
     map<string, unsigned long> contig_lengths;
     map<string, unsigned long>::iterator  it_contig_lens;
     
     fastareader.get_fasta_sequence_info(contig_lengths);
     unsigned long genome_length = 0;
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


unsigned long detect_multireads_blastoutput(const std::string &blastoutput_file, const std::string &format,\
     vector<MATCH> &all_reads, map<std::string, unsigned long> &multireads, bool paired_reads) {

    MatchOutputParser *parser = ParserFactory::createParser(blastoutput_file, format);
    if( parser ==0 ) {
        std::cout << "ERROR : Cannot open a parser to parse the file " << blastoutput_file << std::endl;
    }

    MATCH  match;
    map<std::string, unsigned long> _multireads;
    std::cout << std::endl << "Number of reads processed : " ;
    for(int i =0; ; i++ )  {
       if( !parser->nextline(match) )  break;
       if( i >= _MAX ) break;
        all_reads.push_back(match);
       if(i%10000==0) {
           std::cout << "\n\033[F\033[J";
           std::cout << i ;
       }
       //std::cout << match.query << "   " << match.subject <<  " "  << match.start << " " << match.end << std::endl;
       if( _multireads.find(match.query)==_multireads.end() ) _multireads[match.query] = 0;
       _multireads[match.query] += 1;
    }

    // now store the multireads into the multireads map variable
    unsigned long count;
    for( map<std::string, unsigned long>::iterator it = _multireads.begin(); it != _multireads.end(); it++) {
       if( !paired_reads &&  it->second > 1)  
          multireads[it->first] = it->second;

       if( paired_reads &&  it->second > 2)   {
          count = static_cast<int>( it->second/2);
          multireads[it->first] = count + (it->second%2);
       }
    }
    std::cout << std::endl << "Number of multireads       : " << multireads.size() << std::endl; 
    delete parser;
    return multireads.size();
}


unsigned long process_blastoutput(const std::string & reads_map_file, std::map<string, CONTIG> &contigs_dictionary,\
      const std::string &reads_map_file_format, std::vector<MATCH> &all_reads,   std::map<string, unsigned long> & multireads) {

    MATCH  match;
    std::map<string, bool> read_map;
    //unsigned long count=0;
    //unsigned long length = 0;
    unsigned long read_multiplicity ;

 //   std::map<string, TRIPLETS> temp_contig_dictionary;
    TRIPLETS triplets;
    TRIPLET triplet;

    /*for( map<std::string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
        temp_contig_dictionary[it->first] = triplets;
    }
*/
    
    int i =0;
    std::cout << "Number of hits processed : " ;
    for(vector<MATCH>::iterator it=all_reads.begin();  it!= all_reads.end(); it++ )  {

       if( i >=_MAX ) break;

       if(i%10000==0) {
           std::cout << "\n\033[F\033[J";
           std::cout << i ;
       }
       i++;

       //std::cout << match.query << std::endl;
       //if( read_map.find(match.query)== read_map.begin() ) {
        // std::cout << match.query << std::endl;
       read_map[it->query] = true;
       //}

       if( contigs_dictionary.find(it->subject)==contigs_dictionary.end() ) {
          std::cout << " Missing contig " << it->subject << std::endl;
          std::cerr << "ERROR : Could not find the matched contig in the contig file " << std::endl;
          exit(1);
       }

       read_multiplicity =  (multireads.find(it->query) == multireads.end()) ? 1 : multireads[it->query];

       if(  it->start < it->end ) {
           triplet.start = it->start;
           triplet.end = it->end;
       }
       else {
           triplet.start = it->end;
           triplet.end = it->start;
       }
       triplet.multi = read_multiplicity;
       contigs_dictionary[it->subject].M.push_back(triplet);
    }
    std::cout << i << std::endl;
    std::cout << "Number of mappable reads : " << read_map.size() << std::endl;
    
    return read_map.size();
}


void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,\
                          unsigned long start, unsigned long end, COVERAGE &coverage) {


      if( contigs_dictionary.find(contig) == contigs_dictionary.end()   || contigs_dictionary[contig].L==0 ) {
          coverage.coverage =0 ;
          coverage.numreads =0 ;
          coverage.substring_length = end - start ;
          coverage.uncovered_length = 0;
      }

      float numreads =0;
      float uncovered_length = 0;
      unsigned long p_end = start;
      float _coverage = 0;
      
     for(std::vector<TRIPLET>::iterator it = contigs_dictionary[contig].M.begin(); it != contigs_dictionary[contig].M.end(); it++) {
         uncovered_length  +=  ( p_end > it->start  || it->start > end) ? 0 : it->start - p_end;
         if( it->end > p_end ) p_end = it->end;  //make sure the read start and end are not going beyoing the contig
          
         if( (start <= it->start && it->start <= end) ||  (start <= it->end && it->end <= end)  )
            numreads += 1/static_cast<float>(it->multi);
     }
     uncovered_length += (p_end > end ) ? 0 : (end - p_end);

     unsigned long substring_length = end - start;
     if( substring_length > 0 ) 
        _coverage = ((float)(substring_length - uncovered_length )/(float)substring_length)*100;

     coverage.numreads = numreads;
     coverage.coverage = _coverage;
     coverage.substring_length  = substring_length;
     coverage.uncovered_length =  uncovered_length;
}


unsigned long ORFWise_coverage( map<string, CONTIG> &contigs_dictionary, const string &orf_file,\
                               map<string, float> &orfnames, unsigned long genome_length,\
                               unsigned long &orf_length,  unsigned long num_mappable_reads) {

    MATCH  match;
    COVERAGE coverage;
    MatchOutputParser *parser = ParserFactory::createParser(orf_file, std::string("GFF"));

    std::map<string, vector<MATCH> > orf_dictionary;

    unsigned long _num_orfs=0;

    vector<MATCH> match_vector;

    for(int i =0; ; i++ )  {
       if( !parser->nextline(match) )  break;
       if(i%10000==0) {
           std::cout << "\n\033[F\033[J";
           std::cout << i ;
       }

       if( orf_dictionary.find(match.query) == orf_dictionary.end() ) 
          orf_dictionary[match.query] = match_vector;

       orf_dictionary[match.query].push_back(match);
       _num_orfs += 1;
     //  if( orfnames.find(match.subject) != orfnames.end() ) {
           try {
               substring_coverage(contigs_dictionary, match.query, match.start, match.end, coverage); 
               orf_length += match.end - match.start;
         //   std::cout << coverage.coverage << "  " << coverage.numreads << "  " << coverage.substring_length << "   " << coverage.uncovered_length << std::endl;
               orfnames[match.subject] = (1E9/static_cast<float>(num_mappable_reads))*(static_cast<float>(coverage.numreads)/static_cast<float>(coverage.substring_length));
           }
           catch(...) {
               std::cout << "error\n";
               orfnames[match.subject] = 0;
           }
      // }

       //std::cout << match.query << " " << match.subject << "  " << match.start << " " << match.end << std::endl;

    }

     return _num_orfs;
}

void add_RPKM_value_to_pathway_table(const string &pathways_table_filename, const string &output_file, map<string, float> &orfnames) {
    char buf[1000000];
    std::ifstream input(pathways_table_filename.c_str());
    if(!input.good()){
        std::cerr << "Error opening '"<<pathways_table_filename<<"'. Bailing out." << std::endl;
        return ;
    }   

    string stringCOMMENT("PWY_NAME");
/*
    vector<char *> fields; 

       fields.clear();
       split(line, fields, buf); 
       if( fields.size() <= 5) continue; 
*/
    vector<string> lines;
    std::string line;

    //read in all the lines
    string headerline;
    while( std::getline( input, line ).good() ){
       if( matchString(line, stringCOMMENT, true) ) { headerline = line;  continue;}
       lines.push_back(line);
    }
    input.close();

    std::ostream *output  = &std::cout;
    std::ofstream fout;
    if( output_file.size()!=0) {
        fout.open(output_file.c_str(), std::ifstream::out);
        output = &fout;
    }

    
    float pwy_rpkm;
    vector<char *> fields; 
    *output << headerline << "\tRPKM_COUNT" << std::endl ;
    for(vector<string>::iterator it=lines.begin(); it != lines.end(); it++) {
        split(*it, fields, buf, '\t'); 
        if( fields.size() < 5) continue; 
         
        for(int i =0; i< 5; i++)  {
            if( i > 0 ) *output << '\t';  
            *output << fields[i]  ;
        }

       
       pwy_rpkm = 0;
       for(vector<char *>::iterator it=fields.begin()+5; it!=fields.end() ; it++) {
          pwy_rpkm += orfnames[std::string(*it)];
       }
       sprintf(buf, "\t%0.2f", pwy_rpkm) ;
       *output << buf;

       for(vector<char *>::iterator it=fields.begin()+5; it!=fields.end() ; it++) {
          *output << "\t" << *it;
       }
       *output << std::endl ;
    }
    if( output_file.size() !=0 )  fout.close();
    
}

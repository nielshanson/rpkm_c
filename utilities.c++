#include "utilities.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
   std::cout << "USAGE : "   << arg\
             << "      : -h [ for help ]\n"\
             << "      : --c <contigs_file>  [REQUIRED]\n"\
             << "      : --r  <read_map_file_sam> \n"\
             << "      : --ORFS  <orf_file_gff> [REQUIRED]\n"
             << "      : --ORF-RPKM  <orf_rpkm_file> [default: None]\n"\
             << "      : --stats stats_file [OPTIONAL]\n"\
             << "      : --o <outputfile> \n"\
             << "      : --m <multireads>  [OPTIONAL , default on]\n"\
             << "      : --status [shows status]"\
             << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) { 
   for(int i = 1; i < argc ; i++) {   
       if( strncmp(argv[i], "--c", strlen("--c")) == 0 ) {   
          this->contigs_file = argv[++i];
       }   
       else if( strncmp(argv[i], "-h", strlen("-h")) == 0 ) {   
          print_usage(argv[0]);
          exit(0);
       }   
       else if( strncmp(argv[i], "--stats", strlen("--stats")) == 0 ) {   
          this->stats_file = argv[++i];
       }   
       else if(strncmp(argv[i], "-r", strlen("-r")) == 0 ) {   
          this->read_map_files.push_back(string(argv[++i]));
       }   
       else if(strncmp(argv[i], "--status", strlen("--status")) == 0 ) {   
    //      std::cout << "Status is true";
          this->show_status = true;
       }   
       else if( strncmp(argv[i], "--o", strlen("--o")) == 0 ) {   
          this->output_file = argv[++i];
       }
       else if( strncmp(argv[i], "--ORFS", strlen("--ORFS")) == 0 ) {   
          this->orf_file = argv[++i];
       }
       else if( strncmp(argv[i], "--ORF-RPKM", strlen("--ORF-RPKM")) == 0 ) {   
          this->orf_rpkm_file = argv[++i];
       }
       else if( strncmp(argv[i], "--m", strlen("--m")) == 0  ) {   
     //     std::cout << "Multireads is true";
          this->multi_reads = true;
       }
       else if( (strncmp(argv[i], "--f", strlen("--f")) == 0 ) || ( strncmp(argv[i], "--format", strlen("--format")) == 0) ) {   
          if( (strncmp(argv[i+1], "blastout", strlen("blastout")) == 0 ) ||\
              (strncmp(argv[i+1], "sam-1", strlen("sam-1")) == 0) ||\
              (strncmp(argv[i+1], "sam-2", strlen("sam-2")) == 0) ) {   
              this->reads_map_file_format = argv[++i];
          }
          else {
               std::cout << "ERROR: Choices for -format option must of type blastout, sam-1 or sam-2" << std::endl;
               return false;
          }
       }
       else if( strncmp(argv[i], "--p", strlen("--p")) == 0  ) {   
          this->pathways_table = argv[++i];
       }
    } //for loop for arguments processing

    if( this->contigs_file.size()==0) {
       std::cout << "ERROR: There must be a contigs file" << std::endl;
       return false;
    }

    if( this->read_map_files.size()==0 ) {
       std::cout << "ERROR: There must be a at least one read map file" << std::endl;
       return false;
    }
    
    if( this->orf_file.size()==0) {
       std::cout << "ERROR: There must be an ORF file [ GFF format ]" << std::endl;
       return false;
    }

/*
    if( this->output_file.size()==0) {
       std::cout << "ERROR: There must be a output file prefix" << std::endl;
       return false;
    }
*/
    return true;
    
};

void Options::print_options() {
       std::cout << "Contigs file (FASTA)              : "<< this->contigs_file << std::endl; 
       std::cout << "Read alignment file (SAM format)  : "<< this->read_map_files.size() << std::endl; 
       std::cout << "ORF file [GFF]                    : "<< this->orf_file << std::endl; 
       std::cout << "Pathways table                    : "<< this->pathways_table << std::endl; 
       std::cout << "Output file                       : "<< this->output_file << std::endl; 
       std::cout << "Multi read treatment              : "<< this->multi_reads << std::endl; 
}

 
//this function gets us in a map the name of a sequnce and length
// as in the provided fasta file
void get_fasta_sequence_info(const std::string & fasta_file_name) {

    std::ifstream input(fasta_file_name.c_str());
    if(!input.good()){
        std::cerr << "Error opening '"<<fasta_file_name<<"'. Bailing out." << std::endl;
        return ;
    }
 
    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                std::cout << extract_sequence_name(name) << " : " << content.size() << std::endl;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        name = extract_sequence_name(name);
        std::cout << name << " : " << content.size() << std::endl;
    }
    input.close();
 
}

std::string extract_sequence_name(const std::string &name) {
     char  cstr[1000000];
     std::strcpy(cstr, name.c_str());
     
     char * cptr = cstr;

     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++; 
     (*cptr) ='\0';

     std::string sname(cstr);
     return sname;
}



void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.push_back(s1);
  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       v.push_back(s1+1);
     }
     s1++;
  }
}

std::string get_orf_name(std::string  &strn, std::vector<char *> &v, char *buf) {
    split(strn, v, buf, ';'); 
    if(v.size() == 0)  return std::string("");
    split(std::string(v[0]), v, buf, '=');
    if(v.size() < 2)  return std::string("");
    return std::string(v[1]);
}
bool matchString(const string &str, const string & stringtomatch, bool fromstart) {

    unsigned long pos = str.find(stringtomatch);
    if(fromstart &&  pos ==0 ) return true;

    if( !fromstart && pos >= 0) return true;
    return false;

}

string to_string(unsigned long i) {
    char  c[100];
    char *p = c;
    int j = 0;
    char z = '0';

    while( i > 0 ) {
       if(i< 10) {
         *p='0' + i;
          p++;
          break;
       } 
       else {
           j = i%10;
           i = (i - j)/10;
           *p =  z + j;
           *p++;
        }
    }
    *p = '\0';
    p--;

    return string(c);
}

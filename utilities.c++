#include "utilities.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
   std::cout << "USAGE : " << arg\
             << "-c/--contigs contigs_file --r1/--reads-map1 reads_ma  p_file "\
             << "--r2/--reads-map2 (Only required for paired-end sam files) reads_map_file "\
             << "-O/--ORF  orf_file_gff -p/--pathways pathways_table -f/--format format [default: blastout]  [ -o outputfile ]"\
             << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) { 
   for(int i = 1; i < argc ; i++) {   
       if( strncmp(argv[i], "-c", strlen("-c")) == 0 ) {   
          this->contigs_file = argv[++i];
       }   
       else if( strncmp(argv[i], "--r1", strlen("--r1")) == 0 ) {   
          this->se_reads_map_file =  argv[++i];
       }   
       else if( strncmp(argv[i], "--r2", strlen("--r2")) == 0 ) {   
          this->pe_reads_map_file =  argv[++i];
       }
       else if( strncmp(argv[i], "-o", strlen("-o")) == 0 ) {   
          this->output_file = argv[++i];
       }
       else if( strncmp(argv[i], "-O", strlen("-O")) == 0 ) {   
          this->orf_file = argv[++i];
       }
       else if( (strncmp(argv[i], "-m", strlen("-m")) == 0 ) || ( strncmp(argv[i], "--multireads", strlen("--multireads")) == 0) ) {   
          this->multi_reads = true;
       }
       else if( (strncmp(argv[i], "-f", strlen("-f")) == 0 ) || ( strncmp(argv[i], "--format", strlen("--format")) == 0) ) {   
          if( (strncmp(argv[i+1], "blastout", strlen("blastout")) == 0 ) ||\
              (strncmp(argv[i+1], "sam-1", strlen("sam-1")) == 0) ||\
              (strncmp(argv[i+1], "sam-2", strlen("sam-2")) == 0) ) {   
              this->map_format = argv[++i];
          }
          else {
               std::cout << "ERROR: Choices for -format option must of type blastout, sam-1 or sam-2" << std::endl;
               return false;
          }
       }
       else if( (strncmp(argv[i], "-p", strlen("-p")) == 0 ) || ( strncmp(argv[i], "--pathways", strlen("--pathways")) == 0) ) {   
          this->pathways_table = argv[++i];
       }
    }

    if( this->contigs_file.size()==0) {
       std::cout << "ERROR: There must be a contigs file" << std::endl;
       return false;
    }

    if( this->se_reads_map_file.size()==0 &&  this->pe_reads_map_file.size()==0) {
       std::cout << "ERROR: There must be a at least one read map file" << std::endl;
       return false;
    }
    
    if( this->orf_file.size()==0) {
       std::cout << "ERROR: There must be an ORF file [ GFF format ]" << std::endl;
       return false;
    }

    if( this->output_file.size()==0) {
       std::cout << "ERROR: There must be a output file prefix" << std::endl;
       return false;
    }
    return true;
    
};

void Options::print_options() {
       std::cout << "Contigs file            : "<< this->contigs_file << std::endl; 
       std::cout << "Single end reads file   : "<< this->se_reads_map_file << std::endl; 
       std::cout << "Paired ends reads file  : "<< this->pe_reads_map_file << std::endl; 
       std::cout << "ORF file [GFF]          : "<< this->orf_file << std::endl; 
       std::cout << "Pathways table          : "<< this->pathways_table << std::endl; 
       std::cout << "Output file             : "<< this->output_file << std::endl; 
       std::cout << "File Format             : "<< this->map_format << std::endl; 
       std::cout << "Multi reads             : "<< this->multi_reads << std::endl; 
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

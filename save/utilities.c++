#include "utilities.h"
#include <iostream>
#include <fstream>
 
//this function gets us in a map the name of a sequnce and length
// as in the provided fasta file
void get_fasta_sequence_info(char *fasta_file_name) {

    std::ifstream input(fasta_file_name);
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

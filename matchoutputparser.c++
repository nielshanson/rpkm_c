#include "matchoutputparser.h"

using namespace std;

MatchOutputParser::MatchOutputParser(const std::string &filename, const std::string &format) {
     this->filename = filename;
     this->format = format;
     this->num_unmapped_reads =0;
};

unsigned long MatchOutputParser::get_Num_Unmapped_Reads() {
    return  this->num_unmapped_reads;
}

MatchOutputParser::~MatchOutputParser() {
}


SamFileParser::SamFileParser(const std::string &filename, const std::string &format):MatchOutputParser(filename, format) {
     this->input.open(filename.c_str(), std::ifstream::in);

     if(!this->input.good()){
         std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
         return ;
     }  
     this->count = 0;

}

SamFileParser::~SamFileParser() {
   this->input.close();
}

bool SamFileParser::getMateInfo(unsigned int i, MATCH &match)  {

    unsigned int a = i;
    bool singleton  = 0;
    //bool c = a&1;
   
     
    a = a >> 2;
    match.mapped = !(a&1); 
    singleton = a&1;

    a = a >> 1;
    match.orphan = a&1; 
    singleton = singleton^(a&1);

    a = a >> 3;
     
    if( a&1 )  {
         match.parity = 0; 
         a = a >> 1;
    }
    else {
         a = a >> 1;
         if( a&1 )   match.parity  = 1; 
         else return false;
    }
    a = a >> 1;
    match.multi = a&1; 
    
    a = a >> 3;
    match.chimeric = a&1; 
    match.singleton = singleton;

    return true;

}
bool SamFileParser::nextline(MATCH &match) {
     string line;
     std::string skipPattern("@");
     std::string skipStar("*");

     bool _success = false;
     while( std::getline(this->input, line ).good()) {
        // std::cout << line << std::endl;
         if(matchString(line, skipPattern, true) ) continue;

         fields.clear();
         split(line, fields, this->buf,'\t');

         if(fields.size() < 9)  continue;

/*
         if( matchString(std::string(fields[2]), skipStar, true)) { 
            this->num_unmapped_reads++; continue; 
         }
*/
         _success = true;
         break;
     }  
    
     if( _success )  {  
        //match.query = to_string(this->count++); 
        match.query =  fields[0]; 
        match.subject = std::string(fields[2]);
        match.start = atoi(fields[3]);
        match.end =  match.start +  std::string(fields[9]).size();
        bool status = getMateInfo(static_cast<unsigned int>(atoi(fields[1])), match);
        //if(status==false) return false;


        //std::cout << match.query << "  " << match.mapped << " "  << match.parity << " " << match.orphan << " " << match.chimeric << " " << match.multi << "  " << fields[1] << std::endl;

        return true;
     }
     
    return false;
}

GffFileParser::GffFileParser(const std::string &filename, const std::string &format):MatchOutputParser(filename, format) {
     this->input.open(filename.c_str(), std::ifstream::in);

     if(!this->input.good()){
         std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
         return ;
     }  

}

GffFileParser::~GffFileParser() {
   this->input.close();
}




bool GffFileParser::nextline(MATCH &match) {
     string line;
     string field8; 

     bool _success = false;
     
     while( std::getline(this->input, line ).good() ) {
         split(line, fields, this->buf, '\t');
         if(fields.size() < 9)  continue;
        // std::cout << fields[3] << " - " << fields[4] << " " << _success <<std::endl;
         _success = true;
         break;
     }  
    
     if( _success )  {  
        match.query = std::string(fields[0]); 
        match.start = atoi(fields[3]);
        match.end = atoi(fields[4]);
        //respect the order before calling get_orf_name
        // because the same buffer is used
        field8 = std::string(fields[8]);
        match.subject = get_orf_name(field8, this->tempv, this->buf);
        return true;
     }
     
    return false;
}

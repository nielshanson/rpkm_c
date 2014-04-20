#include "matchoutputparser.h"

using namespace std;

MatchOutputParser::MatchOutputParser(const std::string &filename, const std::string &format) {
     this->filename = filename;
     this->format = format;
};


MatchOutputParser::~MatchOutputParser() {
}


SamFileParser::SamFileParser(const std::string &filename, const std::string &format):MatchOutputParser(filename, format) {
     this->input.open(filename.c_str(), std::ifstream::in);

     if(!this->input.good()){
         std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
         return ;
     }  

}

SamFileParser::~SamFileParser() {
   this->input.close();
}

bool SamFileParser::nextline(MATCH &match) {
     string line;
     std::string skipPattern("@");
     std::string skipStar("*");

     bool _success = false;
     while( std::getline(this->input, line ).good()) {
         if(matchString(line, skipPattern, true) ) continue;

         fields.clear();
         split(line, fields, this->buf);

         if(fields.size() < 9)  continue;

         if( matchString(std::string(fields[2]), skipStar, true)) continue;
         _success = true;
         break;
     }  
    
     if( _success )  {  
        match.query = std::string(fields[0]); 
        match.subject = std::string(fields[2]);
        match.start = atoi(fields[3]);
        match.end =  match.start +  std::string(fields[9]).size();
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

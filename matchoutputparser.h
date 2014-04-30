#ifndef _MATHOUTPUTPARSER
#define _MATHOUTPUTPARSER
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "utilities.h"
#include "types.h"


using namespace std;


class MatchOutputParser {

public:
      std::string filename;
      std::string format;
      std::ifstream input;
      char buf[100000];
      vector<char *> tempv;
      virtual ~MatchOutputParser() = 0;
      vector<char *> fields;
      MatchOutputParser(const std::string &filename, const std::string &format);
      unsigned long get_Num_Unmapped_Reads();
      virtual bool nextline(MATCH &match)=0;
protected:
      unsigned long num_unmapped_reads;
};


//subclass of the MathOutputParser
class SamFileParser: virtual public MatchOutputParser {
private:
      unsigned long count;
public:
      SamFileParser(const std::string &filename, const std::string &format);
      virtual bool nextline(MATCH &match);
      ~SamFileParser();
};

//subclass of the MathOutputParser
class GffFileParser: virtual public MatchOutputParser {
public:
      GffFileParser(const std::string &filename, const std::string &format);
      virtual bool nextline(MATCH &match);
      ~GffFileParser();
};

class ParserFactory {
public:
     static MatchOutputParser * createParser(const std::string &filename, const std::string &format) {
            if(format.find("sam-1") != string::npos  || format.find("sam-2") !=string::npos )  {
                return new SamFileParser( filename, format);
            }

            if(format.find("GFF") != string::npos ) {
                return new GffFileParser( filename, format);
            }
            return 0;
     }
};



#endif //_MATHOUTPUTPARSER

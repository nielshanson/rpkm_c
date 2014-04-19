#ifndef _MATHOUTPUTPARSER
#define _MATHOUTPUTPARSER
#include <string>
#include "types.h"

using namespace std;

class MatchOutputParser {
private:
      std::string filename;
      std::string format;
      
public:
      MatchOutputParser(const std::string &filename, const std::string &format);
      virtual bool nextline(MATCH &match)=0;
};


//subclass of the MathOutputParser
class SamFileParser: public MatchOutputParser {

public:
      SamFileParser(const std::string &filename, const std::string &format);
      bool nextline(MATCH &match);
};




#endif //_MATHOUTPUTPARSER

#ifndef __RPKM_TYPE
#define __RPKM_TYPE
#include <map>
#include <vector>

using namespace std;

typedef struct _TRIPLET {
    unsigned int start, end, multi;
} TRIPLET;

typedef vector<TRIPLET> TRIPLETS;
typedef struct _CONTIG {
    unsigned int L;
    TRIPLETS M;
} CONTIG;


typedef  struct _MATCH {
    std::string query, subject;
    unsigned int start, end;

} MATCH;


#endif //__RPKM_TYPE




#ifndef __RPKM_TYPE
#define __RPKM_TYPE
#include <map>
#include <vector>
#include <iostream>

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
    bool parity; 
    bool mapped;
    bool orphan;
    bool multi;
    bool chimeric;
    bool singleton;
    float  w;
    _MATCH(): w(0) { } 
} MATCH;


template< typename A, typename B, typename C, typename D>
struct QUADRUPLE {
     A first;
     B second;
     C third;
     D fourth;
};


typedef struct _RUN_STATS {
    int num_unmapped_reads ; 
    int num_mapped_reads ; 
    int num_singleton_reads ; 
    int num_reads_1 ; 
    int num_reads_2 ; 
    int num_multireads ;
    int num_total_reads ;
    int num_secondary_hits ;
    int num_distinct_reads_unmapped;
    int num_distinct_reads_mapped;
    
    _RUN_STATS() {
        num_unmapped_reads =0; 
        num_mapped_reads =0; 
        num_singleton_reads =0; 
        num_reads_1 =0; 
        num_reads_2 =0; 
        num_multireads = 0;
        num_total_reads = 0;
        num_secondary_hits = 0;
        num_distinct_reads_unmapped =0 ;
        num_distinct_reads_mapped = 0;
    }

    struct _RUN_STATS&  operator+( const struct _RUN_STATS & stats) {
        this->num_unmapped_reads += stats.num_unmapped_reads;  
        this->num_mapped_reads  += stats.num_mapped_reads  ;
        num_singleton_reads += stats.num_singleton_reads ;
        num_reads_1 += stats.num_reads_1;  
        num_reads_2 += stats.num_reads_2; 
        num_multireads += stats.num_multireads; 
        num_total_reads  += stats.num_total_reads  ;
        num_secondary_hits  += stats.num_secondary_hits;
        num_distinct_reads_unmapped  += stats.num_distinct_reads_unmapped  ;
        num_distinct_reads_mapped  += stats.num_distinct_reads_mapped  ;
        return *this;
    }


    void printStats() {
        std::cout << std::endl;
        std::cout << "Number of lines in SAM file      : " << num_total_reads << std::endl; 
        std::cout << "Number of unmapped lines         : " << num_unmapped_reads << std::endl; 
        std::cout << "Number of mapped lines           : " << num_mapped_reads << std::endl; 
        std::cout << "Number of singletons             : " << num_singleton_reads << std::endl; 
        std::cout << "Number of read 1                 : " << num_reads_1  << std::endl; 
        std::cout << "Number of read 2                 : " << num_reads_2 << std::endl; 
        std::cout << "Number of multireads             : " << num_multireads << std::endl; 
        std::cout << "Number of secondary hits         : " << num_secondary_hits << std::endl; 
        std::cout << "Number distinct reads mapped     : " << num_distinct_reads_mapped << std::endl; 
        std::cout << "Number distinct reads unmapped   : " << num_distinct_reads_unmapped << std::endl; 
    }
} RUN_STATS;

typedef struct _COVERAGE {
   float coverage, numreads;
   unsigned int substring_length, uncovered_length;
} COVERAGE;
#endif //__RPKM_TYPE




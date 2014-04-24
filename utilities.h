#ifndef _UTILITIES
#define _UTILITIES
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

using namespace std;

// Structure for RPKM input options
struct Options {
   
    /* Input files for RPKM */
    string contigs_file; // the contigs file
    string se_reads_map_file; // single-end or BLAST alginments
    string pe_reads_map_file; // paired-end alignments
    string orf_file; // .gff ORF file from PRODIGAL
    string pathways_table; // a table from Pathway Tools with ORFs
    string output_file; // location to write output file (i.e., update pathway table)
    
    /* Flags and settings */
    bool multi_reads; // flag for detecting multiple mapping of reads
    string reads_map_file_format; // aligner type BWA or BLAST, two SAM files or one
    
    // Constructor with default settings 
	Options(){
       contigs_file = "";
       se_reads_map_file = "";
       pe_reads_map_file = "";
       orf_file = "";
       pathways_table = "";
       output_file = "";

       multi_reads = false;
       reads_map_file_format = "blastout";  
	};
    
    void print_usage( char *arg);
    void print_options();
 
    bool SetOptions(int argc, char *argv[]);
	// bool SetOptionCommon( const char *flag, const char *value );
	// bool SetOption( const char *flag, const char *value );
	// bool SetOption2D( const char *flag, const char *value );
	// bool SetOptionEST( const char *flag, const char *value );
	// bool SetOptions( int argc, char *argv[], bool twodata=false, bool est=false );

	void Validate();
	// void ComputeTableLimits( int naan, int typical_len, size_t mem_need );

	void Print();
};

void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');
std::string get_orf_name(std::string & strn, std::vector<char *> &v, char *buf);


bool matchString(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

std::string extract_sequence_name(const std::string &name);

#endif //_UTILITIES

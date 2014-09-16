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
    string orf_file; // .gff ORF file from PRODIGAL
    string pathways_table; // a table from Pathway Tools with ORFs
    string output_file; // location to write output file (i.e., update pathway table)
    string stats_file; // location to write output file (i.e., update pathway table)
    string orf_rpkm_file; // location to write output file (i.e., update pathway table)
    vector<string> read_map_files;
    
    /* Flags and settings */
    bool multi_reads; // flag for detecting multiple mapping of reads
    bool show_status; // shows the counter that countes the number of reads processed, and other info 
                       // on screen
    string reads_map_file_format; // aligner type BWA or BLAST, two SAM files or one
    
    // Constructor with default settings 
    Options(){
       contigs_file = "";
       stats_file = ""; // location to write output file (i.e., update pathway table)
       read_map_files.clear();
       orf_file = "";
       pathways_table = "";
       output_file = "";
       output_file = "";
       orf_rpkm_file = ""; // location to write output file (i.e., update pathway table)

       multi_reads = false;
       show_status = false;
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

string to_string(unsigned long i);

#endif //_UTILITIES


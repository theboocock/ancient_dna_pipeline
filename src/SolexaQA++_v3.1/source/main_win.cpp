/*
Program: SolexaQA++ v.3.1
Calculates quality statistics on Illumina FASTQ sequence files
and creates visual representations of run quality
Murray Cox, Patrick Biggs, Daniel Peterson and Mauro Truglio
Massey University, New Zealand
Email contact <m.truglio@massey.ac.nz>
July 2014

Changelog:

-Version 3.1: corrected a file naming bug when operating on .gz files (Linux and Mac only); improved command line output; added warning message for SRA-generated fastq files not containing the tile number in the header;

  Analysis:     no other changes;
  Dynamictrim:  added the [-a|--anchor] option, which allows trimming from the 3' end only;
  Lengthsort:   New naming system for paired end reads:
                    Read 1 are named "[read1 name].paired"
                    Read 2 are named "[read2 name].paired"
                    Single reads are named "[read1 name].single"
                    Discarded reads are named "[read1 name].discard";

- Version 3.0: Complete rewrite of the algorithm in C++.
  The three components of the SolexaQA v2.x suite (SolexaQA.pl, DynamicTrim.pl and LengthSort.pl)
  are now three options ("analysis", "dynamictrim" and "lengthsort" respectively) of the main program.

  Analysis: general performance improvement; it can now process fastq files with variable
            read lengths and is compatible with IonTorrent files (-t option) and compressed files (.gz) on Linux and Mac
  DynamicTrim: general performance improvement (up to 100% faster); bug fix for truncated graphs on
               variable length reads; supports compressed files (.gz) on Linux and Mac
  LengthSort: paired-end mode includes a new option (-c) to remove non-matching reads from the two fastq
              files before processing them; supports compressed files (.gz) on Linux and Mac; new histogram output.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <string.h>
#include <tuple>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <getopt.h>
#include "boost/regex.hpp"
#include <vector>
#include "boost/assign.hpp"
#include "boost/filesystem.hpp"
#include <libgen.h>
#include <map>
#include <typeinfo>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "R_codes.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;
using std::vector;
using std::string;
using std::cin;
using std::cout;
using std::endl;


using namespace boost::filesystem;
using namespace boost::assign;
using namespace boost;

using namespace std;


static void show_usage(std::string name, std::string type){
        if(type=="general"){
            std::cerr << "\nUsage: " << name << " <command> [options]\n\n"<<"Command: analysis\tquality analysis and graphs generation\n"
                            "\t dynamictrim\ttrim reads using a chosen threshold\n"
                            "\t lengthsort\tsort reads by a chosen length\n\n";

        }
        else if(type=="analysis"){
            std::cerr << "\nUsage: " << name << " analysis input_files [-p|probcutoff 0.05] [-h|phredcutoff 13] [-v|variance] [-m|minmax] [-s|sample 10000] [-b|bwa] [-d|directory path] [--sanger --solexa --illumina] [-t|torrent]\n"
                        << "\nOptions:\n"
                        << "-p|--probcutoff     probability value (between 0 and 1) at which base-calling error is considered too high (default; p = 0.05) *or*\n"
                        << "-h|--phredcutoff    Phred quality score (between 0 and 41) at which base-calling error is considered too high\n"
                        << "-v|--variance       calculate variance statistics\n"
                        << "-m|--minmax         calculate minimum and maximum error probabilities for each read position of each tile\n"
                        << "-s|--sample         number of sequences to be sampled per tile for statistics estimates (default; s = 10000)\n"
                        << "-b|--bwa            use BWA trimming algorithm\n"
                        << "-d|--directory      path to directory where output files are saved\n"
                        << "--sanger            Sanger format (bypasses automatic format detection)\n"
                        << "--solexa            Solexa format (bypasses automatic format detection)\n"
                        << "--illumina          Illumina format (bypasses automatic format detection)\n"
                        << "-t|--torrent        Ion Torrent fastq file\n"
                        << std::endl;
        }
        else if(type=="dynamictrim"){
            std::cerr << "\nUsage: " << name << " dynamictrim input_files [-t|torrent] [-p|probcutoff 0.05] [-h|phredcutoff 13] [-b|bwa] [-d|directory path] [--sanger --solexa --illumina] [-t|torrent]\n"
                        << "\nOptions:\n"
                        << "-p|--probcutoff     probability value (between 0 and 1) at which base-calling error is considered too high (default; p = 0.05) *or*\n"
                        << "-h|--phredcutoff    Phred quality score (between 0 and 41) at which base-calling error is considered too high\n"
                        << "-b|--bwa            use BWA trimming algorithm\n"
                        << "-d|--directory      path to directory where output files are saved\n"
                        << "--sanger            Sanger format (bypasses automatic format detection)\n"
                        << "--solexa            Solexa format (bypasses automatic format detection)\n"
                        << "--illumina          Illumina format (bypasses automatic format detection)\n"
                        << "-a|--anchor         Reads will only be trimmed from the 3' end\n"
                        << "-t|--torrent        Ion Torrent fastq file\n"
                        << std::endl;
        }else if(type=="lengthsort"){
            std::cerr << "\nUsage: " << name << " lengthsort input_files (one single-end or two paired-end FASTQ files) [-l|length 25] [-d|directory path]\n"
                        << "\nOptions:\n"
                        << "-l|--length       length cutoff [defaults to 25 nucleotides]\n"
                        << "-d|--directory    path to directory where output files are saved\n"
                        << "-c|--correct      when running in paired mode, removes unpaired reads from the two fastq files, saves them into two new *.fastq.clean files, and normally processes them."
                        << std::endl;
        }
}


static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ( (x != n) && (x % (n/100+1) != 0) ) {
        return;}
    if (x>n) return;
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;

    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (unsigned int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

static inline void loadbar_bysize(unsigned int x, unsigned int n, unsigned int w = 50)
{
    //cout<<x<<","<<n<<"%"<<endl;
    if ( (x != n) && (x % (n/100) != 0) ) {
        return;}
    if (x>n) return;
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;

    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (unsigned int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

int check_path(const char* path){
    if (access(path, X_OK) != 0){
        return 1;
    }
    else{
        return 0;
    }
}

string check_R(char* orig_path, char* exepath){        // Need to make a Windows case: where is R?
    //cout<<std::system("Refw &")<<endl;
    int i=std::system ("Rscript --version");
    if(i==0){
        return orig_path;
    }
    else{
        cout<<"Subsidiary program R not found in your PATH. Line graphs, histogram and heatmap will not be produced."<<endl;
        return "0";
    }
//    char Rans;
//    string tmp_path;
//    const char* custom_path;
//    int custom_check;
//    namespace fs = boost::filesystem;
//    fs::path p(exepath);
//    fs::path full_path = fs::absolute(p);
//    string config_path=full_path.parent_path().string();
//    config_path=config_path+"/config.ini";
//    std::cout << config_path << std::endl;
//    const char* config_path_chr=config_path.c_str();
//    if (check_path(orig_path)==1 and fs::exists(config_path_chr)==0) {
//        cout<<"Subsidiary program R not found in /usr/bin/. Do you want to specify a custom R path? [y/n]"<<endl;
//        cin>>Rans;
//        if (Rans=='y'){
//            custom_check=1;
//            while(custom_check==1){
//                cout<<"Please insert a valid path to R executable, or quit [q]: ";
//                cin>>tmp_path;
//                if (tmp_path=="q"){
//                    cout<<"R path not specified. Line graphs, histogram and heatmap will not be produced."<<endl;
//                    return "0";
//                }
//                custom_path = tmp_path.c_str();
//                size_t len = strlen(custom_path);
//
//                if (custom_path[len-1] == '/' ||  custom_path[len-1] == '\\'){
//                    tmp_path = tmp_path.substr(0, tmp_path.size()-1);
//                }
//                tmp_path += "/R";
//                custom_path = tmp_path.c_str();
//
//                custom_check=check_path(custom_path);
//                if(custom_check==1){
//                    cout<<"R not found in "<<custom_path<<endl;
//                }
//            }
//            cout<<"R found at "<<custom_path<<endl;
//            ofstream CONFIG (config_path_chr, std::ofstream::out);
//            CONFIG<<custom_path;
//            CONFIG.close();
//
//            return custom_path;
//        }else{
//            cout<<"R path not specified. Line graphs, histogram and heatmap will not be produced."<<endl;
//            return "0";
//        }
//    }else if (check_path(orig_path)==1 and fs::exists(config_path_chr)==1){
//        cout<<"Reading R path from config file (config.ini)"<<endl;
//        string line;
//        ifstream CONFIG(config_path_chr);
//        if(CONFIG.is_open()){
//            while(getline(CONFIG,line)){
//
//                  custom_path=line.c_str();
//            }
//            CONFIG.close();
//        }
//        custom_check=check_path(custom_path);
//        if(custom_check==1){
//            cout<<"R not found in "<<custom_path<<endl;
//            remove(config_path_chr);
//            return check_R(orig_path, exepath);
//
//        }else{
//            return custom_path;
//        }
//
//    }else{
//        cout<<"R found at "<<orig_path<<endl;
//        return orig_path;
//    }

}


// Change ASCII character to Phred/Solexa quality score

map <char, int> fill_q_to_Q(string format){

    map <char, int> q_to_Q;
    int num;
    if (format=="sanger"){
        for(num=28; num<=126; ++num){
            q_to_Q[num]=num-33;
            }
    }else{
        for(num=56; num<=126; ++num){
            q_to_Q[num]=num-64;
        }
    }
//    typedef map<char, int>::const_iterator MapIterator;
//    for (MapIterator iter = q_to_Q.begin(); iter != q_to_Q.end(); iter++)
//    {
//        cout << "Key: " << iter->first << " Value:"<< iter->second <<endl;
//    }
    return q_to_Q;
}

// Change Phred/Solexa quality score to ASCII character

map <int, char> fill_Q_to_q(string format){

    map <int, char> Q_to_q;
    int num;
    if (format=="sanger"){
        for(num=-5; num<=93; ++num){
            Q_to_q[num]=num+33;
            }
    }else{
        for(num=-5; num<=93; ++num){
            Q_to_q[num]=num+64;
        }
    }
//    typedef map<int, char>::const_iterator MapIterator;
//    for (MapIterator iter = Q_to_q.begin(); iter != Q_to_q.end(); iter++)
//    {
//        cout << "Key: " << iter->first << " Value:"<< iter->second <<endl;
//    }
    return Q_to_q;
}

map <float, double> fill_Q_to_p(string format){

    map <float, double> Q_to_p;
    float num;
    if (format=="solexa"){
        for(num=-5; num<=93; ++num){
            Q_to_p[num]=pow(10,-num/10) / (pow(10, -num/10)+1);
            }
    }else{
        for(num=-5; num<=93; ++num){
            Q_to_p[num]=pow(10,-num/10);
        }
    }
//    typedef map<float, double>::const_iterator MapIterator;
//    for (MapIterator iter = Q_to_p.begin(); iter != Q_to_p.end(); iter++)
//    {
//        cout << "Key: " << iter->first << " Value:"<< iter->second <<endl;
//    }
    return Q_to_p;
}

float p_to_Q(string format, double p){

    if (format=="solexa"){
        return -10 * log10(p/(1-p));
    }else{
        return -10 * log10(p);
    }
}


long fsize(FILE* binaryStream){
  long ofs, ofs2;
  int result;

  if (fseek(binaryStream, 0, SEEK_SET) != 0 ||
      fgetc(binaryStream) == EOF)
    return 0;

  ofs = 1;

  while ((result = fseek(binaryStream, ofs, SEEK_SET)) == 0 &&
         (result = (fgetc(binaryStream) == EOF)) == 0 &&
         ofs <= LONG_MAX / 4 + 1)
    ofs *= 2;

  /* If the last seek failed, back up to the last successfully seekable offset */
  if (result != 0)
    ofs /= 2;

  for (ofs2 = ofs / 2; ofs2 != 0; ofs2 /= 2)
    if (fseek(binaryStream, ofs + ofs2, SEEK_SET) == 0 &&
        fgetc(binaryStream) != EOF)
      ofs += ofs2;

  /* Return -1 for files longer than LONG_MAX */
  if (ofs == LONG_MAX)
    return -1;

  return ofs + 1;
}

/* File must be open with 'b' in the mode parameter to fopen() */
/* Set file position to size of file before reading last line of file */
char* fgetsr(char* buf, int n, FILE* binaryStream){
  long fpos;
  int cpos;
  int first = 1;

  if (n <= 1 || (fpos = ftell(binaryStream)) == -1 || fpos == 0)
    return NULL;

  cpos = n - 1;
  buf[cpos] = '\0';

  for (;;)
  {
    int c;

    if (fseek(binaryStream, --fpos, SEEK_SET) != 0 ||
        (c = fgetc(binaryStream)) == EOF)
      return NULL;

    if (c == '\n' && first == 0) /* accept at most one '\n' */
      break;
    first = 0;

    if (c != '\r') /* ignore DOS/Windows '\r' */
    {
      unsigned char ch = c;
      if (cpos == 0)
      {
        memmove(buf + 1, buf, n - 2);
        ++cpos;
      }
      memcpy(buf + --cpos, &ch, 1);
    }

    if (fpos == 0)
    {
      fseek(binaryStream, 0, SEEK_SET);
      break;
    }
  }

  memmove(buf, buf + cpos, n - cpos);

  return buf;
}


std::tuple<int,int,int> check_header_format(const std::string& header){
    vector<string> splitspace;
    vector<string> splitcolon;
    regex sra("^@[ES]R[RA][\\d]+.[\\d]+\\s*");
        // istringstream(std::string(match[1].first, match[1].second))>>number_of_tiles;
    boost::split(splitspace, header, boost::is_any_of(" "));
//    cout<<"analyzing "<<header<<endl;
//    cout<<splitspace[0]<<endl;
    int number_of_tiles=0;

    if (boost::regex_search(splitspace[0],sra) && splitspace.size()>1) {
        boost::split(splitcolon, splitspace[1], boost::is_any_of(":"));

        if(splitcolon.size()==1){
            cout<<"\nERROR: The header of this SRA-generated fastq file does not contain tile and coordinates information, colon separated.\nTry to re run SolexaQA++ with the --torrent option."<<endl;
            return std::make_tuple(0,0,0);
        }
//        cout<<"SRA "<<*(splitcolon.rbegin() + 2)<<endl;
        if(splitcolon.size()<8){
            istringstream buffer(*(splitcolon.rbegin() + 2));
            buffer>>number_of_tiles;
            return std::make_tuple(1,2,number_of_tiles);
        }else{
            istringstream buffer(*(splitcolon.rbegin() + 3));
            buffer>>number_of_tiles;
            return std::make_tuple(1,3,number_of_tiles);

        }
    }else{
        boost::split(splitcolon, splitspace[0], boost::is_any_of(":"));
        if(splitcolon.size()<8){
            istringstream buffer(*(splitcolon.rbegin() + 2));
            buffer>>number_of_tiles;
            return std::make_tuple(0,2,number_of_tiles);

        }else{

            istringstream buffer(*(splitcolon.rbegin() + 3));
            buffer>>number_of_tiles;
            return std::make_tuple(0,3,number_of_tiles);

        }
//            if len(line.split()[0].split(':'))<8:
//               print line.split()[0].split(':')[-3]
//            else:
//               print line.split()[0].split(':')[-4]
    }
    return std::make_tuple(0,0,0);

}

unsigned int FileRead( istream & is, vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

unsigned int CountLines( const vector <char> & buff, int sz ) {
    int newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
    	if ( p[i] == '\n' ) {
    		newlines++;
    	}
    }
    return newlines;
}


bool is_file_exist(const char *fileName){
    std::ifstream infile(fileName);
    return infile.good();
}


// trim sequences using the BWA algorithm
int bwa_trim(double threshold, vector<int> array_ref){

    int length=array_ref.size();
    // only calculate if quality fails near end
    if (array_ref.back()>= threshold) {
        return length;
    }
    vector<double> arg(length, 0);
    for(int i=0; i<length-1; ++i){
        int x=i+1;
        for (int j=x; j< length; ++j) arg[x]+=threshold-array_ref[j];
    }

    // find number of 5' bases to retain
    int index=0;
    float maxval=0;
    for(int i=1; i<length; ++i) {
        if(maxval<arg[i]){
            index=i;
            maxval=arg[i];
        }
    }
    return index;
}

// trim from 3' only
int right_trim(double threshold, vector<int> array_ref){

    int length=array_ref.size();

    // find number of 5' bases to retain
    int index=0;
    for(int i=0; i<length; ++i) {
        //cout<<array_ref[i]<<endl;
        if(array_ref[i]<threshold){
            index=i;
            return index;
        }
    }
    return length;
}


// Extracts read name without the pair indication (e.g. "/1" or "1:x:x:")
string extract_read_name(string temp_header){
    vector<string> splitspace;
    vector<string> splitcolon;
    int sra_check=0;
    regex sra("^@[ES]R[RA][\\d]+.[\\d]+\\s*");
    boost::split(splitspace, temp_header, boost::is_any_of(" "));
    if (boost::regex_search(splitspace[0],sra)) sra_check=1;
    boost::split(splitcolon, splitspace[0+sra_check], boost::is_any_of("/"));
    istringstream buffer(*(splitcolon.begin()));
    string out_header;
    buffer>>out_header;
    return out_header;
}




int correct_reads(char* directory, string first_filename, string second_filename, string filepath1, string filepath2, const char* first_filepath_chr, const char* second_filepath_chr){
    gzFile fp1 = gzopen(first_filepath_chr, "r");
    gzFile fp2 = gzopen(second_filepath_chr, "r");
    cout<<"Checking pairs, please wait..."<<endl;
    int l1, l2;
    set<string> headers_r1, headers_r2, result;
    std::vector<string>::iterator it;
    std::string header1, header2, readname1, readname2;
    kseq_t *seq1 = kseq_init(fp1);
    kseq_t *seq2 = kseq_init(fp2);

    // Open output files for writing matching pairs
    const char* output_clean_1;
    const char* output_clean_2;

    string output_name_1, output_name_2;

    if(directory){
                output_name_1=directory+first_filename+".clean";
                output_clean_1=output_name_1.c_str();
                output_name_2=directory+second_filename+".clean";
                output_clean_2=output_name_2.c_str();
    }else{
                output_name_1=filepath1+".clean";
                output_clean_1=output_name_1.c_str();
                output_name_2=filepath2+".clean";
                output_clean_2=output_name_2.c_str();
    }

    if (is_file_exist(output_clean_1)){
        cout<<"Error: Clean output file "<<output_clean_1<<" already exists."<<endl;
        return -1;
    }

    if (is_file_exist(output_clean_2)){
        cout<<"Error: Clean output file "<<output_clean_2<<" already exists."<<endl;
        return -1;
    }

    ofstream CLEAN1 (output_clean_1, std::ofstream::out);
    ofstream CLEAN2 (output_clean_2, std::ofstream::out);



    // Find unmatched reads
    while ((l1 = kseq_read(seq1)) >= 0){
        if (seq1->comment.s)header1=std::string("@")+seq1->name.s+" "+seq1->comment.s;
        else header1=std::string("@")+seq1->name.s;
        readname1=extract_read_name(header1);
        headers_r1.insert(readname1);
    }

    while ((l2 = kseq_read(seq2)) >= 0){
        if (seq2->comment.s) header2=std::string("@")+seq2->name.s+" "+seq2->comment.s;
        else header2=std::string("@")+seq2->name.s;
        readname2=extract_read_name(header2);
        headers_r2.insert(readname2);

    }

    set_symmetric_difference (headers_r1.begin(), headers_r1.end(), headers_r2.begin(), headers_r2.end(), inserter(result, result.begin()));
    cout << endl << "Unpaired reads:" << endl << "-------------" << endl;
    for (set<string>::const_iterator i = result.begin(); i != result.end(); ++i) {
        cout << *i << endl;
    }

    // Clean unmatched reads and write new *_clean.fastq files
    gzclose(fp1);
    kseq_destroy(seq1);
    fp1 = gzopen(first_filepath_chr, "r");
    seq1 = kseq_init(fp1);



//    while ((l1 = kseq_read(seq1)) >= 0){
//        cout<<seq1->name.s<<seq1->comment.s<<endl;
//        if (seq1->comment.s){
//                header1=std::string("@")+seq1->name.s+" "+seq1->comment.s;
//                cout<<"Comment present"<<endl;
//        }else{
//            header1=std::string("@")+seq1->name.s;
//            cout<<"No comment here"<<endl;
//        }
//        readname1=extract_read_name(header1);
//        //cout<<"Checking"<<readname1<<endl;
//        if (result.find(readname1) != result.end()) cout<<"Cleaned from read 1: "<<readname1<<endl;
//        else{
//            cout<<"Writing"<<header1<<endl;
//            CLEAN1<<header1<<endl;
//            CLEAN1<<seq1->seq.s<<endl;
//            CLEAN1<<"+"<<endl;
//            CLEAN1<<seq1->qual.s<<endl;
//        }
//    }

    while ((l1 = kseq_read(seq1)) >= 0){
        if (seq1->comment.s) header1=std::string("@")+seq1->name.s+" "+seq1->comment.s;
        else header1=std::string("@")+seq1->name.s;
        readname1=extract_read_name(header1);
        if (result.find(readname1) != result.end()) cout<<"Cleaned from read 1: "<<readname1<<endl;
        else{
            CLEAN1<<header1<<endl;
            CLEAN1<<seq1->seq.s<<endl;
            CLEAN1<<"+"<<endl;
            CLEAN1<<seq1->qual.s<<endl;
        }
    }


    gzclose(fp2);
    kseq_destroy(seq2);
    fp2 = gzopen(second_filepath_chr, "r");
    seq2 = kseq_init(fp2);;

    while ((l2 = kseq_read(seq2)) >= 0){
        if (seq2->comment.s) header2=std::string("@")+seq2->name.s+" "+seq2->comment.s;
        else header2=std::string("@")+seq2->name.s;
        readname2=extract_read_name(header2);
        if (result.find(readname2) != result.end()) cout<<"Cleaned from read 2: "<<readname2<<endl;
        else{
            CLEAN2<<header2<<endl;
            CLEAN2<<seq2->seq.s<<endl;
            CLEAN2<<"+"<<endl;
            CLEAN2<<seq2->qual.s<<endl;
        }
    }
    cout<<"Paired reads were written to:"<<endl;
    cout<<output_name_1<<endl;
    cout<<output_name_2<<endl;
    cout<<endl;
    CLEAN1.close();
    CLEAN2.close();
    return 0;
}

string readzipped(gzFile zippedfile){
    cout<<"Reading compressed file, please wait..."<<endl;
    kseq_t *seq;
    int l;
    string temp_header;
    seq = kseq_init(zippedfile);

    if (! zippedfile) {
            cout<<"Error: could not find file"<<endl;
            exit (EXIT_FAILURE);
    }

    while ((l = kseq_read(seq)) >= 0) {
            if (seq->comment.s) temp_header=std::string("@")+seq->name.s+" "+seq->comment.s;
            else temp_header=std::string("@")+seq->name.s;
    }
    return temp_header;
}

int zipped_filesize(gzFile zippedfile){
    cout<<"Reading compressed file, please wait..."<<endl;
    #define LENGTH 0x1000
    int tot_bytes=0;
    while (1) {

        int err;
        int bytes_read;
        unsigned char buffer[LENGTH];
        bytes_read = gzread (zippedfile, buffer, LENGTH - 1);

        buffer[bytes_read] = '\0';
        tot_bytes+=bytes_read;
        if (bytes_read < LENGTH - 1) {
            if (gzeof (zippedfile)) {
                break;
            }
            else {
                const char * error_string;
                error_string = gzerror (zippedfile, & err);
                if (err) {
                    fprintf (stderr, "Error: %s.\n", error_string);
                    exit (EXIT_FAILURE);
                }
            }
        }
    }
    return tot_bytes;
}


int main(int argc, char* argv[]){

    cout<<"\nSolexaQA++ v3.1\n"
        "Released under GNU General Public License version 3\n"
        "C++ version developed by Mauro Truglio (M.Truglio@massey.ac.nz)"<<endl;

    //##### Get input from command line ##########################

    if (argc < 2 || (argv[1]!=string("analysis") && argv[1]!=string("dynamictrim") && argv[1]!=string("lengthsort")) ) {
            // Tell the user how to run the program
            show_usage(argv[0], "general");
            /* "Usage messages" are a conventional way of telling the user
             * how to run a program if they enter the command incorrectly.
             */
            return 1;
        }



    if(argv[1]==string("analysis")){
        int sample = 10000;
        double prob_cutoff=-1, phred_cutoff=-1;
        bool torrent=false, sanger = false, illumina = false, solexa = false, min_max= false, variance=false;
        bool user_defined=false;
        bool bwa=false;
        string format;
        char* directory=NULL; //Rpath;
        bool is_hiseq_32=false, is_hiseq_48=false;
        static map<int, int> hiseq_tile_32 = map_list_of  (1 , 1 )(2 , 2 )(3 , 3 )(4 , 4 )(5 , 5 )(6 , 6 )(7 , 7 )(8 , 8 )(21 , 9 )(22 , 10 )(23 , 11 )(24 , 12 )(25 , 13 )(26 , 14 ) (27 , 15 ) (28 , 16 ) (41 , 17 ) (42 , 18 ) (43 , 19 ) (44 , 20 ) (45 , 21 ) (46 , 22 ) (47 , 23 ) (48 , 24 ) (61 , 25 ) (62 , 26 ) (63 , 27 ) (64 , 28 ) (65 , 29 ) (66 , 30 ) (67 , 31 ) (68 , 32 ) (1101 , 1 ) (1102 , 2 ) (1103 , 3 ) (1104 , 4 ) (1105 , 5 ) (1106 , 6 ) (1107 , 7 ) (1108 , 8 ) (1201 , 9 ) (1202 , 10 ) (1203 , 11 ) (1204 , 12 ) (1205 , 13 ) (1206 , 14 ) (1207 , 15 ) (1208 , 16 ) (2101 , 17 ) (2102 , 18 ) (2103 , 19 ) (2104 , 20 ) (2105 , 21 ) (2106 , 22 ) (2107 , 23 ) (2108 , 24 ) (2201 , 25 ) (2202 , 26 ) (2203 , 27 ) (2204 , 28 ) (2205 , 29 ) (2206 , 30 ) (2207 , 31 ) (2208 , 32);
        static map<int, int> hiseq_tile_48 = map_list_of  (1101 , 1)(1102 , 2)(1103 , 3)(1104 , 4)(1105 , 5)(1106 , 6)(1107 , 7)(1108 , 8)(1201 , 9)(1202 , 10)(1203 , 11)(1204 , 12)(1205 , 13)(1206 , 14)(1207 , 15)(1208 , 16)(1301 , 17)(1302 , 18)(1303 , 19)(1304 , 20)(1305 , 21)(1306 , 22)(1307 , 23)(1308 , 24)(2101 , 25)(2102 , 26)(2103 , 27)(2104 , 28)(2105 , 29)(2106 , 30)(2107 , 31)(2108 , 32)(2201 , 33)(2202 , 34)(2203 , 35)(2204 , 36)(2205 , 37)(2206 , 38)(2207 , 39)(2208 , 40)(2301 , 41)(2302 , 42)(2303 , 43)(2304 , 44)(2305 , 45)(2306 , 46)(2307 , 47)(2308 , 48);
        vector<int> hiseq_keys_32;
        vector<int> hiseq_keys_48;



        for(map<int,int>::iterator it = hiseq_tile_32.begin(); it != hiseq_tile_32.end(); ++it) {
          hiseq_keys_32.push_back(it->first);
            }
        for(map<int,int>::iterator it = hiseq_tile_48.begin(); it != hiseq_tile_48.end(); ++it) {
          hiseq_keys_48.push_back(it->first);
            }

        while (1) {
            int option_index = 0;
            static struct option long_options[] = {
                {"torrent",     no_argument,        0,  't' },
                {"probcutoff",  required_argument,  0,  'p' },
                {"phredcutoff", required_argument,  0,  'h' },
                {"variance",    no_argument,        0,  'v' },
                {"minmax",      no_argument,        0,  'm' },
                {"sample",      required_argument,  0,  's' },
                {"bwa",         no_argument,        0,  'b' },
                {"directory",   required_argument,  0,  'd' },
                {"sanger",      no_argument,        0,  '1' },
                {"solexa",      no_argument,        0,  '2' },
                {"illumina",    no_argument,        0,  '3' },
                {0,             0,                  0,  0 }
            };

            int c;
            c = getopt_long(argc, argv, "p:h:tvms:bd:123",
                     long_options, &option_index);
            if (c == -1)
                break;

            switch (c) {
                   case 't':
                        torrent=true;
                        break;

                   case 'p':
                        prob_cutoff=atof(optarg);
                        cout<<"Running with p="<<prob_cutoff<<endl;
                        break;

                   case 'h':
                        phred_cutoff=atof(optarg);
                        cout<<"Running with h="<<phred_cutoff<<endl;
                        break;

                   case 'v':
                        variance=true;
                        break;

                   case 'm':
                        min_max=true;
                        break;

                   case 's':
                        sample=atoi(optarg);
                        break;

                   case 'b':
                        bwa=true;
                        break;

                   case 'd':
                        directory=optarg;
                        break;

                   case '1':
                        sanger=true;
                        break;

                   case '2':
                        solexa =true;

                        break;

                   case '3':
                        illumina=true;
                        break;

                   default:
                        show_usage(argv[0], "analysis");
                    }
                }





        //##### Get unflagged arguments (files) ###################
        vector<char*> files;
        if (optind < argc) {

            while (optind < argc)
                files.push_back(argv[optind++]);
            }
        files.erase (files.begin());
        if(files.size()==0){
            show_usage(argv[0], "analysis");
            return 1;
        }




        //##### Check if user has defined fastq format ###########

        if((sanger && solexa) || (sanger && illumina) || (solexa && illumina) || (sanger && torrent) || (solexa && torrent) || (illumina && torrent)){
            std::cerr<<"Error: Please select only one of --sanger, --solexa, --illumina or --torrent"<<endl;
            return 1;
        }
        if( sanger || solexa || illumina ){
            user_defined = true;
        }



        //##### Check parameters range ################################

        if( sample < 10000 ){
            std::cerr<<"Warning: Sample sizes less than 10,000 may lead to inaccurate estimates"<<endl;
        }

        if( ( sample > 10000 ) && variance ){
            std::cerr<<"Warning: Running variance method with sample sizes greater than 10,000 may take a long time"<<endl;
        }

        if(phred_cutoff!=-1 && prob_cutoff!=-1){
            std::cerr<<"Error: Please select p OR h cutoff, not both."<<endl;
            return 1;
        }else if (phred_cutoff==-1 && prob_cutoff==-1){
            cout<<"Running with default p=0.05"<<endl;
            prob_cutoff=0.05;
        }else if(phred_cutoff<0 and phred_cutoff!=-1){
            std::cerr<<"Error: Q cutoff must be greater than or equal to 0"<<endl;
            return 1;
        }else if((prob_cutoff!=-1 and prob_cutoff<0) || prob_cutoff>1){
            std::cerr<<"Error: P cutoff must be between 0 and 1"<<endl;
            return 1;
        }

        //##### Check for presence of R ####################################

        string R_path=check_R((char*)"/usr/bin/R", argv[0]);

        //##### Iterate through files ######################################

        for(vector<char*>::iterator it = files.begin(); it != files.end(); ++it) {
            const char* filepath_chr;
            std::string filename;
            std::string filepath;
            bool zipped=false;
            filepath = boost::filesystem::absolute(boost::filesystem::path(*it)).string();
            filename = boost::filesystem::path(*it).filename().string();
            if(filename.substr(filename.find_last_of(".") + 1) == "gz"){
                zipped=true;
            }
            filepath_chr=filepath.c_str();

            // If Windows, zipped files not supported.
            cout<<endl<<"Working on "<<filepath<<endl;
            #ifdef _WIN32
                if(zipped==true){
                    cout<<"ERROR: Compressed files are not supported in Windows. Please decompress your fastq files and run again."<<endl;
                    continue;
                }
            #endif // OS_WINDOWS

            gzFile fp;
            kseq_t *seq;
            int l;
            fp = gzopen(filepath_chr, "r");
            seq = kseq_init(fp);
            int seq_counter=0;


            // ###### Get the last header in order to read the max tile number
            int last_tile=0;
            std::tuple<int,int,int> header_format;
            string format;
            bool sanger_counter=false, solexa_counter=false, solill_counter=false;
            regex sanger_regexp("[!\"#$%&'()*+,-./0123456789:]");
            regex solexa_regexp("[;<=>?]");
            regex solill_regexp("[JKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh]");
            regex all_regexp("[@ABCDEFGHI]");

            int first_tile;
            int number_of_tiles;
            int seq_per_tile=0;
            bool warn_sample_size=false;
            bool four_digits = false;

            if(torrent==false){
                string lastheader;
                if(zipped==false){
                    FILE* f;
                    long sz;
                    if ((f = fopen(filepath_chr, "rb")) == NULL)
                      {
                        printf("Failed to open file \'%s\'\n", filepath_chr);
                        return -1;
                      }
                    sz=fsize(f);

                    if (sz > 0)
                      {
                        int counter=0;
                        char buf[1024];
                        fseek(f, sz, SEEK_SET);
                        while (fgetsr(buf, sizeof(buf), f) != NULL){
                            ++counter;

                            lastheader=buf;
                            if (counter==4) break;
                        }

                      }

                    fclose(f);
                }else{

                    lastheader=readzipped(fp);
                    gzrewind(fp); kseq_rewind(seq); // Go to beginning of file


                }

                lastheader.erase(std::remove(lastheader.begin(), lastheader.end(), '\n'), lastheader.end()); //Remove trailing newlines

                // #### Get the header format using the function check_header. It returns a tuple containing:                                                   ######################
                // #### <position of the main header field when splitting by space, position of tile number when splitting the latter by ":", last tile number> ######################


                header_format=check_header_format(lastheader);
                last_tile=std::get<2>(header_format); //Add 1 for later
                //cout<<"Tile location in header: "<<std::get<0>(header_format)<<","<<std::get<1>(header_format)<<endl;
                number_of_tiles=last_tile;

                // Tile number correction for HiSeq //
                if( last_tile == 68 || last_tile == 2208 ){
                    is_hiseq_32 = true;
                    number_of_tiles = 32;
                }else if( last_tile == 2308 ){
                    is_hiseq_48 = true;
                    number_of_tiles = 48;
                }

                // Tile number correction when tile is in the 4 digits format //
                regex tile_long("\\d\\d\\d\\d");
                if(boost::regex_search(std::to_string(last_tile),tile_long) and !is_hiseq_32 and !is_hiseq_48){
                     //cout<<"Four digits tile"<<endl;
                     four_digits=true;
                     int tile_ord=last_tile%100;
                     int first_multiplier=last_tile/1000%10;
                     int second_multiplier=last_tile/100%10;
                     number_of_tiles=tile_ord*first_multiplier*second_multiplier;
                }

                if( last_tile == 0){
                    cout<<"Error: File "<<filename<<" does not match Solexa ID format (note: this may be caused by an incomplete final entry/empty terminal lines)"<<endl;
                    continue;
                }
                //cout<<"Last tile: "<<last_tile<<endl;
                //cout<<"Tiles number: "<<number_of_tiles<<endl;


                // #### Auto-detect score format and number of reads per tile. The loop ends when the first tile ends   ####################################
                // #### and at least 10000 sequences (or all the sequences) have been read.                             ####################################


            }else{number_of_tiles=1;}

            while ((l = kseq_read(seq)) >= 0) {

                    if(torrent==false){
                        string temp_header;

                        if (seq->comment.s) temp_header=std::string("@")+seq->name.s+" "+seq->comment.s;
                        else temp_header=std::string("@")+seq->name.s;

                        vector<string> splitspace;
                        vector<string> splitcolon;

                        boost::split(splitspace, temp_header, boost::is_any_of(" "));
                        int colonpos=std::get<0>(header_format);
                        boost::split(splitcolon, splitspace[colonpos], boost::is_any_of(":"));
                        int tilenum;
                        istringstream buffer(*(splitcolon.rbegin() + std::get<1>(header_format)));
                        buffer>>tilenum;
                        if (seq_counter==0) first_tile=tilenum; // Save first tile number

                        if (tilenum!=first_tile and seq_per_tile==0){ // If tile number changes for the first time
                            seq_per_tile=seq_counter;
                            if(seq_per_tile<sample){
                                warn_sample_size=true;
                                }
                            if(seq_counter>10000) break;

                            }
                    }
                    if (seq->qual.l and (seq_counter<=10000 or seq_per_tile==0) and user_defined==false){

                            // Detect format for each quality line
                            if (regex_search(seq->qual.s, sanger_regexp)) sanger_counter = true;
                            if (regex_search(seq->qual.s, solexa_regexp)) solexa_counter = true;
                            if (regex_search(seq->qual.s, solill_regexp)) solill_counter = true;

                        }
                    else if (seq_counter>10000 and seq_per_tile!=0){
                        break;
                    }

                    ++seq_counter;
            }


            if (seq_per_tile==0) seq_per_tile=seq_counter;

            if (seq_counter==0){
                std::cerr<<"Error: file "<<filepath_chr<<" looks empty."<<endl;
                continue;
            }


            if (user_defined==false){
                format="";
            }else if(user_defined==true){
                if (sanger) format="sanger";
                if (solexa) format="solexa";
                if (illumina) format="illumina";
            }

            // ########### Final decision on format, based on the true/false results obtained while traversing reads + maps initialization based on format ######################
            if (format==""){

                if( sanger_counter==true ){
                    format = "sanger";
                }else if( sanger_counter==false && solexa_counter==true ){
                    format = "solexa";
                }else if( sanger_counter==false && solexa_counter==false && solill_counter==true ){
                    format = "illumina";
                }else{
                        std::cerr<<"Error: File format cannot be determined"<<endl;
                        return 1;

                }
            }

            //################# Print format information #########################################################################################################################
            if( format=="sanger"){
                if(user_defined==true){
                    cout<<"User defined format: Sanger FASTQ format"<<endl;
                }else{
                    cout<<"Automatic format detection: Sanger FASTQ format"<<endl;
                }
            }else if( format=="solexa"){
                if(user_defined==true ){
                    cout<<"User defined format: Solexa FASTQ format, Illumina pipeline 1.2 or less"<<endl;
                }else{
                    cout<<"Automatic format detection: Solexa FASTQ format, Illumina pipeline 1.2 or less"<<endl;
                }
            }else if( format=="illumina" ){
                if( user_defined==true ){
                    cout<<"User defined format: Illumina FASTQ format, Illumina pipeline 1.3+"<<endl;
                }else{
                    cout<<"Automatic format detection: Illumina FASTQ format, Illumina pipeline 1.3+"<<endl;
                }
            }

            map <char, int> q_to_Q=fill_q_to_Q(format);
            map <float, double> Q_to_p=fill_Q_to_p(format);
            map <int, char> Q_to_q=fill_Q_to_q(format);


            // ################# Check for phred cutoff value #####################################################################################################################
            if (phred_cutoff!=-1){
                prob_cutoff = Q_to_p[int(phred_cutoff)];
                printf("Info: Phred value of %f converted to a probability value of %f\n", phred_cutoff, prob_cutoff);
            }

            // ################# Determine bwa trimming threshold #################################################################################################################
            double threshold=0;
            if (bwa==true){
                if (phred_cutoff!=-1){
                    threshold=phred_cutoff;
                }else{
                    threshold= p_to_Q(format, prob_cutoff);
                }
                cout<<"BWA threshold "<<threshold<<endl;
            }

            // ########## Remove trailing "/" from directory ##############
            if (directory && directory[strlen(directory)-1] == '/'){
                while(directory[strlen(directory)-1] == '/'){

                    char* correctdir= new char[strlen(directory)-1];
                    strcpy(correctdir, directory);
                    correctdir[strlen(directory)-1] = '\0';
                    directory=correctdir;

                }
            }
            if (directory){
                size_t end_index = strlen(directory);
                char* correctdir= new char[end_index+2];
                strcpy(correctdir, directory);
                correctdir[end_index] = '/';
                correctdir[end_index+1] = '\0';
                directory=correctdir;
            }
            // ######### Open output quality file  #######################
            const char* quality_file;
            string quality_filename_str;
            if(directory){
                quality_filename_str=directory+filename+".quality";
                quality_file=quality_filename_str.c_str();
            }else{
                quality_filename_str=filepath+".quality";
                quality_file=quality_filename_str.c_str();
            }

            if (is_file_exist(quality_file)){
                cout<<"Error: Quality file "<<quality_file<<" already exists."<<endl;
                continue;
            }

            ofstream QUALITY (quality_file, std::ofstream::out);
            if(not QUALITY.is_open()){
                cout<<"Error: Failure opening "<<quality_file<<" for writing."<<endl;
                continue;
            }
            QUALITY << fixed << showpoint;
            QUALITY << setprecision(5);

            // ######### Open output matrix file  #######################

            const char* matrix_file;
            string matrix_filename_str;
            if(directory){
                matrix_filename_str=directory+filename+".matrix";
                matrix_file=matrix_filename_str.c_str();
            }else{
                matrix_filename_str=filepath+".matrix";
                matrix_file=matrix_filename_str.c_str();
            }

            if (is_file_exist(matrix_file)){
                cout<<"Error: Matrix file "<<matrix_file<<" already exists."<<endl;
                continue;
            }

            ofstream MATRIX (matrix_file, std::ofstream::out);
            if(not MATRIX.is_open()){
                cout<<"Error: Failure opening "<<matrix_file<<" for writing."<<endl;
                continue;
            }
            MATRIX << fixed << showpoint;
            MATRIX << setprecision(5);

            // ######### Open output segments file  #######################

            const char* segments_file;
            string segments_filename_str;
            if(directory){
                segments_filename_str=directory+filename+".segments";
                segments_file=segments_filename_str.c_str();
            }else{
                segments_filename_str=filepath+".segments";
                segments_file=segments_filename_str.c_str();
            }

            if (is_file_exist(segments_file)){
                cout<<"Error: Segments file "<<segments_file<<" already exists."<<endl;
                continue;
            }

            ofstream SEGMENT (segments_file, std::ofstream::out);
            if(not SEGMENT.is_open()){
                cout<<"Error: Failure opening "<<segments_file<<" for writing."<<endl;
                continue;
            }

            SEGMENT << fixed << showpoint;
            SEGMENT << setprecision(5);


            map <int, vector<double>> mean_bytiles;
            map <int, vector<double>> mean_count_bytiles;
            map <int, vector<double>> min_bytiles;
            map <int, vector<double>> max_bytiles;
            map <int, vector<double>> var_sum_bytiles;
            map <int, vector<double>> var_count_bytiles;
            map <int, double> trim_histogram;

            // Each tile has its own set of arrays
            for( int i = 0; i < number_of_tiles+1; ++i ){

                vector <double> new_array1;
                vector <double> new_array2;
                vector <double> new_array3(9999, 1);                 // Filled by 1
                vector <double> new_array4(9999, 0);                 // Filled by 0
                vector <double> new_array5;
                vector <double> new_array6;


                mean_bytiles[i]=new_array1;
                mean_count_bytiles[i]=new_array2;
                min_bytiles[i]=new_array3;
                max_bytiles[i]=new_array4;
                var_sum_bytiles[i]=new_array5;
                var_count_bytiles[i]=new_array6;
            }




            gzrewind(fp); kseq_rewind(seq); // Go to beginning of file

            seq_counter = 0;
            float percentage;
            if (warn_sample_size==true) {
                percentage=1;
                cout<<"Warning: Desired sample size ("<<sample<<") greater than number of sequences per tile ("<<seq_per_tile<<")"<<endl;
                }
            else percentage=float(sample)/seq_per_tile;
            int tot_seqs=seq_per_tile*number_of_tiles*percentage;

            //cout<<"Percentage "<<percentage<<endl<<"Tot seqs "<<tot_seqs<<endl;

            // #### Loop through sequences ########################################################################

            int previous_first_multiplier=0;
            int previous_second_multiplier=0;
            int previous_tile_number=0;
            int place_holder=0;
            int max_read_length=0;

            while ((l = kseq_read(seq)) >= 0) {
                double r = ((double) rand() / (RAND_MAX));
                if (r<= percentage){
                    ++seq_counter;
                    vector<int> qual;
                    vector<double> prob;
                    int tile_number;
                    string temp_header;
                    if(torrent==false){
                        if (seq->comment.s) temp_header=std::string("@")+seq->name.s+" "+seq->comment.s;
                        else temp_header=std::string("@")+seq->name.s;
                        vector<string> splitspace;
                        vector<string> splitcolon;

                        boost::split(splitspace, temp_header, boost::is_any_of(" "));
                        int colonpos=std::get<0>(header_format);
                        boost::split(splitcolon, splitspace[colonpos], boost::is_any_of(":"));
                        istringstream buffer(*(splitcolon.rbegin() + std::get<1>(header_format)));
                        buffer>>tile_number;

                        if(is_hiseq_32){
                            tile_number=hiseq_tile_32[tile_number];
                        }
                        else if (is_hiseq_48){
                            tile_number=hiseq_tile_48[tile_number];
                        }
                        else if(!is_hiseq_32 and !is_hiseq_48 and four_digits==true){
                                int tile_n=tile_number%100;
                                int first_multiplier=tile_number/1000%10;
                                int second_multiplier=tile_number/100%10;
                                if((first_multiplier==previous_first_multiplier and second_multiplier==previous_second_multiplier) or seq_counter==1){
                                    tile_number=tile_n+place_holder;
                                }else if(seq_counter!=1){
                                    place_holder=previous_tile_number;
                                    tile_number=tile_n+place_holder;
                                }
                                previous_first_multiplier=first_multiplier;
                                previous_second_multiplier=second_multiplier;
                                previous_tile_number=tile_number;
                        }
                    }else{
                        tile_number=1;
                    }
                    char current_Q;
                    double current_p;
                    int char_counter=0;
                    for(char* it2 = seq->qual.s; *it2; ++it2) {
                        current_Q=q_to_Q[*it2];
                        current_p=Q_to_p[current_Q];
                        qual.push_back(current_Q);
                        prob.push_back(current_p);
                        if (mean_bytiles[tile_number].size()<(unsigned)char_counter+1) mean_bytiles[tile_number].push_back(0);
                        if (mean_bytiles[0].size()<(unsigned)char_counter+1) mean_bytiles[0].push_back(0);
                        if (mean_count_bytiles[tile_number].size()<(unsigned)char_counter+1) mean_count_bytiles[tile_number].push_back(0);
                        if (mean_count_bytiles[0].size()<(unsigned)char_counter+1) mean_count_bytiles[0].push_back(0);
                        if (min_bytiles[tile_number].size()<(unsigned)char_counter+1) min_bytiles[tile_number].push_back(99);
                        if (max_bytiles[tile_number].size()<(unsigned)char_counter+1) max_bytiles[tile_number].push_back(0);



                        mean_bytiles[tile_number][char_counter]+=current_p;
                        mean_count_bytiles[tile_number][char_counter]+=1;

                        // if option -m activated, evaluate each probability against current min and max, store if better

                        if (min_max and current_p<min_bytiles[tile_number][char_counter]) min_bytiles[tile_number][char_counter]=current_p;
                        if (min_max and current_p>max_bytiles[tile_number][char_counter]) max_bytiles[tile_number][char_counter]=current_p;
                        // tile 0 represents global probability
                        mean_bytiles[0][char_counter]+=current_p;
                        mean_count_bytiles[0][char_counter]+=1;
                        if (min_max and current_p<min_bytiles[0][char_counter]) min_bytiles[0][char_counter]=current_p;
                        if (min_max and current_p>max_bytiles[0][char_counter]) max_bytiles[0][char_counter]=current_p;

                        char_counter+=1;
                    }

                    int longest_segment=0, current_start=0, cutoff_hit=0;
                    if(bwa==true) longest_segment=bwa_trim(threshold, qual);
                    else{
                        for(unsigned int i=0; i<(unsigned)prob.size(); ++i){
                            if(prob[i]>=prob_cutoff){
                                cutoff_hit = 1;
                                int current_segment_length=i-current_start;
                                current_start=i+1;
                                if(current_segment_length > longest_segment){
                                    longest_segment=current_segment_length;
                                }
                            }
                        }
                        if(cutoff_hit==0) longest_segment=-1; // I still don't know the max read length, so I'm temporarily storing counts to key -1
                    }

                    trim_histogram[longest_segment]+=1;

                    if(prob.size()>(unsigned)max_read_length){
                        max_read_length=prob.size();
                    }
                    loadbar(seq_counter, tot_seqs);
                }

            }

            loadbar(seq_counter, seq_counter);
            cout<<"\nWriting results...."<<endl;

            for(int j = 0; j <= max_read_length; ++j){
                if(!trim_histogram[j]) trim_histogram[j]=0; // Fill the map trim_histogram with zeros in empty positions
            }
            trim_histogram[max_read_length]=trim_histogram[-1]; // Now I know the max_read_length, so moving counts from trim_histogram[-1] to trim_histogram[max_read_length]

            // #### End of loop. Calculate means and store in maps ###################################################################################
            for( int j = 0; j < max_read_length; j++ ){
                    for( int i = 0; i <= number_of_tiles; i++){
                            if (mean_bytiles[i].empty()==false){
                                if( mean_bytiles[i][0] ){
                                    double mean_pertile = mean_bytiles[i][j] / mean_count_bytiles[i][j];
                                    mean_bytiles[i][j] = mean_pertile;
                                    }
                            }
                    }

            }

            // #### If calculating variance, loop through file again to find differences from mean ####################################################
            if(variance==true){


                gzrewind(fp); kseq_rewind(seq); // Go to beginning of file

                for(int j=0; j<max_read_length; ++j){   // Now I know max_read_length, I can initialize var_sum_bytiles and var_count_bytiles to 0.
                    for(int i=0; i<=number_of_tiles; ++i){
                        var_sum_bytiles[i].push_back(0);
                        var_count_bytiles[i].push_back(0);
                    }
                }
                seq_counter = 0;

                previous_first_multiplier=0;
                previous_second_multiplier=0;
                previous_tile_number=0;
                place_holder=0;

                while ((l = kseq_read(seq)) >= 0) {
                    double r = ((double) rand() / (RAND_MAX));
                    if (r<= percentage){
                        ++seq_counter;
                        vector<int> qual;
                        vector<double> prob;
                        int tile_number;
                        string temp_header;

                        if (seq->comment.s) temp_header=std::string("@")+seq->name.s+" "+seq->comment.s;
                        else temp_header=std::string("@")+seq->name.s;

                        vector<string> splitspace;
                        vector<string> splitcolon;

                        boost::split(splitspace, temp_header, boost::is_any_of(" "));
                        int colonpos=std::get<0>(header_format);
                        boost::split(splitcolon, splitspace[colonpos], boost::is_any_of(":"));
                        istringstream buffer(*(splitcolon.rbegin() + std::get<1>(header_format)));
                        buffer>>tile_number;

                        if(is_hiseq_32){
                            tile_number=hiseq_tile_32[tile_number];
                        }
                        else if (is_hiseq_48){
                            tile_number=hiseq_tile_48[tile_number];
                        }
                        else if(!is_hiseq_32 and !is_hiseq_48 and four_digits==true){
                                int tile_n=tile_number%100;
                                int first_multiplier=tile_number/1000%10;
                                int second_multiplier=tile_number/100%10;
                                if((first_multiplier==previous_first_multiplier and second_multiplier==previous_second_multiplier) or seq_counter==1){
                                    tile_number=tile_n+place_holder;
                                }else if(seq_counter!=1){
                                    place_holder=previous_tile_number;
                                    tile_number=tile_n+place_holder;
                                }
                                previous_first_multiplier=first_multiplier;
                                previous_second_multiplier=second_multiplier;
                                previous_tile_number=tile_number;
                        }

                        char current_Q;
                        double current_p;
                        for(char* it2 = seq->qual.s; *it2; ++it2) {
                            current_Q=q_to_Q[*it2];
                            current_p=Q_to_p[current_Q];
                            qual.push_back(current_Q);
                            prob.push_back(current_p);
                        }

                        // Calculate sum of squared differences from the mean
                        for(unsigned int i=0; i<(unsigned)prob.size(); ++i){
                            if (mean_bytiles[i].empty()==false){
                                double diff_from_mean = mean_bytiles[tile_number][i]-prob[i];
                                var_sum_bytiles[tile_number][i]+=pow(diff_from_mean, 2);
                                var_count_bytiles[tile_number][i]+=1;

                                // Tile 0 represents global probability
                                var_sum_bytiles[0][i]+= pow(diff_from_mean, 2);
                                var_count_bytiles[0][i]+=1;
                            }
                        }


                    }
                }
            }

            // #### Print output to quality file ######################################################


            // print header line of quality file
            for( int i = 0; i <= number_of_tiles; ++i ){
                if (mean_bytiles[i].empty()==false){
                    if( mean_count_bytiles[i][0] ){
                        // assign column ID (G if global, otherwise tile #)
                        string col_id;
                        if(i==0){
                            col_id="global";
                        }else{
                            col_id="tile_"+std::to_string(i);
                        }
                        QUALITY<<"mean_"<<col_id<<"\t";
                        if(variance==true) QUALITY<<"var_"<<col_id<<"\t";
                        if(min_max==true) QUALITY<<"min_"<<col_id<<"\tmax_"<<col_id<<"\t";
                    }
                }
            }

            QUALITY<<endl;

            for(int j=0; j<max_read_length; ++j){
                for(int i=0; i<=number_of_tiles; ++i){
                    if (mean_bytiles[i].empty()==false){
                        if(mean_count_bytiles[i][0]>0){

                            // Print mean, variance, min and max
                            QUALITY<<mean_bytiles[i][j]<<"\t";

                            if(variance==true){
                                double var_pertile = var_sum_bytiles[i][j] / var_count_bytiles[i][j];
                                QUALITY<<var_pertile<<"\t";
                            }

                            if(min_max==true){
                                double Min = min_bytiles[i][j];
                                double Max = max_bytiles[i][j];

                                QUALITY<<Min<<"\t"<<Max<<"\t";
                            }


                        }
                    }
                }
                QUALITY<<endl;
            }



            // #### Print output to matrix file for heatmaps ######################################################

            // Print header line with the column headers "Tile" and each read position
            MATRIX<<".\t";
            for(int i=1; i<=max_read_length; ++i){
                MATRIX<<i<<"\t";
            }
            MATRIX<<endl;

            // Print stored mean probability values for each read position of each tile

            // If needed, translate back tiles to format XXXX (e.g from 1 to 1101), otherwise use normal number.
            int tiles_per_lane, surface=1, lane=1, current_tile=1, flag=0;
            if (!is_hiseq_32 and !is_hiseq_48 and four_digits==true){ //In case there has been a translation from XXXX to x
                tiles_per_lane=previous_tile_number/(previous_first_multiplier*previous_second_multiplier);
            }

            for(int i=1; i<=number_of_tiles; ++i){
                if (mean_bytiles[i].empty()==false){
                    bool all_zeros=true;
                    for( int j = 0; j < max_read_length; ++j ){
                        if( mean_count_bytiles[i][j] > 0 ){
                                all_zeros = false;
                        }
                    }
                    if(all_zeros==true) continue;
                    int real_tile=0;
                    // change tile names if hiseq
                    if( is_hiseq_32==true ){
                        vector<int> sorted_hiseq_keys_32=hiseq_keys_32;
                        std::sort(sorted_hiseq_keys_32.begin(), sorted_hiseq_keys_32.end());
                        std::reverse(sorted_hiseq_keys_32.begin(),sorted_hiseq_keys_32.end());
                        for (auto &entry : sorted_hiseq_keys_32){
                            if( hiseq_tile_32[entry] == i and entry < 69 ){
                                real_tile = entry;
                            }
                        }
                        MATRIX<<"tile_"<<real_tile<<"  ("<<std::setprecision(3)<<mean_count_bytiles[i][0]/mean_count_bytiles[0][0]<<")\t";

                    }else if( is_hiseq_48==true ){
                        vector<int> sorted_hiseq_keys_48=hiseq_keys_48;
                        std::sort(sorted_hiseq_keys_48.begin(), sorted_hiseq_keys_48.end());
                        std::reverse(sorted_hiseq_keys_48.begin(),sorted_hiseq_keys_48.end());
                        for (auto &entry : sorted_hiseq_keys_48){
                            if( hiseq_tile_48[entry] == i){
                                real_tile = entry;
                            }
                        }
                        MATRIX<<"tile_"<<real_tile<<"  ("<<std::setprecision(3)<<mean_count_bytiles[i][0]/mean_count_bytiles[0][0]<<")\t";

                    //change back tile names if original format was XXXX
                    }else if(!is_hiseq_32 and !is_hiseq_48 and four_digits==true){
                        if (i%tiles_per_lane!=0){
                            current_tile=i%tiles_per_lane;
                            }
                        else{
                            current_tile=tiles_per_lane;
                            flag=1;
                        }

                        MATRIX<<"tile_"<<surface<<lane<<std::setw(2)<<std::setfill('0')<<current_tile<<"  ("<<std::setprecision(3)<<mean_count_bytiles[i][0]/mean_count_bytiles[0][0]<<")\t";

                        if (flag==1){
                             lane+=1;
                             if (lane==previous_second_multiplier+1){
                                        surface+=1;
                                        lane=1;
                             }
                             flag=0;
                        }

                    }else{
                        MATRIX<<"tile_"<<i<<"  ("<<std::setprecision(3)<<mean_count_bytiles[i][0]/mean_count_bytiles[0][0]<<")\t";

                    }

                    // print means by tiles
                    for(int j = 0; j < max_read_length; ++j){

                        if( mean_count_bytiles[i][j] > 0 ){

                            MATRIX<<std::setprecision(5)<<mean_bytiles[i][j]<<"\t";

                        }else{
                            MATRIX<<"\t";
                        }
                    }

                    MATRIX<<endl;
                }
            }


            // #### Print output to segments file for histograms ######################################################

            double segment_sum=0;
            for(int j = 0; j <= max_read_length; ++j){
                segment_sum+=trim_histogram[j];
            }
            SEGMENT<<"read_length\tproportion_of_reads"<<endl;
            double proportion;
            for(int j = 0; j <= max_read_length; ++j){
                proportion = trim_histogram[j]/segment_sum;
                SEGMENT<<j<<"\t"<<proportion<<endl;
            }



    //        cout<<"Max read: "<<max_read_length<<endl;
    //        cout<<"Reads: "<<seq_counter<<endl;



            // #### Close files #########################################################################

            SEGMENT.close();
            MATRIX.close();
            QUALITY.close();
            kseq_destroy(seq);
            gzclose(fp);


            // #### R graphs ############################################################################

            if (R_path!="0"){
                cout<<"Generating graphs...."<<endl;

                // Heatmap from matrix
                R_heatmap(directory, filepath, filename, matrix_filename_str);

                // Quality graph
                int extra_columns = 0;
                if( variance==true ){
                    extra_columns += 1;
                }
                if( min_max==true ){
                    extra_columns += 2;
                }
                R_quality(directory, filepath, filename, quality_filename_str, extra_columns);

                // Histogram from segments
                R_histogram(directory, filepath, filename, segments_filename_str, prob_cutoff);

                // Cumulative from segments
                R_cumulative(directory, filepath, filename, segments_filename_str, prob_cutoff);

            }
            cout<<"Analysis complete"<<endl;

        }
    } // End of "analysis" program


    if(argv[1]==string("dynamictrim")){

        char* directory=NULL; //Rpath;
        double prob_cutoff=-1, phred_cutoff=-1;
        string format;
        char* poor_quality_char=(char *)"B";
        bool bwa=false, user_defined=false, torrent=false, sanger = false, illumina = false, solexa = false, anchor = false;

        while (1) {
            int option_index = 0;
            static struct option long_options[] = {
                {"anchor",      no_argument,        0,  'a' },
                {"torrent",     no_argument,        0,  't' },
                {"probcutoff",  required_argument,  0,  'p' },
                {"phredcutoff", required_argument,  0,  'h' },
                {"bwa",         no_argument,        0,  'b' },
                {"directory",   required_argument,  0,  'd' },
                {"sanger",      no_argument,        0,  '1' },
                {"solexa",      no_argument,        0,  '2' },
                {"illumina",    no_argument,        0,  '3' },
                {0,             0,                  0,  0 }
            };

            int c;
            c = getopt_long(argc, argv, "p:h:tvms:bad:123",
                     long_options, &option_index);
            if (c == -1)
                break;

            switch (c) {
                   case 'a':
                        anchor=true;
                        break;

                   case 't':
                        torrent=true;
                        break;

                   case 'p':
                        prob_cutoff=atof(optarg);
                        cout<<"Running with p="<<prob_cutoff<<endl;
                        break;

                   case 'h':
                        phred_cutoff=atof(optarg);
                        cout<<"Running with h="<<phred_cutoff<<endl;
                        break;

                   case 'b':
                        bwa=true;
                        break;

                   case 'd':
                        directory=optarg;
                        break;

                   case '1':
                        sanger=true;
                        break;

                   case '2':
                        solexa =true;

                        break;

                   case '3':
                        illumina=true;
                        break;

                   default:
                        show_usage(argv[0], "dynamictrim");
                    }
        }



        //##### Get unflagged arguments (files) ###################
        vector<char*> files;
        if (optind < argc) {

            while (optind < argc)
                files.push_back(argv[optind++]);
            }
        files.erase (files.begin());
        if(files.size()==0){
            show_usage(argv[0], "dynamictrim");
            return 1;
            }

        //##### Check if user has defined fastq format ###########

        if( (sanger && solexa) || (sanger && illumina) || (solexa && illumina) || (sanger && torrent) || (solexa && torrent) || (illumina && torrent)){
            std::cerr<<"Error: Please select only one of --sanger, --solexa, --illumina or --torrent"<<endl;
            return 1;
        }
        if( sanger || solexa || illumina || torrent ){
            user_defined = true;
            if (sanger==true || torrent==true) format="sanger";
            if (solexa==true) format="solexa";
            if (illumina==true) format="illumina";
        }

        if (user_defined==false){
                format="";
        }

        //##### Check parameters range ################################

        if( phred_cutoff!=-1 && prob_cutoff!=-1){
            std::cerr<<"Error: Please select p OR h cutoff, not both."<<endl;
            return 1;
        }else if (phred_cutoff==-1 && prob_cutoff==-1){
            cout<<"Info: Using default quality cutoff of p=0.05 (change with -p or -h flag)"<<endl;
            prob_cutoff=0.05;
        }else if(phred_cutoff<0 and phred_cutoff!=-1){
            std::cerr<<"Error: Q cutoff must be greater than or equal to 0"<<endl;
            return 1;
        }else if((prob_cutoff!=-1 and prob_cutoff<0) || prob_cutoff>1){
            std::cerr<<"Error: P cutoff must be between 0 and 1"<<endl;
            return 1;
        }

        //##### Check for presence of R ####################################

        string R_path=check_R((char *)"/usr/bin/R", argv[0]);


        //##### Iterate through files ######################################        #!!! Windows warning: filename() splits by '/' and not \\

        for(vector<char*>::iterator it = files.begin(); it != files.end(); ++it) {
            const char* filepath_chr;
            std::string filename;
            std::string filepath;
            int sz;
            filepath = boost::filesystem::absolute(boost::filesystem::path(*it)).string();
            filename = boost::filesystem::path(*it).filename().string();
            bool zipped=false;
            if(filename.substr(filename.find_last_of(".") + 1) == "gz"){
                zipped=true;
            }
            //cout<<"FILENAME: "<<filename<<endl;
            filepath_chr=filepath.c_str();
            cout<<endl<<"Working on "<<filepath<<endl;

            // If Windows, zipped files not supported.
            #ifdef _WIN32
                if(zipped==true){
                    cout<<"Compressed files are not supported in Windows. Please decompress your fastq files and run again."<<endl;
                    continue;
                }
            #endif // OS_WINDOWS

            gzFile fp;
            kseq_t *seq;
            int l;
            fp = gzopen(filepath_chr, "r");
            seq = kseq_init(fp);
            int seq_counter=0;


            // #### Auto-detect score format.  ####################################
            bool sanger_counter=false, solexa_counter=false, solill_counter=false;
            regex sanger_regexp("[!\"#$%&'()*+,-./0123456789:]");
            regex solexa_regexp("[;<=>?]");
            regex solill_regexp("[JKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh]");
            regex all_regexp("[@ABCDEFGHI]");
            if(zipped==false){
                FILE* f;
                if ((f = fopen(filepath_chr, "rb")) == NULL)
                  {
                    printf("Failed to open file \'%s\'\n", filepath_chr);
                    return -1;
                  }
                sz=fsize(f);
                fclose(f);
            }else{
                sz=zipped_filesize(fp);
                gzrewind(fp); kseq_rewind(seq); // Go to beginning of file

            }
            if(torrent==false && user_defined==false){
                while ((l = kseq_read(seq)) >= 0) {
                    if (seq->qual.l and seq_counter<=10000 and user_defined==false){

                            // Detect format for each quality line
                            if (regex_search(seq->qual.s, sanger_regexp)) sanger_counter = true;
                            if (regex_search(seq->qual.s, solexa_regexp)) solexa_counter = true;
                            if (regex_search(seq->qual.s, solill_regexp)) solill_counter = true;

                        }
                    else if (seq_counter>10000){
                        break;
                    }

                    ++seq_counter;

                }
                // ########### Final decision on format, based on the true/false results obtained while traversing reads + maps initialization based on format ######################
                if (format==""){

                    if( sanger_counter==true ){
                        format = "sanger";
                    }else if( sanger_counter==false && solexa_counter==true ){
                        format = "solexa";
                    }else if( sanger_counter==false && solexa_counter==false && solill_counter==true ){
                        format = "illumina";
                    }else{
                            std::cerr<<"Error: File format cannot be determined"<<endl;
                            return 1;

                    }
                }



            }
            //################# Print format information #########################################################################################################################
            if( format=="sanger"){
                poor_quality_char=(char *)"!";
                if(user_defined==true){
                    cout<<"User defined format: Sanger FASTQ format"<<endl;
                }else{
                    cout<<"Automatic format detection: Sanger FASTQ format"<<endl;
                }
            }else if(format=="solexa"){
                poor_quality_char=(char *)";";
                if(user_defined==true ){
                    cout<<"User defined format: Solexa FASTQ format, Illumina pipeline 1.2 or less"<<endl;
                }else{
                    cout<<"Automatic format detection: Solexa FASTQ format, Illumina pipeline 1.2 or less"<<endl;
                }
            }else if(format=="illumina" ){
                poor_quality_char=(char *)"@";
                if( user_defined==true ){
                    cout<<"User defined format: Illumina FASTQ format, Illumina pipeline 1.3+"<<endl;
                }else{
                    cout<<"Automatic format detection: Illumina FASTQ format, Illumina pipeline 1.3+"<<endl;
                }
            }

            map <char, int> q_to_Q=fill_q_to_Q(format);
            map <float, double> Q_to_p=fill_Q_to_p(format);
            map <int, char> Q_to_q=fill_Q_to_q(format);
            char ascii_cutoff;
            if(phred_cutoff!=-1){
                ascii_cutoff = Q_to_q[phred_cutoff];
                prob_cutoff = Q_to_p[int(phred_cutoff)];
            }else{
                phred_cutoff=p_to_Q(format, prob_cutoff);
                ascii_cutoff = Q_to_q[phred_cutoff];
            }


            // ################# Determine bwa trimming threshold ########################################################################################
            double threshold=0;
            if (bwa==true or anchor==true){
                if (phred_cutoff!=-1){
                    threshold=phred_cutoff;
                }else{
                    threshold= p_to_Q(format, prob_cutoff);
                }
                //cout<<"Quality threshold "<<threshold<<endl;
            }


            // ########## Remove trailing "/" from directory ##############
            if (directory && directory[strlen(directory)-1] == '/'){
                while(directory[strlen(directory)-1] == '/'){
                    char* correctdir= new char[strlen(directory)-1];
                    strcpy(correctdir, directory);
                    correctdir[strlen(directory)-1] = '\0';
                    directory=correctdir;
                }
            }
            if (directory){
                size_t end_index = strlen(directory);
                char* correctdir= new char[end_index+2];
                strcpy(correctdir, directory);
                correctdir[end_index] = '/';
                correctdir[end_index+1] = '\0';
                directory=correctdir;
            }


            // ######### Open output .trimmed file  #######################
            const char* output_file;
            string output_filename_str;
            if(directory){
                output_filename_str=directory+filename+".trimmed";
                output_file=output_filename_str.c_str();
            }else{
                output_filename_str=filepath+".trimmed";
                output_file=output_filename_str.c_str();
            }

            if (is_file_exist(output_file)){
                cout<<"Error: Output file "<<output_file<<" already exists."<<endl;
                continue;
            }

            ofstream OUTPUT (output_file, std::ofstream::out);
            if(not OUTPUT.is_open()){
                cout<<"Error: Failure opening "<<output_file<<" for writing."<<endl;
                continue;
            }


            // ######### Open output segments file  #######################

            const char* segments_file;
            string segments_filename_str;
            if(directory){
                segments_filename_str=directory+filename+"_trimmed.segments";
                segments_file=segments_filename_str.c_str();
            }else{
                segments_filename_str=filepath+"_trimmed.segments";
                segments_file=segments_filename_str.c_str();
            }

            if (is_file_exist(segments_file)){
                cout<<"Error: Segments file "<<segments_file<<" already exists."<<endl;
                continue;
            }

            ofstream SEGMENT (segments_file, std::ofstream::out);
            if(not SEGMENT.is_open()){
                cout<<"Error: Failure opening "<<segments_file<<" for writing."<<endl;
                continue;
            }

            SEGMENT << fixed << showpoint;
            SEGMENT << setprecision(5);


            // ######### Step through input and trim #############################################################
            gzrewind(fp); kseq_rewind(seq); // Go to beginning of file
            map<int, int> segment_hist;
            map <int, int> Hash;
            int segment_count=0, original_length, max_length=0;
            double segment_sum = 0;
            seq_counter=0;
            std::string header;
            int progress_unit=int(sz/100);
            int percent_progress=int(sz/100);
            int progress_counter=0;


            while ((l = kseq_read(seq)) >= 0) {
                progress_counter+=strlen(seq->qual.s)+strlen(seq->seq.s)+strlen(seq->name.s);
                if(seq->comment.s!=NULL){
                    progress_counter+=strlen(seq->comment.s);
                    }
                if (progress_counter>percent_progress){
                    loadbar_bysize(percent_progress, sz);
                    percent_progress+=progress_unit;
                }
                vector<int> qual (strlen(seq->qual.s), -1);
                if (seq->comment.s) header=std::string("@")+seq->name.s+" "+seq->comment.s;
                else header=std::string("@")+seq->name.s;
                original_length=strlen(seq->qual.s);
                if(original_length>max_length) max_length=original_length;
                int best_start_index=0, best_length=0, current_start=0, bad_first=0;
                bool cutoff_hit=false;

                // BWA //
                if(bwa==true){
                    for(unsigned int i=0; i<(unsigned)strlen(seq->qual.s); ++i){
                        qual[i] = q_to_Q[seq->qual.s[i]];
                    }
                    if(qual[0]==-1){
                        bad_first=1;
                        best_length=0;
                    }else if(qual[0]<threshold){
                        bad_first=1;
                        best_length=bwa_trim(threshold, qual);
                    }else{
                        best_length= bwa_trim(threshold, qual);

                    }
                }else if(anchor==true){
                    // Trim from 3' only
                    for(unsigned int i=0; i<(unsigned)strlen(seq->qual.s); ++i){
                        qual[i] = q_to_Q[seq->qual.s[i]];
                    }

                    if(qual[0]<threshold){
                        bad_first=1;
                        best_length=0;
                    }else{
                        best_length=right_trim(threshold, qual);
                    }
                }else{
                    // Normal trimming
                    int current_segment_length;
                    // Loop through each position in quality string
                    for(int i=0; i<original_length; ++i){
                        if(seq->qual.s[i]<= ascii_cutoff){ // If quality lower than cutoff
                            cutoff_hit=true;
                            // Determine length of good segment that just ended
                            current_segment_length=i-current_start;
                            // If this segment is the longest so far
                            if (current_segment_length > best_length){
                                best_length=current_segment_length;
                                best_start_index=current_start;
                            }
                            // Reset start
                            current_start=i+1;
                        }else if(i==original_length-1){
                            current_segment_length=i+1-current_start;
                            if (current_segment_length > best_length){
                                best_length=current_segment_length;
                                best_start_index=current_start;
                            }

                        }

                    }

                    // if quality cutoff is never exceeded, set the marker for the end of the good segment
                    // to the end of the read
                    if( cutoff_hit==false ){
                        best_length = original_length;
                    }

                }

                // Update the segment_hist map, and all the other variables that hold segment stats
                if(!segment_hist[best_length]) segment_hist[best_length]=0;

                segment_hist[best_length]+=1;
                segment_sum+=best_length;
                segment_count+=1;
                if(Hash[best_length]){
                    Hash[best_length]+=1;
                }else{
                    Hash[best_length]=1;
                }
                // Actual trimming based on best segment found before
                string outseq, outqual;
                string tempseq = string(seq->seq.s);
                string tempqual = string(seq->qual.s);

                if(bwa==true){
                    if(best_length<=1 and bad_first==1){
                        outseq="N";
                        outqual=poor_quality_char;
                    }else{
                        outseq=tempseq.substr(0, best_length);
                        outqual=tempqual.substr(0, best_length);
                    }
                }else{
                    if(best_length <= 0){
                        outseq = "N";
                        outqual = poor_quality_char;
                    }else{
                        outseq = tempseq.substr(best_start_index, best_length);
                        outqual = tempqual.substr(best_start_index, best_length);
                    }

                }
                OUTPUT<<header<<endl<<outseq<<endl<<"+"<<endl<<outqual<<endl;
            }
            loadbar(100, 100);
            cout<<"\nWriting files...."<<endl;

            // Calculate mean segment length
            double segment_mean = segment_sum/segment_count;

            // Set index at halfway through segment counts
            double halfway_index = double(segment_count) / 2;
            // Set variables needed to find median
            int current_sum   = 0;
            int current_index = 0;
            int median_index1 = -1;
            int median_index2 = -1;

            while(median_index1==-1 or median_index2==-1){ // while median indexes are not defined
                // add segment count to current sum for each segment length from array
                if(segment_hist[current_index]) current_sum+= segment_hist[current_index];

                //if current sum of segment counts has surpassed halfway index
                if(current_sum>halfway_index){

                    if(median_index1==-1) median_index1=current_index;
                    if(median_index2==-1) median_index2=current_index;

                //else if current sum of segment counts is exactly equal to the halfway index
                }else if(current_sum==halfway_index and median_index1==-1){
                    median_index1=current_index;
                }
                ++current_index;
            }
            --current_index;

            double segment_median;

            // if number of segments is odd, store index2 as median segment length
            if( segment_count % 2 == 1) segment_median = median_index1;

            // if number of segments is even, store average of index1 and index2 as median segment length
            else segment_median=(median_index1+median_index2)/2;
            cout<<"\nInfo: "<<filename<<": mean segment length = "<<std::fixed<<std::setprecision(1)<<segment_mean<<", median segment length = "<< std::setprecision(0)<<segment_median<<endl;



            //########## Segments file generator ##########################################################################################################################################################

            SEGMENT<<"Read_length\tProportion_of_reads\n";
            double percentage;
            for(int i=0; i<=max_length; ++i){
                if(Hash[i]){
                    percentage=double(Hash[i])/segment_count;
                    SEGMENT<<i<<"\t"<<percentage<<endl;
                }else{
                    SEGMENT<<i<<"\t0"<<endl;
                }
            }

            // ######### Close files ######################################################################################################################################################################
            OUTPUT.close();
            SEGMENT.close();
            gzclose(fp);


            // #### R graphs ############################################################################
            if (R_path!="0"){
                // Histogram from segments
                R_histogram(directory, filepath, filename, segments_filename_str, prob_cutoff);

            }

        } //end files loop

    }

    if(argv[1]==string("lengthsort")){
        int length=25;
        bool paired = false, correct=false;
        char* directory=NULL;

        while (1) {
            int option_index = 0;
            static struct option long_options[] = {
                {"length",     required_argument,  0,  'l' },
                {"directory",  required_argument,  0,  'd' },
                {"correct",    no_argument,        0,  'c' },
                {0,             0,                 0,  0 }
            };

            int c;
            c = getopt_long(argc, argv, "cl:d:",long_options, &option_index);
            if (c == -1)
                break;

            switch (c) {
                   case 'l':
                        length=atoi(optarg);
                        break;

                   case 'd':
                        directory=optarg;
                        break;

                   case 'c':
                        correct=true;
                        break;

                   default:
                        show_usage(argv[0], "lengthsort");
                    }
        }


        //##### Get unflagged arguments (files) ###################
        vector<char*> files;
        if (optind < argc) {

            while (optind < argc)
                files.push_back(argv[optind++]);
            }
        files.erase (files.begin());
        if(files.size()==0 or files.size()>2){
            show_usage(argv[0], "lengthsort");
            return 1;
        }

        if(files.size()==2){
            paired=true;
        }

        if(correct==true and paired==false){
                cerr<<"The -c option can only be used with paired fastq files.\n";
                return -1;
        }


        //##### Check for presence of R ####################################
        string R_path=check_R((char *)"/usr/bin/R", argv[0]);


        // ########## Remove trailing "/" from directory ##############
        if (directory && directory[strlen(directory)-1] == '/'){
            while(directory[strlen(directory)-1] == '/'){
                char* correctdir= new char[strlen(directory)-1];
                strcpy(correctdir, directory);
                correctdir[strlen(directory)-1] = '\0';
                directory=correctdir;
            }
        }
        if (directory){
            size_t end_index = strlen(directory);
            char* correctdir= new char[end_index+2];
            strcpy(correctdir, directory);
            correctdir[end_index] = '/';
            correctdir[end_index+1] = '\0';
            directory=correctdir;
        }

        //###### Open file(s): input(s) and output(s) #######################################

        int count_d=0, count_s=0, count_p1=0, count_p2=0;
        const char* first_filepath_chr;
        std::string first_filename;
        std::string first_filepath;
        std::string single, discard;
        const char* single_chr;
        const char* discard_chr;
        std::string second_filename;
        std::string second_filepath;
        int l1,l2, sz;
        gzFile fp1;
        kseq_t *seq1;
        gzFile fp2;
        kseq_t *seq2;

        // Input SE
        first_filepath = boost::filesystem::absolute(boost::filesystem::path(files[0])).string();
        first_filename = boost::filesystem::path(files[0]).filename().string();
        first_filepath_chr=first_filepath.c_str();
        fp1 = gzopen(first_filepath_chr, "r");
        seq1 = kseq_init(fp1);

        bool zipped=false;
        if(first_filename.substr(first_filename.find_last_of(".") + 1) == "gz"){
            zipped=true;
        }

        // If Windows, zipped files not supported.
        #ifdef _WIN32
            if(zipped==true){
                cout<<"Compressed files are not supported in Windows. Please decompress your fastq files and run again."<<endl;
                return -1;
            }
        #endif // OS_WINDOWS

        if(zipped==false){
            FILE* FIRST;
            if ((FIRST = fopen(first_filepath_chr, "rb")) == NULL){
                printf("Failed to open file \'%s\'\n", first_filepath_chr);
                return -1;
            }
            sz=fsize(FIRST);
        }else{
            sz=zipped_filesize(fp1);
            gzrewind(fp1); kseq_rewind(seq1); // Go to beginning of file
        }
        // Progress bar setup

        int progress_unit=int(sz/100);
        int percent_progress=int(sz/100);
        int progress_counter=0;
        //


        // Output for SE
        if(directory){
                    single=directory+first_filename+".single";
                    single_chr=single.c_str();
                    discard=directory+first_filename+".discard";
                    discard_chr=discard.c_str();
        }else{
                    single=first_filepath+".single";
                    single_chr=single.c_str();
                    discard=first_filepath+".discard";
                    discard_chr=discard.c_str();
        }
        if (is_file_exist(single_chr)){
                cout<<"Error: Single file "<<single_chr<<" already exists."<<endl;
                return -1;
        }

        if (is_file_exist(discard_chr)){
                cout<<"Error: Discard file "<<discard_chr<<" already exists."<<endl;
                return -1;
        }


        ofstream SINGLE (single_chr, std::ofstream::out);
        ofstream DISCARD (discard_chr, std::ofstream::out);


        if(paired==true){
            const char* second_filepath_chr;
            std::string second_filename;
            std::string second_filepath;
            std::string paired1, paired2;
            const char* paired1_chr;
            const char* paired2_chr;
            bool zipped2=false;
            // Input PE
            second_filepath = boost::filesystem::absolute(boost::filesystem::path(files[1])).string();
            second_filename = boost::filesystem::path(files[1]).filename().string();
            second_filepath_chr=second_filepath.c_str();
            if(second_filename.substr(second_filename.find_last_of(".") + 1) == "gz"){
                zipped2=true;
            }
             // If Windows, zipped files not supported.
            #ifdef _WIN32
                if(zipped2==true){
                    cout<<"Compressed files are not supported in Windows. Please decompress your fastq files and run again."<<endl;
                    return -1;
                }
            #endif // OS_WINDOWS

            FILE* SECOND;
            if ((SECOND = fopen(second_filepath_chr, "rb")) == NULL){
                printf("Failed to open file \'%s\'\n", second_filepath_chr);
                return -1;
            }

            fp2 = gzopen(second_filepath_chr, "r");
            seq2 = kseq_init(fp2);

            // Output PE
            if(directory){
                        paired1=directory+first_filename+".paired";
                        paired1_chr=paired1.c_str();
                        paired2=directory+second_filename+".paired";
                        paired2_chr=paired2.c_str();
            }else{
                        paired1=first_filepath+".paired";
                        paired1_chr=paired1.c_str();
                        paired2=second_filepath+".paired";
                        paired2_chr=paired2.c_str();
            }


            if (is_file_exist(paired1_chr)){
                cout<<"Error: Paired file "<<paired1_chr<<" already exists."<<endl;
                return -1;
            }
            if (is_file_exist(paired2_chr)){
                cout<<"Error: Paired file "<<paired2_chr<<" already exists."<<endl;
                return -1;
            }

            ofstream PAIRED1 (paired1_chr, std::ofstream::out);
            ofstream PAIRED2 (paired2_chr, std::ofstream::out);

            std::string header1, header2;

            if(correct==true){
                if (correct_reads(directory, first_filename, second_filename, first_filepath, second_filepath, first_filepath_chr, second_filepath_chr)==-1) return -1;
                gzclose(fp1);
                gzclose(fp2);
                first_filepath=first_filepath+".clean";
                first_filepath_chr=first_filepath.c_str();
                second_filepath=second_filepath+".clean";
                second_filepath_chr=second_filepath.c_str();
                fp1 = gzopen(first_filepath_chr, "r");
                seq1 = kseq_init(fp1);
                fp2 = gzopen(second_filepath_chr, "r");
                seq2 = kseq_init(fp2);

            }
            while ((l1 = kseq_read(seq1)) >= 0 and (l2 = kseq_read(seq2)) >= 0) {
                // Progress bar update
                progress_counter+=strlen(seq1->qual.s)+strlen(seq1->seq.s)+strlen(seq1->name.s);
                if(seq1->comment.s!=NULL){
                    progress_counter+=strlen(seq1->comment.s);
                    }
                if (progress_counter>percent_progress){
                    loadbar_bysize(percent_progress, sz);
                    percent_progress+=progress_unit;
                }
                //

                if (seq1->comment.s) header1=std::string("@")+seq1->name.s+" "+seq1->comment.s;
                else header1=std::string("@")+seq1->name.s;
                //cout<<"Header1: "<<header1<<endl;

                if (seq2->comment.s) header2=std::string("@")+seq2->name.s+" "+seq2->comment.s;
                else header2=std::string("@")+seq2->name.s;
                //cout<<"Header2: "<<header2<<endl;
                string readname1=extract_read_name(header1);
                string readname2=extract_read_name(header2);
                if (readname1!=readname2){
                    cerr<<"Error: reads in files "<<first_filename<<" and "<<second_filename<<" do not seem to be paired.\nPlease re-run LengthSort using the -c flag to remove unpaired reads from fastq files."<<endl;
                    SINGLE.close();
                    DISCARD.close();
                    PAIRED1.close();
                    PAIRED2.close();
                    remove(paired1);
                    remove(paired2);
                    remove(single);
                    remove(discard);

                    return -1;
                }
                if(strlen(seq1->seq.s)>=(unsigned)length and strlen(seq2->seq.s)>=(unsigned)length){
                    PAIRED1<<header1<<endl<<seq1->seq.s<<endl<<"+"<<endl<<seq1->qual.s<<endl;
                    count_p1+=1;
                    PAIRED2<<header2<<endl<<seq2->seq.s<<endl<<"+"<<endl<<seq2->qual.s<<endl;
                    count_p2+=1;
                }else if(strlen(seq1->seq.s)<(unsigned)length and strlen(seq2->seq.s)<(unsigned)length){
                    DISCARD<<header1<<endl<<seq1->seq.s<<endl<<"+"<<endl<<seq1->qual.s<<endl;
                    count_d+=1;
                    DISCARD<<header2<<endl<<seq2->seq.s<<endl<<"+"<<endl<<seq2->qual.s<<endl;
                    count_d+=1;
                }else if(strlen(seq1->seq.s)<(unsigned)length and strlen(seq2->seq.s)>=(unsigned)length){
                    DISCARD<<header1<<endl<<seq1->seq.s<<endl<<"+"<<endl<<seq1->qual.s<<endl;
                    count_d+=1;
                    SINGLE<<header2<<endl<<seq2->seq.s<<endl<<"+"<<endl<<seq2->qual.s<<endl;
                    count_s+=1;
                }else if(strlen(seq1->seq.s)>=(unsigned)length and strlen(seq2->seq.s)<(unsigned)length){
                    SINGLE<<header1<<endl<<seq1->seq.s<<endl<<"+"<<endl<<seq1->qual.s<<endl;
                    count_s+=1;
                    DISCARD<<header2<<endl<<seq2->seq.s<<endl<<"+"<<endl<<seq2->qual.s<<endl;
                    count_d+=1;
                }
            }
            loadbar(100, 100);
            cout<<"\nWriting files..."<<endl;
            gzclose(fp2);
            PAIRED1.close();
            PAIRED2.close();

        }else if(paired==false){
            std::string header1;

            while ((l1 = kseq_read(seq1)) >= 0){

                // Progress bar update
                progress_counter+=strlen(seq1->qual.s)+strlen(seq1->seq.s)+strlen(seq1->name.s)+strlen(seq1->comment.s);
                if (progress_counter>percent_progress){
                    loadbar_bysize(percent_progress, sz);
                    percent_progress+=progress_unit;
                }
                //

                if (seq1->comment.s) header1=std::string("@")+seq1->name.s+" "+seq1->comment.s;
                else header1=std::string("@")+seq1->name.s;

                if(strlen(seq1->seq.s)>=(unsigned)length){
                    SINGLE<<header1<<endl<<seq1->seq.s<<endl<<"+"<<endl<<seq1->qual.s<<endl;
                    count_s+=1;
                }else{
                    DISCARD<<header1<<endl<<seq1->seq.s<<endl<<"+"<<endl<<seq1->qual.s<<endl;
                    count_d+=1;
                }

            }
            loadbar(100, 100);
            cout<<"\nWriting files..."<<endl;

        }

        gzclose(fp1);
        SINGLE.close();
        DISCARD.close();

        // ##### Write summary file ##################################//
        std::string summary;
        const char* summary_chr;
        if(directory){
            summary=directory+first_filename+".summary.txt";
            summary_chr=summary.c_str();
        }else{
            summary=first_filepath+".summary.txt";
            summary_chr=summary.c_str();
        }
        ofstream SUMMARY (summary_chr, std::ofstream::out);

        if(paired==true){
            SUMMARY<<"paired1\t"<<count_p1<<endl<<"paired2\t"<<count_p2<<endl<<"single\t"<<count_s<<endl<<"discard\t"<<count_d<<endl;
        }else{
            SUMMARY<<"single\t"<<count_s<<endl<<"discard\t"<<count_d<<endl;
        }
        SUMMARY.close();

        // #### R graph ############################################################################
        if (R_path!="0"){
            R_length_summary(directory, first_filepath, first_filename, summary, paired, length);
        }


    }
return 0;
}

// print vector
//         for (auto c : qual)
//             std::cout << c << ' ';


// Print a map of type i->vector
//        for(auto it = mean_bytiles.cbegin(); it != mean_bytiles.cend(); ++it)
//        {
//            std::cout << it->first << " "<< it->second[0]<< "\n";
//        }

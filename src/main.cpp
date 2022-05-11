//
//  main.cpp
//  demultiplex_satay
//
//  Copyright © 2022 Jordan Berg. All rights reserved.
//

#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <ctime>
#include <map>
#include <string.h>
#include <vector>
#include "fastqreader.h"
#include "cmdline.h"
//#include "custom.h"

using namespace std;

#define DEMULTIPLEX_SATAY_VER "0.0.1"
typedef map<string, string> BarcodeMap;
typedef map<string, vector<string>> IndexMap;

// Timer functions
clock_t START_TIMER;

clock_t start() {
    return START_TIMER = clock();
}

void stop(clock_t start = START_TIMER) {
    cout
        << endl
        << "Elapsed time:                                  "
        << (clock() - start) / (double)CLOCKS_PER_SEC
        << "s"
        << endl
        << endl;
}

// Source: http://stackoverflow.com/a/236803/248823
void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

// Read in header-less, index-less 2d tsv/txt where column 1 is barcode/index and column 2 is sample label
void readBarcodes(string file_url, vector<string>* index_array, vector<string>* sample_array) {
    cout 
        << "Reading:                                       " 
        << file_url 
        << endl;

    ifstream infile (file_url);
    string line;
    vector<string> a_array, b_array;

    // Read each line of tab-delimited file
    while (getline(infile, line)) {
        vector<string> row_values;
        split(line, '\t', row_values);
        a_array.push_back(row_values.at(0));
        b_array.push_back(row_values.at(1));
    }

    // Check array sizes match
    if (a_array.size() != b_array.size()) {
        cout
            << "Barcode columns lengths must be equal"
            << endl
            << "Column 1: " << a_array.size()
            << "Column 2: " << b_array.size()
            << endl;
        return 1;
    } else {
        cout 
            << "Read in " << a_array.size() << " barcodes"
            << endl;
    }

    // return
    *index_array = a_array;
    *sample_array = b_array;
}

// Source: https://stackoverflow.com/a/33075485/9571488
string reverseComplement(string seq) {
    auto lambda = [](const char c) {
        switch (c) {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        default:
            throw domain_error("Invalid nucleotide.");
        }
    };
    transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
    return seq;
}

bool fuzzyMatch(string s1, string s2, int threshold) {

    if (s1.length() != s2.length()) {
        return false;
    }

    int mismatch;
    mismatch = 0;
    for (i = 0; i < s1.length(); ++i) {
        if (s1.at(i) != s2.at(i)) {
            mismatch++;
        }
    }
    if (mismatch <= threshold) {
        return true;
    } else {
        return false;
    }
}


// =============================== //
// ------------ MAIN ------------- //
// =============================== //
int main(int argc, char* argv[]) {

    // Parse user arguments
    if (argc == 1) {
        cerr
            << "demultiplex_satay: demultiplexing for pooled Satay screen reads"
            << endl
            << "version "
            << DEMULTIPLEX_SATAY_VER
            << endl;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr
            << "demultiplex_satay "
            << DEMULTIPLEX_SATAY_VER
            << endl;
        return 0;
    }

    cmdline::parser cmd;
    //input file - sequencing reads
    cmd.add<string>("reads", 'r', "Input fastq file name", false, ""); 
    //input file - read indices
    cmd.add<string>("index", 'i', "Index fastq file name", false, ""); 
    //input file - read indices
    cmd.add<string>("barcodes", 'b', "File name for barcodes (.tsv or .txt)", false, ""); 
    //threshold for fuzzy searching of read indices
    cmd.add<int>("fuzzy-threshold", 'f', "Length of UMI, default is 1", false, 1); 

    cmd.parse_check(argc, argv);

    if (argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }

    string reads_file;
    reads_file = cmd.get<string>("reads");

    string index_file;
    index_file = cmd.get<string>("index");

    string barcode_file;
    barcode_file = cmd.get<string>("barcodes");

    int fuzzy_threshold;
    fuzzy_threshold = cmd.get<int>("fuzzy-threshold");

    // Print user inputs
    cout
        << "\ndemultiplex_satay v"
        << DEMULTIPLEX_SATAY_VER << endl << endl;
    cout
        << "Provided sequencing reads file name:           "
        << reads_file << endl;
    cout
        << "Provided index reads file name:                "
        << index_file << endl;
    cout
        << "Provided barcode file name:                    "
        << barcode_file << endl;
    cout
        << "Provided fuzzy mapping threshold:              "
        << fuzzy_threshold << endl << endl;
    cout
        << "Demultiplexing reads...\n" << endl;

    // Check the inputs and exit if one doesn't work
    if (reads_file == index_file) {
        cout
            << "Error: Reads and index file names can not be identical"
            << endl;
        return 1;
    }
    if (reads_file == "") {
        cout
            << "Error: Reads file name cannot be blank"
            << endl;
        return 1;
    }
    if (index_file == "") {
        cout
            << "Error: Index file name cannot be blank"
            << endl;
        return 1;
    }
    if (barcode_file == "") {
        cout
            << "Error: Barcode file name cannot be blank"
            << endl;
        return 1;
    }
    if (fuzzy_threshold < 0) {
        cout
            << "Fuzzy mapping threshold cannot be less than 0"
            << endl;
        return 1;
    }

    // Begin processing file
    
    start(); // start elapsed time

    // Read barcode file
    vector<string> barcode_index, barcode_sample;
    readBarcodes(barcode_file, &barcode_index, &barcode_sample);
    int barcode_number;
    barcode_number = barcode_index.size();

    // Populate barcode map
    BarcodeMap barcode_dictionary;
    for (int i = 0; i < barcode_index.size(); ++i) {
        barcode_dictionary[barcode_index[i]] = barcode_sample[i];
    }

    //fstream output(out_file, std::fstream::in | std::fstream::out | std::fstream::trunc); // prepare output file
    // Read index fastq file
    FastqReader reader1 (index_file); // initialize input FASTQ file
    Read* r1 = NULL;

    // Initialize index matrix
    IndexMap index_dictionary;

    // Process reads from FASTQ file to identify read sample indices
    while (true) {
        r1 = reader1.read();
        if (r1 == NULL) {
            break;
        } else {
            string _name, _index;
            r1->toString();
            _name = r1 -> mName;
            _index = r1 -> mSeq.mStr;
                 
            // Fuzzy search index against barcodes for sample labels
            string matched_index, sample_name;

            // Check if reverse complement in barcode dictionary
            string rev_comp;
            rev_comp = reverseComplement(_index);
            
            // Check if fuzzy _index in barcode_dictionary
            string fuzzy_string;
            vector<string> fuzzy_indices, fuzzy_samples;
            for (BarcodeMap::iterator i = barcode_dictionary.begin(); i != barcode_dictionary.end(); ++i) {
                if (fuzzyMatch(_index, i, fuzzy_threshold)) {
                    fuzzy_indices.push_back(i);
                    fuzzy_samples.push_back(barcode_dictionary[i]);
                }
            }

            // Check if fuzzy reverse complement in barcode dictionary
            string fuzzy_rev_comp;
            vector<string> fuzzy_rev_comp_indices, fuzzy_rev_comp_samples;
            for (BarcodeMap::iterator i = barcode_dictionary.begin(); i != barcode_dictionary.end(); ++i) {
                if (fuzzyMatch(rev_comp, i, fuzzy_threshold)) {
                    fuzzy_indices.push_back(i);
                    fuzzy_samples.push_back(barcode_dictionary[i]);
                }
            }

            // Update values
            if (barcode_dictionary.contains(_index)) {
                matched_index = _index;
                sample_name = barcode_dictionary[_index];
            } else if (barcode_dictionary.contains(rev_comp)) {
                matched_index = rev_comp;
                sample_name = barcode_dictionary[rev_comp];
            } else if () {


            } else if () {


            } else {

            }


            vector<string> strVec;
            strVec.push_back(_index);
            strVec.push_back(matched_index);
            strVec.push_back(sample_name);

            // Add read to index record
            //index_dictionary[_name] = strVec;
        }
        delete r1;
    }

    //output.close(); // close output file

    stop(); // stop and print elapsed time
    cout.flush();
    return 0;
}

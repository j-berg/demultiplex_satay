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
#include <cctype>
#include <string.h>
#include <stdio.h>
#include <vector>
#include "fastqreader.h"
#include "cmdline.h"

using namespace std;

// Global variables
const string DEMULTIPLEX_SATAY_VER = "0.0.1";
const string UNASSIGNED_VALUE      = "unassigned";
const string FASTQ_ID_DELIMITER    = " ";
const string FASTQ_DELIMITER       = ".";
const string FASTQ_SUFFIX          = ".fastq";
const int UPDATE_FREQUENCY         = 1000000;

typedef map<string, string> BarcodeMap;
typedef map<string, vector<string>> IndexMap;


// Timer functions
clock_t START_TIMER;

clock_t start() {
    return START_TIMER = clock();
}

void stop(clock_t start = START_TIMER) {
    cout.precision(4);
    cout
        << endl
        << "Processing complete."
        << endl
        << endl
        << "Elapsed time:                                  "
        << (clock() - start) / (double)CLOCKS_PER_SEC
        << "s"
        << endl
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
int readBarcodes(string file_url, vector<string>* index_array, vector<string>* sample_array) {
    cout 
        << "Reading barcodes file..." 
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
    return 0;
}

// Source: https://stackoverflow.com/a/33075485/9571488
string complementSeq(string seq) {
    string c;
    
    for (int i = 0; i < seq.length(); ++i) {
        char base = seq.at(i);
        if (base == 'A') {
            c = c + "T";
        } else if (base == 'T') {
            c = c + "A";
        } else if (base == 'C') {
            c = c + "G";
        } else if (base == 'G') {
            c = c + "C";
        } else {
            c = c + base;
        }
    }

    return c;
}

string reverseSeq(string seq) {
    reverse(seq.begin(), seq.end());
    return seq;
}

string reverseComplement(string seq) {

    string s_seq_c, s_seq_r;
    s_seq_c = complementSeq(seq);
    s_seq_r = reverseSeq(s_seq_c);
    return s_seq_r;
}

bool fuzzyMatch(string s1, string s2, int threshold) {

    if (s1.length() != s2.length()) {
        return 1;
    }

    int mismatch = 0;
    for (int i = 0; i < s1.length(); ++i) {
        if (s1.at(i) != s2.at(i)) {
            ++mismatch;
        }
    }
    if (mismatch <= threshold) {
        return 0;
    } else {
        return 1;
    }
}

// Adapted from: https://stackoverflow.com/a/41369185/9571488
bool isNotAlnum(char c) {
    if (isalnum(c) == 0) {
        return 1;
    } else if (c == '_') {
        return 1;
    } else if (c == '-') {
        return 1;
    } else {
        return 0;
    }
}

string stripSpecial(string s) {
    string stripped;
    for (int i = 0; i < s.length(); ++i) {
        if (isalnum(s.at(i))) {
            stripped = stripped + s.at(i);
        } else if (s.at(i) == '_') {
            stripped = stripped + s.at(i);
        } else if (s.at(i) == '-') {
            stripped = stripped + s.at(i);
        } else {
            
        }
    }
    return stripped;
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
    cmd.add<string>("reads", 'r', "Input fastq file name", true); 
    //input file - read indices
    cmd.add<string>("index", 'i', "Index fastq file name", true); 
    //input file - read indices
    cmd.add<string>("barcodes", 'b', "File name for barcodes (.tsv or .txt)", true); 
    //threshold for fuzzy searching of read indices
    cmd.add<int>("fuzzy-threshold", 'f', "Length of UMI, default is 1", false, 1); 

    // Parse arguments
    cmd.parse_check(argc, argv);

    if (argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }

    string reads_file = cmd.get<string>("reads");
    string index_file = cmd.get<string>("index");
    string barcode_file = cmd.get<string>("barcodes");
    int fuzzy_threshold = cmd.get<int>("fuzzy-threshold");

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
        << "--------------------------------------------------------------"
        << endl
        << endl;
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
    vector<string> barcodeindex, barcode_sample;
    int success = readBarcodes(barcode_file, &barcodeindex, &barcode_sample);
    if (success == 1) {
        return 1;
    }
    int barcode_number = barcodeindex.size();

    // Populate barcode map
    BarcodeMap barcode_dictionary;
    for (int i = 0; i < barcodeindex.size(); ++i) {
        barcode_dictionary[barcodeindex[i]] = stripSpecial(barcode_sample[i]);
    }

    // Read index fastq file
    FastqReader reader1 (index_file); // initialize input FASTQ file
    Read* r1 = NULL;

    // Initialize index matrix
    IndexMap index_dictionary;

    // Process reads from index FASTQ file to identify read sample indices
    int counter_index = 0;
    int counter_matched_index = 0;
    cout 
        << endl
        << "Reading index file..."
        << endl;
    while (true) {
        
        r1 = reader1.read();

        if (r1 == NULL) {
            break;
        } else {
            ++counter_index;
            if (counter_index % UPDATE_FREQUENCY == 0) {
                cout 
                    << counter_index
                    << " reads processed"
                    << endl;
            }

            r1 -> toString();
            string full_name = r1 -> mName;
            string name = full_name.substr(0, full_name.find(FASTQ_ID_DELIMITER));
            string index = r1 -> mSeq.mStr;

            // Fuzzy search index against barcodes for sample labels
            string matched_index, sample_name;

            // Check if reverse complement in barcode dictionary
            string rev_comp = reverseComplement(index);

            // Check if fuzzy index in barcode_dictionary
            string fuzzy_string;
            vector<string> fuzzy_indices, fuzzy_samples;
            for (auto it = barcode_dictionary.begin(); it != barcode_dictionary.end(); ++it) {
                auto key = it -> first;
                auto value = it -> second;
                if (fuzzyMatch(index, key, fuzzy_threshold) == 0) {
                    fuzzy_indices.push_back(key);
                    fuzzy_samples.push_back(value);
                }
            }

            // Check if fuzzy reverse complement in barcode dictionary
            string fuzzy_rev_comp;
            vector<string> fuzzy_rev_comp_indices, fuzzy_rev_comp_samples;
            for (auto it_f = barcode_dictionary.begin(); it_f != barcode_dictionary.end(); ++it_f) {
                auto key_f = it_f -> first;
                auto value_f = it_f -> second;
                if (fuzzyMatch(rev_comp, key_f, fuzzy_threshold) == 0) {
                    fuzzy_rev_comp_indices.push_back(key_f);
                    fuzzy_rev_comp_samples.push_back(value_f);
                }
            }

            // Update values
            if (barcode_dictionary.count(index) > 0) {
                matched_index = index;
                sample_name = barcode_dictionary[index];
                ++counter_matched_index;
            } else if (barcode_dictionary.count(rev_comp) > 0) {
                matched_index = rev_comp;
                sample_name = barcode_dictionary[rev_comp];
                ++counter_matched_index;
            } else if (fuzzy_indices.size() == 1) {
                // add if fuzzy match array is 1
                matched_index = fuzzy_indices.at(0);
                sample_name = fuzzy_samples.at(0);
                ++counter_matched_index;
            } else if (fuzzy_rev_comp_indices.size() == 1) {
                // add if fuzzy match array is 1
                matched_index = fuzzy_rev_comp_indices.at(0);
                sample_name = fuzzy_rev_comp_samples.at(0);
                ++counter_matched_index;
            } else {
                matched_index = UNASSIGNED_VALUE;
                sample_name = UNASSIGNED_VALUE;
            }

            vector<string> strVec;
            strVec.push_back(index);
            strVec.push_back(matched_index);
            strVec.push_back(sample_name);
            index_dictionary[name] = strVec;
        }

        delete r1;
    }

    cout.precision(4);
    cout 
        << "Read "
        << counter_index
        << " total index reads"
        << endl 
        << "Matched "
        << counter_matched_index
        << " index reads"
        << endl
        << "Sample index match rate: "
        << (double)counter_matched_index / (double)counter_index * 100.00
        << "%"
        << endl;

    // Parse out samples from main FASTQ read file

    // Get list of barcode_dictionary + UNASSIGNED_VALUE
    // Open files
    string output_prefix = reads_file.substr(0, reads_file.find_last_of(FASTQ_DELIMITER));

    // Delete existing files of same names
    for (auto it_d = barcode_dictionary.begin(); it_d != barcode_dictionary.end(); ++it_d) {
        auto key = it_d -> first;
        auto value = it_d -> second;
        string file_name = output_prefix + "_" + value + FASTQ_SUFFIX;
        char* file_name_c = const_cast<char*>(file_name.c_str());

        // Delete if exists
        ifstream infile(file_name);
        if (infile.good()) {
            remove(file_name_c);
        }
    }
    
    string this_file, output_destination;

    // Read sequence fastq file
    FastqReader reader2 (reads_file); // initialize input FASTQ file
    Read* r2 = NULL;

    // Process reads from index FASTQ file to identify read sample indices
    int counter_read = 0;
    int counter_matched_read = 0;
    cout 
        << endl
        << "Reading sequence read file..."
        << endl;

    while (true) { 

        r2 = reader2.read();

        if (r2 == NULL) {
            break;
        } else {
            ++counter_read;
            if (counter_read % UPDATE_FREQUENCY == 0) {
                cout 
                    << counter_read
                    << " reads processed"
                    << endl;
            }

            // Check read ID - before space
            // Dictate output file
            // Append read record to file
            r2->toString();
            string full_name = r2 -> mName;
            string name = full_name.substr(0, full_name.find(FASTQ_ID_DELIMITER));

            if (index_dictionary.count(name) > 0) {
                output_destination = output_prefix + "_" + index_dictionary[name].at(2) + FASTQ_SUFFIX; // Get matched sample name
                fstream output(output_destination, fstream::in | fstream::out | fstream::app); // prepare output file
                output << r2->toString();
                output.close();
                if (index_dictionary[name].at(2) != UNASSIGNED_VALUE) {
                    ++counter_matched_read;
                }
            } else {
                output_destination = output_prefix + "_" + UNASSIGNED_VALUE + FASTQ_SUFFIX;
                fstream output(output_destination, fstream::in | fstream::out | fstream::app); // prepare output file
                output << r2->toString();
                output.close();
            }

        }  

        delete r2;
    }

    cout.precision(4);
    cout 
        << "Read "
        << counter_read
        << " total sequencing reads"
        << endl 
        << "Matched "
        << counter_matched_read
        << " sequencing reads"
        << endl
        << "Sequencing read match rate: "
        << (double)counter_matched_read / (double)counter_read * 100.00
        << "%"
        << endl
        << endl;

    // Exit
    stop(); // stop and print elapsed time
    cout.flush();
    return 0;
}

#include "io.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream> 
#include <stdexcept>
#include <unordered_set>
#include <vector>

bool checkSeq( const std::string& line )
{
    for ( const auto& l : line ) {
        if ( l != 'A' && l != 'T' && l != 'C' && l != 'G' ) return false;
    }
    return true;
}

fastaRecord readFasta(const std::string& filename)
{
    std::ifstream file(filename);
    if ( !file ) throw std::runtime_error("Failed to open FASTA file: " + filename);
    fastaRecord records;
    std::string header, line;
    std::ostringstream seq_stream{};
    bool has_header = false; 
    
    while ( getline(file, line) ) {
        if ( line.empty() ) continue;
        if ( line[0] == '>' ) {
            if (has_header) {
                records.rec[header] += seq_stream.str();
                seq_stream.str("");
                seq_stream.clear();
            }
            header = line;
            // If the current header is not the most recent entry it must be a duplicate
            if ( records.rec.find(header) != records.rec.end() ) {
                throw std::runtime_error("The header '" + header + "' appears more than once.");
            }
            has_header = true;
        } 
        else {
            if (!has_header) throw std::runtime_error("The file must begin with a header line.");
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq_stream << line;
        }
    }
    // Store the final sequence
    if (has_header) {
        records.rec[header] = seq_stream.str();
    }

    return records;
}

// Output the probe sequences in FASTA format
void panelOut(const fastaRecord& probePanel, const std::filesystem::path& outdir)
{
    std::ofstream file(outdir.string() + "/probes.fa");
    std::unordered_set<std::string> seen;
    seen.reserve(probePanel.rec.size());
    int dup = 0;

    for (const auto& [id, seq] : probePanel.rec) {
        // Write to file if insert into unordered set is successful (seq is unique).
        if (seen.insert(seq).second) {
            if (!checkSeq(seq)) {
                std::cout << id << " contains a non-nucleotide character. Skipping.\n";
                continue;
            } 
            file << id << "\n"
                 << seq << "\n";
        } else {
            dup++;
        }
    }
    file.close();
    std::cout << dup << " duplicate probe(s).\n";
    std::cout << "\nUnique probe sequences written to " << outdir.string() + "/probes.fa\n";

    return;
}
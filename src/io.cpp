#include "io.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

bool checkSeq( const std::string& line )
{
    for ( const auto& l : line )
    {
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

    while ( getline(file, line) )
    {
        if ( line.empty() ) continue;
        if ( line[0] == '>' ) 
        {
            header = line;
            // If the current header is not the most recent entry it must be a duplicate
            if ( records.rec.find(header) != records.rec.end() )
            {
                throw std::runtime_error("The header '" + header + "' appears more than once.");
            }
        } 
        else 
        {
            if (header.empty()) throw std::runtime_error("The file must begin with a header line.");
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            if (!checkSeq(line)) throw std::runtime_error("Invalid character in " + header + " sequence.");
            records.rec[header] += line;
        }
    }

    return records;
}

// Output the probe sequences in FASTA format
void panelOut(const fastaRecord& probePanel, const std::filesystem::path& outdir)
{
    std::ofstream file(outdir.string() + "/probes.fa");
    for (const auto& [id, seq] : probePanel.rec)
    {
        file << id << "\n"
             << seq << "\n";
    }
    file.close();
    std::cout << "Probe sequences written to " << outdir.string() + "/probes.fa\n";

    return;
}
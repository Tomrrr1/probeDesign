#include "io.h"
#include "probes.h"
#include <iostream>

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 7)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fasta_file> <probe_len> <offset> <mode> <regex> <outdir>\n";
        return -1;
    }
    std::string filename = argv[1];
    int offset = std::stoi(argv[2]);
    int probe_len = std::stoi(argv[3]);
    char mode = argv[4][0];
    std::string reg = argv[5];
    std::filesystem::path outdir = argv[6];
    
    try
    {
        if (!std::filesystem::is_directory(outdir))
        {
            throw std::runtime_error(outdir.string() + " does not exist. Please enter a valid path");
        }
        fastaRecord faRecords = readFasta(filename);
        std::cout << "Total records read = " << faRecords.rec.size() << '\n';
        fastaRecord updateFa = reverseSeq(faRecords, reg);
        fastaRecord probePanel = designProbe(updateFa, probe_len, offset, mode);
        panelOut(probePanel, outdir);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return -1;
    }
    
    return 0;
}
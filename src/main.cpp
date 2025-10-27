#include "io.h"
#include "probes.h"
#include <boost/program_options.hpp>
#include <cstring>
#include <iostream>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    std::string fasta, regex;
    std::filesystem::path outdir;
    int probe_len, offset;
    char mode{};

    po::options_description desc("\nAllowed options");
    desc.add_options()
        ("in,i", po::value<std::string>(&fasta), "Path to a multi-sequence fasta file.")
        ("length,l", po::value<int>(&probe_len)->default_value(100), "Probe length.")
        ("offset,b", 
            po::value<int>(&offset)->default_value(0), 
            "Ignore the first and/or last N nucleotides of the target sequence in probe creation.")
        ("mode,m", po::value<char>(&mode)->required(), 
            "Probe creation mode:\n"
                "   5: \tsingle 5' probe per sequence at [offset, offset + probe_len).\n"
                "   3: \tsingle 3' probe per sequence ending offset bases from the 3' end.\n"
                "   a: \tboth 5' and 3' probes. If probes overlap the sequences are skipped.\n"
                "   t: \ttile each sequence with probes, skipping regions with extreme GC content.")
        ("regex,r", 
            po::value<std::string>(&regex), 
            "Regular expression. Sequences with a matching header are reverse-complemented before probe design.")
        ("out,o", 
            po::value<std::filesystem::path>(&outdir)->default_value("./"), 
            "Directory where 'probes.fa' will be written.")
        ("help,h", "print programme options.");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) 
    {
        std::cout << desc << "\n";
        return -1;
    }
    po::notify(vm);
    
    try
    {
        if (!std::filesystem::is_directory(outdir))
        {
            throw std::runtime_error(outdir.string() + " does not exist. Please enter a valid path");
        }
        fastaRecord faRecords = readFasta(fasta);
        std::cout << "Total records read = " << faRecords.rec.size() << '\n';
        fastaRecord updateFa = reverseSeq(faRecords, regex);
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
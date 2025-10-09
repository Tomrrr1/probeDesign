// A simple C++ programme for designing DNA capture probes given a multi-sequence FASTA file
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>

struct fastaRecord
{
    std::map<std::string, std::string> rec;
};

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
            records.rec[header] += line;
        }
    }

    return records;
}

// Reverse complement sequences that have an identifying regex in their name
fastaRecord reverseSeq(fastaRecord fa, std::string reg)
{
    std::regex exp(reg);
    // Manual iterator so that we can safely modify the container in the loop
    for ( auto it = fa.rec.begin(); it != fa.rec.end(); ) 
    {
        std::string id = it->first;
        std::string seq = it->second;
        if ( std::regex_search(id, exp) )
        {
            std::string rev;
            for (auto it = seq.rbegin(); it != seq.rend(); it++)
            {
                char base = *it;
                switch (base)
                {
                    case 'A':
                        rev.append(1, 'T');
                        break;
                    case 'T': 
                        rev.append(1, 'A');
                        break;
                    case 'C':
                        rev.append(1, 'G');
                        break;
                    case 'G':
                        rev.append(1, 'C');
                        break;
                }
            }
            auto nh = fa.rec.extract(it++); // nh owns the node so it can be modified
            nh.key() = id + "_REORIENTED";
            nh.mapped() = rev;
            fa.rec.insert(std::move(nh));

        } 
        else
        { 
            it++;
            continue;
        }
    }
    
    return fa;
}

// Design probes of length X matching the 5' and 3' end of each sequence
fastaRecord designProbe(const fastaRecord& fa, int probe_length, int end_offset, char probe_pos)
{
    if (end_offset < 0) throw std::runtime_error("end_offset must be a non-negative integer.");
    fastaRecord probePanel;
    std::string new_id5, new_id3, probe5, probe3;

    for (const auto& [id, seq] : fa.rec )
    {
        if ( seq.size() > INT_MAX ) 
        {
            throw std::runtime_error("The size of each sequence must be less than INT_MAX.");
        }
        if (probe_length > seq.length())
        {
            std::cerr << "probe_length is longer than " << id 
                      << ". Skipping probe creation." << "\n";
            continue;
        }
        switch (probe_pos)
        {
            case '5':
            {
                new_id5 = id + "_5_PROBE";
                probe5 = seq.substr(end_offset, probe_length);
                probePanel.rec[new_id5] = probe5;
                break;
            }
            case '3':
            {
                new_id3 = id + "_3_PROBE";
                probe3 = seq.substr(static_cast<int>(seq.size()) - end_offset - probe_length, probe_length);
                probePanel.rec[new_id3] = probe3;
                break;
            }
            case 'a':
            {
                if (2*probe_length > seq.size() - 2*end_offset)
                {
                    std::cerr << "Unable to design two non-overlapping probes for  " << id
                              << ". Skipping probe creation." << "\n";
                    continue;
                }
                new_id5 = id + "_5_PROBE";
                probe5 = seq.substr(end_offset, probe_length);
                new_id3 = id + "_3_PROBE";
                probe3 = seq.substr(static_cast<int>(seq.size()) - end_offset - probe_length, probe_length);
                probePanel.rec[new_id5] = probe5;
                probePanel.rec[new_id3] = probe3;
                break;
            }
            default:
            {
                std::string pp_str{probe_pos};
                throw std::runtime_error(pp_str + " is not a valid probe_pos. Select one of '5', '3' or 'a'");
            }
        }
    }

    return probePanel;
}

// Output the probe sequences in fasta format
void panelOut(fastaRecord probePanel, std::filesystem::path outdir)
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


int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 7)
    {
        std::cerr << "Usage: " << argv[0] 
                  << " <fasta_file> <end_offset> <probe_length> <probe_pos> <sequence_tag> <outdir>\n";
        return -1;
    }
    std::string filename = argv[1];
    int end_offset = std::stoi(argv[2]);
    int probe_length = std::stoi(argv[3]);
    char probe_pos = argv[4][0]; // Either 5 (probe just at 5'), 3 (probe at 3') or a (probe at both ends).
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
        fastaRecord probePanel = designProbe(updateFa, end_offset, probe_length, probe_pos);
        panelOut(probePanel, outdir);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return -1;
    }
    return 0;
}

#include "probes.h"
#include <iostream>
#include <regex>

// Function to take the reverse complement of a FASTA sequence based on an identifying header tag
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
            for ( auto it = seq.rbegin(); it != seq.rend(); it++ )
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

// Function to tile the sequence with probes
fastaRecord probeTile(const std::string& seq, const std::string& id, int probe_len, int offset)
{
    fastaRecord probePanel;
    int l = 0;
    int r = probe_len - 1;
    while ( r < (int)seq.length() - probe_len - offset )
    {
        std::string sub = seq.substr(l, probe_len);
        double gc = gcContent(sub);
        if (gc > 0.4 && gc < 0.6)
        {
            std::string new_id = id + "_" + std::to_string(l) + "_" + std::to_string(r);
            probePanel.rec[new_id] = sub;
            l += probe_len - 1;
            r += probe_len - 1;
        } 
        else
        {
            l++;
            r++;
        }
    }

    return probePanel;
}

// Helper function to compute the GC content of the current window
double gcContent(const std::string& seq)
{
    double n_g = 0.0;
    double n_c = 0.0;
    for ( char base : seq )
    {
        switch (base)
        {
            case 'G':
                n_g+=1.0;
                break;
            case 'C':
                n_c+=1.0;
                break;
            default:
                break;
        }
    }

    return (n_g + n_c) / seq.length();
}

// Design probes
fastaRecord designProbe(const fastaRecord& fa, int probe_len, int offset, char mode)
{
    if (offset < 0) throw std::runtime_error("offset must be a non-negative integer.");
    fastaRecord probePanel;
    std::string new_id5, new_id3, probe5, probe3;

    for (const auto& [id, seq] : fa.rec )
    {
        if ( seq.size() > INT_MAX )
        {
            throw std::runtime_error("The number of nucleotides in each sequence must be less than INT_MAX.");
        }
        if (probe_len > (int)seq.size())
        {
            std::cerr << "probe_len is longer than " << id 
                      << ". Skipping probe creation." << "\n";
            continue;
        }
        switch (mode)
        {
            case '5':
            {
                new_id5 = id + "_5_PROBE";
                probe5 = seq.substr(offset, probe_len);
                probePanel.rec[new_id5] = probe5;
                break;
            }
            case '3':
            {
                new_id3 = id + "_3_PROBE";
                probe3 = seq.substr((int)seq.size() - offset - probe_len, probe_len);
                probePanel.rec[new_id3] = probe3;
                break;
            }
            case 'a':
            {
                if (2*probe_len > (int)seq.size() - 2*offset)
                {
                    std::cerr << "Unable to design two non-overlapping probes for  " << id
                              << ". Skipping probe creation." << "\n";
                    continue;
                }
                new_id5 = id + "_5_PROBE";
                probe5 = seq.substr(offset, probe_len);
                new_id3 = id + "_3_PROBE";
                probe3 = seq.substr((int)seq.size() - offset - probe_len, probe_len);
                probePanel.rec[new_id5] = probe5;
                probePanel.rec[new_id3] = probe3;
                break;
            }
            case 't':
            {
                fastaRecord tiled = probeTile(seq, id, probe_len, offset);
                probePanel.rec.insert(tiled.rec.begin(), tiled.rec.end());
                break;
            }
            default:
            {
                std::string pp_str{mode};
                throw std::runtime_error(pp_str + " is not a valid mode. Select one of '5', '3', 'a' or 't'");
            }
        }
    }

    return probePanel;
}
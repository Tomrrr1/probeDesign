// io.h
#ifndef IO_H
#define IO_H

#include <filesystem>
#include <map>
#include <string>

struct fastaRecord
{
    std::map<std::string, std::string> rec;
};

bool checkSeq( const std::string& line );
fastaRecord readFasta(const std::string& filename);
void panelOut(const fastaRecord& probePanel, const std::filesystem::path& outdir);

#endif
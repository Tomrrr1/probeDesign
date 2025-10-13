#ifndef PROBES_H
#define PROBES_H

#include "io.h"
#include <string>

fastaRecord reverseSeq(fastaRecord fa, std::string reg);
fastaRecord probeTile(const std::string& seq, const std::string& id, int probe_len, int offset);
double gcContent(const std::string& seq);
fastaRecord designProbe(const fastaRecord& fa, int probe_len, int offset, char probe_pos);

#endif
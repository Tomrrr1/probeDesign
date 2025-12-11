#ifndef PROBES_H
#define PROBES_H

#include "io.h"
#include <string>

fastaRecord probeTile(const std::string& seq, const std::string& id, int probe_len, int spacing, int offset);
double gcContent(const std::string& seq);
fastaRecord designProbe(const fastaRecord& fa, int probe_len, int offset, char mode, int spacing);

#endif
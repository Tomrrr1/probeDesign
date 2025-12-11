# probeDesign

A C++ programme for designing DNA probes for hybridisation-capture experiments. 

# Installation

To get and compile the latest version, run:

    git clone https://github.com/Tomrrr1/probeDesign.git
    cd probeDesign
    make BOOST_ROOT=/path/to/boost/

# Usage

```
probeDesign <in> <length> <offset> <spacing> <mode> <out>
```

- **in**: Path to a multi‑sequence FASTA file.
- **length**: Probe length.
- **offset**: Ignore the first and/or last N nucleotides in probe creation.
- **spacing**: The spacing between probes. Only used in 't' mode.
- **mode**: One of:
  - 5 — single 5' probe per sequence at `[offset, offset + probe_len)`.
  - 3 — single 3' probe per sequence ending `offset` bases from the 3' end.
  - a — both 5' and 3' probes. If probes overlap, the sequences are skipped.
  - t — tile each sequence with probes, skipping regions with extreme GC content.
- **out**: directory where `probes.fa` will be written.

# nanopore-pipeline

This script strings together the following softwares:

- ONT Albacore Sequencing Pipeline Software (version 2.3.1)
- [Porechop (0.2.0)](https://github.com/rrwick/Porechop)
- [NanoPlot (1.13.0)](https://github.com/wdecoster/NanoPlot)
- [LAST-926](http://last.cbrc.jp/)
- [minimap2 (2.10-r761)](https://github.com/lh3/minimap2)
- [NanoSV (1.1.2)](https://github.com/mroosmalen/nanosv)

The script is based off [David Coffey's variant](https://github.com/davidcoffey/MinION) with the following changes:

1. Nanoplot, LAST, minimap2, and NanoSV all run on every barcode separated bin created by porechop.

2. Added options to choose LAST vs minimap2 vs both.

3. New organization: a folder is created with the name of the run (taken as input) and the timestamp of when the script was started. Directories within are also labeled accordingly.

4. The settings for each software have been optimized to run over 12 cores with 250gb memory.

5. Parallel processing

## Parallelization

The biggest improvements made to this pipeline have been in running the programs after porechop's demultiplexing in parallel. In experimenting with the effects of running post-porechop programs, runtimes have been greatly reduced. All of the following runs were performed on a 3 Gbp read with the programs NanoPlot, LAST/NanoSV, and minimap2/NanoSV:

- The [basic pipeline script](scripts/script_sv_barcode_separated_v2.sh), which runs with no background processes on every demuxed barcode, takes **235m25.799s**.
- An [individually parallelized version](scripts/previous-versions/script_sv_parallelized_v1.1.sh), where each program is run separately but all the barcodes are run in parallel within each software, takes **69m0.263s**.
- A [fully parallel version](scripts/script_sv_parallelized_v2.0.sh), where the programs are all run in parallel and are each individually parallelized among themselves, takes **56m46.236s**.

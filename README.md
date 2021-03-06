# sc_read_kallisto_wrapper

License
========

sc_read_kallisto_wrapper
https://github.com/vibbits/sc_read_kallisto_wrapper
Copyright (c) BITS VIB 2017

This program is free software; 
you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; 
either version 3 of the License, or any later version.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of the VIB nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

Full license text can be found in LICENSE.txt

Usage
========

The sc_read_kallisto_wrapper takes as input a BCL file generated by an
Illumina sequencer run on a 10XGENOMICS Chromium single cell 3' RNA-seq
library, extracts the reads, maps them to a transcriptome using kallisto
and formats the equivalence class counts as input for Seurat.

It is necessary to install the latest versions of CellRanger, bcl2fastq and kallisto.

The script sc_read_kallisto_wrapper.pl must, before it can be used, be
modified. There is a block of code "These must be adapted appropriately". The
shebang "#!/usr/bin/perl" at the top might have to be modified as well.

Type ./sc_read_kallisto_wrapper.pl -h to obtain information about how to format
the command line.
There is a stand-alone HTML file sc_read_kallisto_wrapper.html with more information about what the wrapper does.
You can find the same information in the [wiki](https://github.com/vibbits/sc_read_kallisto_wrapper/wiki)

**What if I have no BCL file ?**  
  
If you have no access to the BCL file because the lab sent you directly the fastQ files, you can use instead sc_read_kallisto_simplewrapper.pl. Make sure the fastQ files respect the naming convention. The files must have names  
  
**test_sample_S1_L001_R1_001.fastq.gz**  
where  
* **test_sample** is the sample name  
* **S1** is S followed by a sample number  
* **R1** is R or I followed by a number

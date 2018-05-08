"""
=============================
Metagenome filtering pipeline
=============================

:Author: Matt Jackson
:Release: $Id$
:Date: |today|
:Tags: Python

Filters NGS data from metagenomics experiments to remove one or both of rRNA and host contaminant reads.
Optionally retaining filtered reads.

Overview
========

The pipeline uses the SortMeRNA package to filter rRNAs using a given database.
Host contaminant reads are removed by mapping to a given genome and subsequently filtering
mapped reads.

  
Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 


Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <filename>.<suffix>

The ``suffix`` determines the file type.  The
following suffixes/file types are possible:

fasta,fastq

fasta.1,fastq.1

fasta.gz,fastq.qz

fasta.1.gz,fastq.1.gz

Where unnumbered files relate to single end reads or inter-leaved paired-end files.
Interleaved files should be detected automatically.
Files containg 1 in the suffix will be assumed paired end and should be in the same 
directory as the read 2 files, which share the same file name except for the read number.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
------------

On top of the default CGAT setup and the MetaAssemblyKit pipeline, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The main output is fasta/fastq files filtered to remove the chosen rRNA/host reads.
Optionally, the filtered reads can also be retained.

Glossary
========

.. glossary::

Code
====

"""

#load modules
from ruffus import *
import os
import sys

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
       ["%s.ini" % __file__[:-len(".py")],
       "../pipeline.ini",
       "pipeline.ini" ] )
PARAMS = P.PARAMS

#add PipelineMetaAssemblyKit from ini file
sys.path.insert(0, PARAMS["General_metaassembly_path"])
print(sys.path)
import PipelineMetaAssemblyKit


#get all files within the directory to process
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz", "*.fasta.1",
                 "*.fna", "*.fna.gz", "*.fna.1.gz", "*.fna.1",
                 "*.fa", "*.fa.gz", "*.fa.1.gz", "*.fa.1", 
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz","*.fastq.1")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)")


###################################################
# Check the input file, count number of reads
###################################################
@follows(mkdir("filesummaries.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"filesummaries.dir/\1.seqsummary")
def checkFile(infile, outfile):
    seqdat=PipelineMetaAssemblyKit.SequencingData(infile)
    outf=open(outfile,'w')
    outf.write("name\t{}\nformat\t{}\ncompressed\t{}\npaired\t{}\ninterleaved\t{}\n".format(
        seqdat.filename,seqdat.fileformat,seqdat.compressed,seqdat.paired,seqdat.interleaved))
    seqdat.readCount()
    outf.write("read_count\t{}\n".format(seqdat.readcount))
    outf.close()

@follows(checkFile)
def full():
    pass
    
if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))

.. _rna_pipeline:

RNA Pipeline Quickstart
=======================

1) Download and install the gemtools
------------------------------------
If that goes well, you will have a "gemtools" command line tool available.
Please check the `GEMTools homepage <http://gemtools.github.io/>`_ for download
and installation instructions.

2) Create a genome index
------------------------
In order to run the pipeline for any given genome, you need to create a gem 
index for that particular genome. The gem-indexer takes a single single fasta 
file as input. Assume you have a genome.fa file, you can create the .gem index
by calling::

    gemtools index -i genome.fa -t 8

Here ``-t 8`` indicates to use 8 threads/cpus, put this to the number of cpus 
available to you. The command will create 3 files just next to the fasta file::

    genome.gem  -- the gem index
    genome.hash -- hash version of the genome. You probably don't need this 
                   and the creation can be avoided using --no-hash
    genome.log  -- log file of the indexing process

3) Transcriptome index
---------------------------------
If you have a genome annotation in GTF format available, you should
create a transriptome index from that GTF file. Say you have an ``annotation.gtf``, 
file. You can create a transcriptome index like this::

    gemtools t-index -a annotation.gtf -i genome.gem -t 8

The transcriptome indexer takes the annotation in GTF format and the previously
generated .gem index as input. The ``-t 8`` again indicates that 8 threads should
be used by the process.

This creates 5 files in the current folder:

    annotation.gtf.junctions      -- the splice junctions 
    annotation.gtf.junctions.fa   -- the transcriptome
    annotation.gtf.junctions.gem  -- the transcriptome index
    annotation.gtf.junctions.keys -- keys to translate from transcriptome to genome
    annotation.gtf.junctions.log  -- indexer log file

The pipeline needs all the files except the log file and its searching for them
automatically just next to the annotation fiel (if you do not specify the paths
explicitly). For a quick start, just keep the files next to the annotation.

4) Run the pipeline
-------------------
With the index and the optional transcript index you can run the pipeline.
For this example, lets assume you have paired-end reads in two files 
``reads_1.fastq.gz`` and ``reads_2.fastq.gz``.

Run this command to get an overview of what will happen::

    gemtools rna-pipeline -i genome.gem -a annotation.gtf -f reads_1.fastq -q 33 -t 8 --dry

The tool will complain if anything is missing, otherwise it will print an
overview of the pipeline steps. Note that we assume here that everything is in
the same folder for the sake of simplicity. The main paramters are::


    -i genome.gem     -- the gem genome index
    -a annotation.gtf -- the annotation. You can just skip this if you don't 
                         have one. Otherwise the tool will search for the transcriptome 
                         index just next to the annotation. If you put it 
                         somewhere else you have to explicitly set the paths. Call
                         gemtools rna-pipeline --help for an overview of the available
                         options.

    -f reads_1.fastq  -- the input file. NOTE that by default we assume you have 
                         paired end reads and we search for a second file called 
                         reads_2.fastq.gz (works also for other variants like 
                         reads.0.fastq, reads_0.fq). If it doesn't find the 
                         second file, specify it explicitly like 
                         -f reads_1.fastq.gz reads_2.fastq.gz. If you do NOT 
                         HAVE PAREID reads, just specify the input file and 
                         add --single-end to the paramters.

    -q 33             -- this is the quality offset. Should be 33 or 64, but 
                         you have to specify this as figuring this out automatically
                         can be expensive.
    -t 8              -- threads again. This is important as it significantly 
                         speeds up the runs!

Now, if everything looks good, start the same command but remove the ``--dry`` 
to actually start the run. With the default configuration you will get these 
output files::

    reads.map.gz    -- the gem aligned reads
    reads.bam       -- alignments in bam format
    reads.junctions -- the denovo junction sites found
    reads.stats*    -- two sets of stats, *all* for all the mappings found, 
                       *best* considering only the best mappings. You get two 
                       files each. .txt is the human readable form, .json is 
                       in JSON format so you can easily read the stats with 
                       more or less any modern programming language.

5) Run a stats report on the output
-----------------------------------
This step is optional, but allows a quick result check::

    gemtools report -i reads.stats.all.json -p

The report command will create a .zip file with a graphical report (-p 
indicates paired end reads). Unzip the file and open the index.html in 
your web browser.

Notes
-----
- the pipeline is restartable. When it failed at some point you can just try to restart the run, it will skip the already completed steps.
- you can save pipeline configuration and reuse them later. This is handy if you have more datasets. For example::

    gemtools rna-pipeline -i genome.gem -a annotation.gtf -q 33 -t 8 --save config.gt

  will not run anything but save the configuration to a ``config.gt`` file. 
  You can reuse the configuration like this::

    gemtools rna-pipeline --load config.gt -f reads_A_1.fastq.gz
    gemtools rna-pipeline --load config.gt -f reads_B_1.fastq.gz
    gemtools rna-pipeline --load config.gt -f reads_C_1.fastq.gz
    ....

  All the configuration is taken from the config.gt file and you override the 
  input file from the command line. Note that the configuration file is stored in
  JSON format. You can read this easily and also create configuration files 
  automatically.

- ``gemtools`` has a lot of options. Take a look with ``gemtools rna-pipeline --help``. Here are a few important ones::

    --name         -- specify an output name, otherwise the input file name 
                      is used as a template
    --compress-all -- if you have very limited disk space, you can tell the 
                      pipeline to compress all intermediate files on the fly. 
                      This costs performance though!
    --direct-input -- in case of limited disk space, add this to stream the 
                      initial data directly into gem rather then creating a 
                      dedicated input file.

*LAST BUT VERY IMPORTANT NOTE*

The pipeline expects to find ``samtools`` installed on the system. Try to get
the lates samtools from their github repository
(https://github.com/samtools/samtools -- clone or download and call make to
build it). The latest version is multi-threaded (i.e. ``samtools view --help`` will
show a ``-@`` paramter). Also, see if you have ``pigz`` installed in the system 
you try to run gemtools on. ``pigz`` is a parallel compressor and the pipeline
makes use of it if it is available. It will speed up compression steps a lot!

Running into problems?
----------------------
Please do not hesitate to contact us if you run into any problems or you would
like to see other features implemented. Please consider using the `GEMTools issue
tracker <https://github.com/gemtools/gemtools/issues?>`_.

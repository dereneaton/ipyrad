#!/usr/bin/env python

"""
In silico digest of a fasta genome file to fastq format data files.
"""

import os
import gzip


class DigestGenome(object):
    """
    Digest a fasta genome file with one or two restriction enzymes to create
    pseudo-fastq files to treat as samples in a RAD assembly. 
    
    # Example:
    dg = ipa.digest_genome(
        fasta="genome.fa", 
        workdir="digested_genomes",
        name="quinoa",
        re1="AATCGG",
        re2="CCGG",
        ncopies=5,
        readlen=150,
        paired=True,        
        )
    dg.run()

    """
    def __init__(
        self, 
        fasta, 
        name="digested", 
        workdir="digested_genomes",
        re1="CTGCAG", 
        re2=None, 
        ncopies=1,
        readlen=150, 
        paired=True, 
        min_size=None, 
        max_size=None,
        ):

        self.fasta = fasta
        self.name = name
        self.workdir = workdir
        self.re1 = re1
        self.re2 = re2
        self.ncopies = ncopies
        self.readlen = readlen
        self.paired = paired
        self.min_size = min_size
        self.max_size = max_size

        # use readlen as min_size if not entered
        if not self.min_size:
            self.min_size = self.readlen
        if not self.max_size:
            self.max_size = 9999999


    def run(self):

        # counter
        iloc = 0

        # open output file for writing
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        handle1 = os.path.join(self.workdir, self.name + "_R1_.fastq.gz")
        handle2 = os.path.join(self.workdir, self.name + "_R2_.fastq.gz")
        if self.paired:
            out1 = gzip.open(handle1, 'w')
            out2 = gzip.open(handle2, 'w')

        # load genome file
        if self.fasta.endswith(".gz"):
            fio = gzip.open(self.fasta)
        else:
            fio = open(self.fasta)

        # read in genome and parse scaffolds
        scaffolds = fio.read().decode().split(">")[1:]

        # iterate over scaffolds
        for scaff in scaffolds:
            name, seq = scaff.split("\n", 1)
            seq = seq.replace("\n", "").upper()

            # digest scaffold into fragments and discard scaff ends
            bits = ["1{}1".format(i) for i in seq.split(self.re1)[1:-1]]

            # digest each fragment into second cut fragment
            if self.re2:
                bits1 = bits
                bits = []
                for fragment in bits1:
                    fbits = fragment.split(self.re2)
                    for fbit in fbits:
                        if fbit:

                            # only keep bit if it has both cut sites
                            if (fbit[0] + fbits[-1]).count("1") == 1:
                                bits.append(fbit.strip("1"))

            # filter fragments
            filtered_bits = []
            for bit in bits:
                blen = len(bit)
                if (blen >= self.min_size) & (blen <= self.max_size):
                    filtered_bits.append(bit)

            # turn fragments into (paired) reads
            read1s = []
            read2s = []
            for fragment in filtered_bits:
                r1 = fragment[:self.readlen]
                r2 = fragment[-self.readlen:]
                if self.paired:
                    read1s.append(r1)
                    read2s.append(r2)
                else:
                    read1s.append(r1)
                    read1s.append(r2)

            # write reads to a file
            fastq_r1s = []
            fastq_r2s = []            
            for ridx in range(len(read1s)):
                read1 = read1s[ridx]
                for copy in range(self.ncopies):
                    fastq = "@{name}_loc{loc}_rep{copy} 1:N:0:\n{read}\n+\n{qual}"
                    fastq = fastq.format(**{
                        'name': self.name,
                        'loc': iloc, 
                        'copy': copy,
                        'read': read1, 
                        'qual': "B" * len(read1),
                    })
                    fastq_r1s.append(fastq)

                if self.paired:
                    read2 = read2s[ridx]
                    for copy in range(self.ncopies):
                        fastq = "@{name}_loc{loc}_rep{copy} 2:N:0:\n{read}\n+\n{qual}"
                        fastq = fastq.format(**{
                            'name': self.name,
                            'loc': iloc, 
                            'copy': copy,
                            'read': read2,
                            'qual': "B" * len(read2),
                        })
                        fastq_r2s.append(fastq)
                iloc += 1

            # write all bits of scaffold to disk
            if fastq_r1s:
                out1.write(("\n".join(fastq_r1s)).encode() + b"\n")
                if self.paired:
                    out2.write("\n".join(fastq_r2s).encode() + b"\n")

        # close handles
        fio.close()
        out1.close()
        if self.paired:
            out2.close()


#!/usr/bin/env python

"""Converter class for writing output files formats from HDF5 input.

"""

# pylint: disable=too-many-return-statements

from typing import TypeVar, List, Dict, Iterator
from pathlib import Path
import h5py
import numpy as np
from numba import njit
from ipyrad.assemble.utils import (
    BTS, chroms2ints, IPyradError, NEXHEADER, NEXCLOSER, STRDICT,
)

Assembly = TypeVar("Assembly")


class Converter:
    """Functions for converting hdf5 arrays into output files."""
    def __init__(self, data: Assembly):
        self.data = data
        self.snps_database: Path = data.output_formats["snps_database"]
        self.seqs_database: Path = data.output_formats["seqs_database"]
        self.snames: List[str] = []
        self.pnames: Dict[str,str] = {}

        # fill ordered names from h5 database 
        with h5py.File(self.snps_database, 'r') as io5:
            self.snames = list(io5["snpsmap"]["names"])

        # get padded names based on lengths.
        longname = max(len(i) for i in self.snames)
        self.pnames = {i: i.ljust(longname + 5) for i in self.snames}

    def run(self, oformat: str) -> Path:
        """Write outputs one the selected output format."""
        if oformat == "p":               # phylip format.
            return self.write_phy()
        if oformat == "n":               # nexus format.
            return self.write_nex()
        if oformat == "G":               # gphocs format.
            return self.write_gphocs()
        if oformat == "s":               # snps + snpsmap outputs
            snpsfile = self.write_snps()
            snpsmapfile = self.write_snps_map()
            return snpsfile, snpsmapfile
        if oformat == "u":               # usnps
            return self.write_usnps()
        if oformat == "k":               # str (structure)
            return self.write_str()
        if oformat == "g":               # genos format
            return self.write_geno()
        if oformat == "t":               # treemix format
            return self.write_treemix()
        raise IPyradError(f"output_format {oformat} not recognized.")

    def write_phy(self):
        """Write entire sequence matrix, names padded, with header.
    
        If PE data the PE inserts have been removed.

        10 100
        taxon_aa       AAAAAAAAAATTTTTTTCCCCCCGGGGGGG
        taxon_bbb      NNNNNNNNNNNNNNNNNCCCCCCGGGGGGG
        taxon_cccc     AAAAAAAAAATTTTTTTCCCCCCGGGGGGG                
        taxon_dd       AAAAAAAAAATTTTTTTNNNNNNNNNNNNN
        """
        # write from hdf5 array
        outpath = Path(self.data.stepdir) / f"{self.data.name}.phy"
        with h5py.File(self.seqs_database, 'r') as io5:

            # load seqarray 
            seqarr = io5['phy']

            with open(outpath, 'w', encoding="utf-8") as out:
                out.write("{} {}".format(*seqarr.shape))               # pylint: disable=no-member
                for idx, name in enumerate(self.snames):
                    seq = seqarr[idx, :].tobytes().decode().upper()    # pylint: disable=no-member
                    pname = self.pnames[name]
                    out.write(f"{pname}{seq}\n")
        return outpath


    def _iter_nex_block(self, chunksize: int = 80) -> Iterator[np.ndarray]:
        """Yields chunks of columns of phy array for nexus interleave."""
        with h5py.File(self.seqs_database, 'r') as io5:

            # load seqarray (this could be chunked, this can be >50Gb)
            seqarr = io5['phy']
            shape = seqarr.shape  # pylint: disable=no-member

            # iterate by blocks until all columns fetched
            for bidx in range(0, shape[1], chunksize):
                # grab a big block of data (nsamples, 100_000)
                yield seqarr[:, bidx:bidx + chunksize]

    def _iter_nex_blocks(self, nblocks: int = 100) -> Iterator[List[str]]:
        """Yields lists of locus strings for writing nex interleaved."""
        blocks = []
        bidx = 0
        for chunk in self._iter_nex_block():
            block = []
            for sidx, name in enumerate(self.snames[::-1]):
                pname = self.pnames[name]
                seq = chunk[sidx].tobytes().decode()
                block.append(f"{pname}{seq}")
            blocks.append("\n".join(block))
            bidx += 1
            if bidx >= nblocks:
                yield blocks
                blocks = []
                bidx = 0
        if blocks:
            yield blocks

    def write_nex(self):
        """Write interleaved sequence matrix, names padded, with header.
        
        #nexus
        begin data;
          dimensions ntax=10 nchar=100;
          format datatype=dna missing=N gap=- interleave=yes;
          matrix
            taxon_aa       AAAAAAAAAATTTTTTTCCCCCCGGGGGGG
            taxon_bbb      NNNNNNNNNNNNNNNNNCCCCCCGGGGGGG
            taxon_cccc     AAAAAAAAAATTTTTTTCCCCCCGGGGGGG                
            taxon_dd       AAAAAAAAAATTTTTTTNNNNNNNNNNNNN

            taxon_aa       AAAAAAAAAATTTTTTTCCCCCCGGGGGGG
            taxon_bbb      NNNNNNNNNNNNNNNNNCCCCCCGGGGGGG
            taxon_cccc     AAAAAAAAAATTTTTTTCCCCCCGGGGGGG                
            taxon_dd       AAAAAAAAAATTTTTTTNNNNNNNNNNNNN
        end;
        begin charset;
          ...
        end;
        """
        # open output file to write to
        outpath = Path(self.data.stepdir) / f"{self.data.name}.nex"
        with open(outpath, 'w', encoding="utf-8") as out:
            out.write(NEXHEADER.format(*shape))

            # print intermediate result and clear
            for block in self._iter_nex_blocks():
                out.write("\n".join(block))
                
            # closer
            out.write(NEXCLOSER)
                
            # # add partition information from maparr
            # maparr = io5["phymap"][:, 2]
            # charsetblock = []
            # charsetblock.append("BEGIN SETS;")
            
            # # the first block
            # charsetblock.append("charset {} = {}-{};".format(
            #     0, 1, maparr[0],
            # ))

            # # all other blocks
            # # nexus is 1-indexed. mapparr is dtype np.uint64, so adding
            # # a standard int results in np.float64, so cast uint64 first
            # for idx in range(0, maparr.shape[0] - 1):
            #     charsetblock.append("charset {} = {}-{};".format(
            #         idx + 1, maparr[idx] + np.uint64(1) , maparr[idx + 1]
            #         )
            #     )
            # charsetblock.append("END;")
            # out.write("\n".join(charsetblock))
        return outpath


    def write_snps(self):
        """
        Write all SNPs as a phylip format file. Users should generally
        use ipa.snps_extracter for more fine tuned SNP sampling/filtering.
        """
        # write from hdf5 array
        outfile = os.path.join(self.data.stepdir, f"{self.data.name}.snps.phy")
        with open(outfile, 'w') as out:
            with h5py.File(self.snps_database, 'r') as io5:
                seqarr = io5['snps']
                nsamples = len(self.snames)

                # write dims
                out.write("{} {}\n".format(nsamples, seqarr.shape[1]))

                # write to disk one row at a time
                for idx in range(io5['snps'].shape[0]):
                    seq = seqarr[idx, :].view("S1")
                    out.write(
                        "{}{}".format(
                            self.data.pnames[self.snames[idx]],
                            b"".join(seq).decode().upper() + "\n",
                        )
                    )
        return outfile


    def write_usnps(self):
        """
        TODO: replace with ipa.snps_extracter
        """
        with open(self.data.outfiles.usnps, 'w') as out:
            with h5py.File(self.snps_database, 'r') as io5:
                # load seqarray
                snparr = io5['snps'][:]
                # snp positions 
                maparr = io5["snpsmap"][:, :2]
                maparr[:, 1] = range(maparr.shape[0])

                # get n unlinked snps
                subs = subsample(maparr)
                nsnps = subs.size

                # option to skip ref
                if self.exclude_ref:
                    nsamples = len(self.data.snames) - 1
                    rstart = 1
                else:
                    nsamples = len(self.data.snames)
                    rstart = 0

                # write dims
                out.write("{} {}\n".format(nsamples, nsnps))

                # write to disk one row at a time
                for idx in range(rstart, snparr.shape[0]):

                    # get all SNPS from this sample
                    seq = snparr[idx, subs].view("S1")

                    out.write(
                        "{}{}".format(
                            self.data.pnames[self.data.snames[idx]],
                            b"".join(seq).decode().upper() + "\n",
                        )
                    )

                ## Write the other unlinked formats
                self.write_ustr(snparr[:, subs])
                genos = io5['genos'][:]
                self.write_ugeno(genos[subs, :])


    def write_ustr(self, snparr):
        """
        TODO: replace with ipa.structure() 
        """
        with open(self.data.outfiles.ustr, 'w') as out:
            # option to skip ref
            if self.exclude_ref:
                rstart = 1
            else:
                rstart = 0
            for idx in range(rstart, snparr.shape[0]):
                # get all SNPS from this sample
                seq = snparr[idx, :].view("S1")
                # get sample name
                name = self.data.pnames[self.data.snames[idx]]
                # get row of data
                snps = snparr[idx, :].view("S1")
                # expand for ambiguous bases
                snps = [BTS[i.upper()] for i in snps]
                # convert to numbers and write row for each resolution
                sequence = "\t".join([STRDICT[i[0]] for i in snps])
                out.write(
                    "{}\t\t\t\t\t{}\n"
                    .format(name, sequence))
                ## Write out the second allele if it exists
                if self.data.params.max_alleles_consens > 1:
                    sequence = "\t".join([STRDICT[i[1]] for i in snps])
                    out.write(
                        "{}\t\t\t\t\t{}\n"
                        .format(name, sequence))


    def write_ugeno(self, genos):
        """
        The .geno output is used in ...
        """
        with open(self.data.outfiles.ugeno, 'w') as out:
            # option to skip ref
            if self.exclude_ref:
                rstart = 1
            else:
                rstart = 0

            genos = genos[:, rstart:]
            snpgenos = np.zeros(genos.shape[:2], dtype=np.uint8)
            snpgenos.fill(9)

            # fill (0, 0)
            snpgenos[np.all(genos == 0, axis=2)] = 2

            # fill (0, 1) and (1, 0)
            snpgenos[np.sum(genos, axis=2) == 1] = 1

            # fill (1, 1)
            snpgenos[np.all(genos == 1, axis=2)] = 0

            # write to file
            np.savetxt(out, snpgenos, delimiter="", fmt="%d")


    def write_snps_map(self):
        """
        Write a map file with linkage information for SNPs file
        This dump is mostly unecessary, since ipa tools can read
        the snpsmap from HDF5, but maybe some users will want to see
        it for their own validation.
        """
        counter = 1
        outfile = os.path.join(self.data.stepdir, f"{self.data.name}.snpsmap.tsv")
        with open(outfile, 'w') as out:
            with h5py.File(self.snps_database, 'r') as io5:
                # access array of data
                maparr = io5["snpsmap"]

                ## write to map file in chunks of 10000
                for start in range(0, maparr.shape[0], 10000):
                    outchunk = []

                    # grab chunk
                    rdat = maparr[start:start + 10000, :]

                    # get chroms
                    if self.data.is_ref:
                        revdict = chroms2ints(self.data, 1)
                        for i in rdat:
                            outchunk.append(
                                "{}\t{}:{}\t{}\t{}\n"
                                .format(
                                    i[0], 
                                    # 1-index to 0-index fix (1/6/19)
                                    revdict[i[3] - 1], i[4],
                                    i[2] + 1,
                                    counter,
                                    #i[4], 
                                )
                            )
                            counter += 1
                    else:    
                        # convert to text for writing
                        for i in rdat:
                            outchunk.append(
                                "{}\tloc{}_snp{}_pos{}\t{}\t{}\n"
                                .format(
                                    i[0], 
                                    i[0] - 1, i[4] - 1, i[2],
                                    i[2] + 1, 
                                    counter,
                                    #i[4],
                                )
                            )
                            counter += 1

                    # write chunk to file
                    out.write("".join(outchunk))
                    outchunk = []
        return outfile
                    

    def write_str(self):
        """
        STRUCTURE OUTPUT
        TODO: replace with ipa.structure... (without requiring structure
        to be installed).
        write data from snps database, resolve ambiguous bases and numeric.
        """
        with open(self.data.outfiles.str, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:
                snparr = io5["snps"]

                # option to skip ref
                if self.exclude_ref:
                    rstart = 1
                else:
                    rstart = 0                

                for idx in range(rstart, len(self.data.snames)):
                    # get sample name
                    name = self.data.pnames[self.data.snames[idx]]
                    # get row of data
                    snps = snparr[idx, :].view("S1")
                    # expand for ambiguous bases
                    snps = [BTS[i.upper()] for i in snps]
                    # convert to numbers and write row for each resolution
                    sequence = "\t".join([STRDICT[i[0]] for i in snps])
                    out.write(
                        "{}\t\t\t\t\t{}\n"
                        .format(name, sequence))
                    ## Write out the second allele if it exists
                    if self.data.params.max_alleles_consens > 1:
                        sequence = "\t".join([STRDICT[i[1]] for i in snps])                            
                        out.write(
                            "{}\t\t\t\t\t{}\n"
                            .format(name, sequence))


    def write_gphocs(self):
        """
        b/c it is similar to .loci we just parse .loci and modify it.
        """
        with open(self.data.outfiles.gphocs, 'w') as out:
            indat = iter(open(self.data.outfiles.loci, 'r'))

            # write nloci header
            out.write("{}\n".format(
                self.data.stats_dfs.s7_loci["sum_coverage"].max()))

            # read in each locus at a time
            idx = 0
            loci = []
            locus = []
            while 1:
                try:
                    line = next(indat)
                except StopIteration:
                    indat.close()
                    break

                # end of locus
                if line.endswith("|\n"):
                    
                    # write stats and locus to string and store
                    nsamp = len(locus)
                    slen = len(locus[0].split()[-1])
                    locstr = ["locus{} {} {}\n".format(idx, nsamp, slen)]
                    loci.append("".join(locstr + locus))

                    # reset locus
                    idx += 1
                    locus = []

                else:
                    locus.append(line)

                if not idx % 10000:
                    out.write("\n".join(loci))
                    loci = []
                    
            # write to file
            if loci:
                out.write("\n".join(loci))


    def write_geno(self):
        """
        .geno output is used by ...
        """
        with open(self.data.outfiles.geno, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:

                # option to skip ref
                if self.exclude_ref:
                    rstart = 1
                else:
                    rstart = 0

                genos = io5["genos"][:, rstart:]
                snpgenos = np.zeros(genos.shape[:2], dtype=np.uint8)
                snpgenos.fill(9)
                
                # fill (0, 0)
                snpgenos[np.all(genos == 0, axis=2)] = 2
                
                # fill (0, 1) and (1, 0)
                snpgenos[np.sum(genos, axis=2) == 1] = 1
                
                # fill (1, 1)
                snpgenos[np.all(genos == 1, axis=2)] = 0
                
                # write to file
                np.savetxt(out, snpgenos, delimiter="", fmt="%d")


    def write_treemix(self):
        """
        # We pass in 'binary="ls"' here to trick the constructor into not
        # raising an error if treemix isn't installed. HAX! (lol)
        """
        import ipyrad.analysis as ipa
        tmx = ipa.treemix(
            data=self.data.snps_database,
            name=self.data.name,
            workdir=self.data.dirs.outfiles,
            imap={i: j[1] for (i, j) in self.data.populations.items()},
            minmap={i: j[0] for (i, j) in self.data.populations.items()},
            binary="ls",
            )
        tmx.write_treemix_file()


    def write_migrate(self):
        """
        Calls the migrate ipa utility to write 
        """
        import ipyrad.analysis as ipa
        mig = ipa.migrate_n(
            data=self.data.outfiles.loci,
            name=self.data.name,
            workdir=self.data.dirs.outfiles,
            imap={i: j[1] for (i, j) in self.data.populations.items()},
            minmap={i: j[0] for (i, j) in self.data.populations.items()},
            )
        mig.write_seqfile()


# -------------------------------------------------------------
# jitted subsample func
# -------------------------------------------------------------

@njit
def subsample(snpsmap):
    """Subsample snps, one per locus, using snpsmap."""
    sidxs = np.unique(snpsmap[:, 0])
    subs = np.zeros(sidxs.size, dtype=np.int64)
    idx = 0
    for sidx in sidxs:
        sites = snpsmap[snpsmap[:, 0] == sidx, 1]
        site = np.random.choice(sites)
        subs[idx] = site
        idx += 1
    return subs


if __name__ == "__main__":

    import ipyrad as ip
    DATA = ip.load_json("/tmp/TEST5.json")
    
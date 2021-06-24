#!/usr/bin/env python

"""
Helper classes for step7 (write_outputs.py) for filtering loci for
the .loci output.
"""

import pandas as pd
import numpy as np
from numba import njit
from ipyrad.assemble.write_outputs_edge_trimming import Edges


AMBIGARR = np.array(list(b"RSKYWM")).astype(np.uint8)


class ChunkProcess:
    """
    Runs the remote step7 func apply_filters_and_trimming() to filter
    and trim loci and return stats.
    """
    def __init__(self, data, chunksize, chunkfile):
        
        self.data = data
        self.chunksize = chunksize
        self.chunkfile = chunkfile

        # returned results
        self.filters = pd.DataFrame(
            index=range(self.chunksize),
            columns=[
                "dups", 
                "minsamp",
                "maxind", 
                "maxvar", 
                "maxshared", 
            ],
            data=False,
        )
        # (R1>, <R1, R2>, <R2), currently only trims external edges.
        self.edges = np.zeros((self.chunksize, 4), dtype=np.uint16)
        # stats storage.
        self.stats = {
            'nbases': 0,
            'sample_cov': {i: 0 for i in self.data.samples},
            'locus_cov': {i: 0 for i in range(1, len(self.data.samples) + 1)},
            'var_sites': {i: 0 for i in range(2000)},
            'pis_sites': {i: 0 for i in range(2000)},
            'var_props': {},
            'pis_props': {},
        }
        # updated each iteration
        self.lidx: int = None
        self.outlist = []


    def run(self):
        """
        Iterates over loci in the chunkfile trimming edges and 
        calculating stats for each.
        """
        # Iterate over loci
        io_chunk = open(self.chunkfile, 'rt')
        io_loci = enumerate(iter(io_chunk.read().split("//\n//\n")))

        # keep track of var and pis proportions for histogram at end
        histos = {"pis": [], "var": []}

        while 1:
            # get the next locus in the chunkfile and skips the reference
            # sample if present in the locus and hackers.exclude_reference
            # otherwise reference is first if ref datatype.
            try:
                self.lidx, chunk = next(io_loci)
                names, nidxs, seqs, dup = parse_locus(self.data, chunk)
            except StopIteration:
                break

            # check duplicates filter
            if dup:
                self.filters.loc[self.lidx, "dups"] = True
                continue

            # check minsamp filter (optionally w/ popmins)
            if self.filter_minsamp_pops(names):
                self.filters.loc[self.lidx, "minsamp"] = True
                continue

            # trim edges and store the trimming lengths
            edges = Edges(self.data, seqs)
            edges.get_edges()
            self.edges[self.lidx] = edges.edges

            # check overzealous edge trimming (saved as a minsamp effect)
            if edges.bad:
                self.filters.loc[self.lidx, "minsamp"] = True
                continue

            # apply the edge trimming 
            trimseq = seqs[:, self.edges[self.lidx][0]:self.edges[self.lidx][3]]

            # check if any samples have been overly trimmed (no data)
            keepmask = np.zeros(len(names), dtype=np.bool_)
            nseq = trimseq.copy()
            nseq[nseq == 45] = 78
            for sidx, _ in enumerate(names):
                if np.all(nseq[sidx] == 78):
                    keepmask[sidx] = False
                else:
                    keepmask[sidx] = True
            trimseq = trimseq[keepmask, :]
            names = [i for i, j in zip(names, keepmask) if j]
            nidxs = [i for i, j in zip(nidxs, keepmask) if j]

            # if trimming dropped all samples, then apply minsamp filter
            # NOTE: pre v.1.0 we dropped the whole locus if any sample
            # was removed here. That conservative approach may be safer...
            if sum(keepmask == 0) == len(names):
                self.filters.loc[self.lidx, "minsamp"] = True
                continue

            # create mask to hide sites from stats counts that are in insert
            masked = None

            # [denovo]: store shift of left edge start position from 
            # alignment, this position is needed for pulling depths in VCF.
            if not self.data.is_ref:

                # what is the leftmost consens edge (not -)
                ishift = [
                    np.where(trimseq != 45)[0].min() 
                    for i in range(trimseq.shape[0])
                ]

                # fill nidxs with nidxs and shift info
                inidxs = []
                for idx, (i, j) in enumerate(zip(nidxs, ishift)):

                    # add to ishift if trimmed region contains indels
                    indshift = (trimseq[idx, j:edges.edges[0]] == 45).size
                    inidxs.append("{}-{}".format(i, j + indshift))
                nidxs = inidxs

                # mask insert in denovo data
                trimseq[:, edges.edges[1]:edges.edges[2]] = 78   # N
                # self.aseqs[:, edges.edges[1]:edges.edges[2]] = 110  # n

            # [ref]: nidx string will be updated in to_locus() with edg
            # for is-ref we need to mask the insert between pairs, and
            # for this we need to not check the 'reference' sample.
            else:
                if self.data.is_pair and self.data.params.min_samples_locus > 1:
                    if self.data.hackers.exclude_reference:
                        inserts = np.all(seqs[1:, :] == 78, axis=0)
                    else:
                        inserts = np.all(seqs[:] == 78, axis=0)
                    masked = seqs[:, np.invert(inserts)]

            # apply filters on edge trimmed reads
            if self.filter_maxindels(trimseq):
                continue

            # get snpstring on trimmed reads
            snparr = self.get_snpsarrs(trimseq)

            # check for max % SNPs
            var_prop, pis_prop, filt = self.filter_maxvars(trimseq, masked, snparr)
            if filt:
                continue

            # store proportions variable for histograms
            histos['pis'].append(pis_prop)
            histos['var'].append(var_prop)

            # check for paralogs (shared SNPs)
            if self.filter_maxshared(trimseq):
                continue

            # store stats, locus passed all filters.
            for name in names:
                self.stats['sample_cov'][name] += 1
            # refcount = int(self.data.hackers.exclude_reference)
            self.stats['locus_cov'][trimseq.shape[0]] += 1
            self.stats['var_sites'][snparr.sum()] += 1
            self.stats['pis_sites'][snparr[:, 1].sum()] += 1
            if masked is not None:
                self.stats['nbases'] += masked.shape[1]
            else:
                self.stats['nbases'] += trimseq.shape[1]

            # convert to .loci format string.
            locus = self.to_locus(trimseq, names, nidxs, snparr, edges.edges)
            self.outlist.append(locus)

        # drop keys where values = 0
        self.stats['var_sites'] = {
            i: j for (i, j) in self.stats['var_sites'].items() if j}
        self.stats['pis_sites'] = {
            i: j for (i, j) in self.stats['pis_sites'].items() if j}

        # If no loci survive filtering then don't write the files
        if self.outlist:
            # write the chunk to tmpdir
            with open(self.chunkfile + ".loci", 'wt') as outchunk:
                outchunk.write("\n".join(self.outlist) + "\n")

            # save filters df for only the loci that were filtered.
            mask = self.filters.sum(axis=1).astype(np.bool_).values
            self.filters.loc[mask, :].to_csv(self.chunkfile + ".csv")

            # calculate histograms for stats
            nice_bins = [0, 0.001] + [round(i, 2) for i in np.linspace(0.01, 0.25, 25)]
            mags, bins = np.histogram(histos['var'], bins=nice_bins)
            self.stats['var_props'] = dict(zip(bins, mags))
            mags, bins = np.histogram(histos['pis'], bins=nice_bins)
            self.stats['pis_props'] = dict(zip(bins, mags))
        io_chunk.close()


    def to_locus(self, trimseq, names, nidxs, snparr, edges):
        """
        write chunk to a .loci format string.
        """
        # store as a list 
        locus = []

        # convert snparrs to snpstrings
        snpstring = "".join([
            "-" if snparr[i, 0] else "*" if snparr[i, 1] else " " 
            for i in range(len(snparr))
        ])

        # get nidx string for getting vcf depths to match SNPs
        if self.data.is_ref:
            # get ref position from nidxs 
            refpos = ":".join(nidxs[0].rsplit(":", 2)[-2:])

            # trim ref position string for edge trims
            chrom, pos = refpos.split(":")
            ostart, end = pos.split("-")
            start = int(ostart) + edges[0]
            end = start + (edges[3] - edges[0])

            # get consens hit indexes and start positions
            nidbits = []
            for bit in nidxs[1:]:
                # handle multiple consens merged
                bkey = []
                for cbit in bit.split(";"):
                    cidx, _, pos = cbit.split(":")

                    # start pos of sample is its consens start pos + ostart
                    # where ostart is the ref position start after trim. So 
                    # how far ahead of ref start does the consens read start.
                    posplus = int(pos.split("-")[0]) - int(ostart)
                    bkey.append("{}:{}".format(cidx, posplus))
                nidbits.append("-".join(bkey))                       

            # put ref back into string and append consens hits
            refpos = "{}:{}-{}".format(chrom, start, end)
            nidbits = [refpos] + nidbits
            nidxstring = ",".join(nidbits)

        # denovo stores start read start position in the nidx string
        else:
            nidxstring = ",".join(nidxs)

        # if not paired data (with an insert)
        for idx, name in enumerate(names):
            locus.append(
                "{}{}".format(
                    self.data.pnames[name],
                    trimseq[idx, :].tostring().decode())
            )
        locus.append("{}{}|{}|".format(
            self.data.snppad, snpstring, nidxstring))
        return "\n".join(locus)


    def filter_minsamp_pops(self, names):
        """
        Returns True if the locus is filtered by min_samples_locus
        """
        minsamp = self.data.params.min_samples_locus
        # default: no population information
        if not self.data.populations:
            if len(names) < minsamp:
                return True
            return False

        # use populations 
        minfilters = []
        for pop in self.data.populations:
            samps = self.data.populations[pop][1]
            minsamp = self.data.populations[pop][0]
            if len(set(samps).intersection(set(names))) < minsamp:
                minfilters.append(pop)
        if any(minfilters):
            return True
        return False
      

    def filter_maxindels(self, trimseq):
        """
        Max size of internal indels. Denovo vs. Ref, single versus paired.
        """
        # get max indels for read1, read2
        inds = maxind_numba(trimseq)
        if inds > self.data.params.max_indels_locus:
            self.filters.loc[self.lidx, "maxind"] = True
            return True
        return False


    def get_snpsarrs(self, trimseq):
        """
        Count nsnps with option to exclude reference sample from count
        """
        # flag = self.data.is_ref and self.data.hackers.exclude_reference
        flag = 0
        snpsarr = np.zeros((trimseq.shape[1], 2), dtype=np.bool_)
        return snpcount_numba(trimseq, snpsarr, int(flag))


    def filter_maxvars(self, trimseq, masked, snpstring):
        """
        Filter on the max allowed variants while masking insert.
        """
        nsnps = snpstring.sum()
        # mask insert area
        if masked is not None:
            maxsnps = masked.shape[1] * self.data.params.max_snps_locus
            if nsnps > maxsnps:
                self.filters.loc[self.lidx, "maxvar"] = True
                return 0, 0, True
            var_prop = snpstring[:, 0].sum() / masked.shape[1]
            pis_prop = snpstring[:, 1].sum() / masked.shape[1]

        # use full locus
        else:
            maxsnps = trimseq.shape[1] * self.data.params.max_snps_locus
            if nsnps > maxsnps:
                self.filters.loc[self.lidx, "maxvar"] = True
                return 0, 0, True
            var_prop = snpstring[:, 0].sum() / trimseq.shape[1]
            pis_prop = snpstring[:, 1].sum() / trimseq.shape[1]
        return var_prop, pis_prop, False


    def filter_maxshared(self, trimseq):
        """
        Paralog detection based on shared heterozygosity across samples.
        """
        #print(trimseq.shape)
        nhs = count_maxhet_numba(trimseq)
        maxhs = self.data.params.max_shared_h_locus * trimseq.shape[0]
        if nhs > maxhs:
            self.filters.loc[self.lidx, "maxshared"] = True
            return True
        return False


def parse_locus(data, chunk):
    """
    grabs the next locus from the chunk file
    """
    names = []
    nidxs = []
    seqs = []

    # advance locus to next, parse names and seqs
    lines = chunk.strip().split("\n")
    for line in lines:
        if line[0] == ">":
            name, nidx = line[1:].rsplit("_", 1)
            names.append(name)
            nidxs.append(nidx)
        else:
            seqs.append(list(line.encode()))

    # filter to include only samples in this assembly (including 'reference')
    mask = [i in data.samples for i in names]

    # select samples in locus to keep
    names = [i for (i, j) in zip(names, mask) if j]

    # if duplicates return True
    if len(set(names)) < len(names):
        return names, nidxs, seqs, True

    # get nidx (refpos/locid string) and sequences of unmasked samples
    nidxs = [i for (i, j) in zip(nidxs, mask) if j]
    seqs = np.array([i for (i, j) in zip(seqs, mask) if j], dtype=np.uint8)
    return names, nidxs, seqs, False


# -------------------------------------------------------------
# jitted subsample func
# -------------------------------------------------------------

@njit
def subsample(snpsmap):
    """
    Subsample snps, one per locus, using snpsmap
    """
    sidxs = np.unique(snpsmap[:, 0])
    subs = np.zeros(sidxs.size, dtype=np.int64)
    idx = 0
    for sidx in sidxs:
        sites = snpsmap[snpsmap[:, 0] == sidx, 1]
        site = np.random.choice(sites)
        subs[idx] = site
        idx += 1
    return subs


# -------------------------------------------------------------
# jitted ChunkProcessor functions (njit = nopython mode)
# -------------------------------------------------------------

@njit
def maxind_numba(block):
    "count the max size of internal indels"
    inds = 0
    for row in range(block.shape[0]):
        where = np.where(block[row] != 45)[0]
        left = np.min(where)
        right = np.max(where)
        obs = np.sum(block[row, left:right] == 45)
        if obs > inds:
            inds = obs
    return inds


@njit
def snpcount_numba(block, snpsarr, rowstart):
    """
    Used to count the number of unique bases in a site for snpstring.
    """
    # iterate over all loci
    for site in range(block.shape[1]):

        # make new array
        catg = np.zeros(4, dtype=np.int16)

        # a list for only catgs
        ncol = block[rowstart:, site]
        for idx in range(ncol.shape[0]):
            if ncol[idx] == 67:    # C
                catg[0] += 1
            elif ncol[idx] == 65:  # A
                catg[1] += 1
            elif ncol[idx] == 84:  # T
                catg[2] += 1
            elif ncol[idx] == 71:  # G
                catg[3] += 1
            elif ncol[idx] == 82:  # R
                catg[1] += 1       # A
                catg[3] += 1       # G
            elif ncol[idx] == 75:  # K
                catg[2] += 1       # T
                catg[3] += 1       # G
            elif ncol[idx] == 83:  # S
                catg[0] += 1       # C
                catg[3] += 1       # G
            elif ncol[idx] == 89:  # Y
                catg[0] += 1       # C
                catg[2] += 1       # T
            elif ncol[idx] == 87:  # W
                catg[1] += 1       # A
                catg[2] += 1       # T
            elif ncol[idx] == 77:  # M
                catg[0] += 1       # C
                catg[1] += 1       # A

        # get second most common site
        catg.sort()

        # if invariant e.g., [0, 0, 0, 9], then nothing (" ")
        if not catg[2]:
            pass
        # store that site is variant as synapomorphy or autapomorphy
        else:           
            if catg[2] > 1:
                snpsarr[site, 1] = True
            else:
                snpsarr[site, 0] = True
    return snpsarr


@njit
def count_maxhet_numba(block):
    """
    Fast jit'd counter of the max number of hetero sites in a locus
    """
    counts = np.zeros(block.shape[1], dtype=np.int16)
    for fidx in range(block.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(block[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts.max()






if __name__ == "__main__":

    pass
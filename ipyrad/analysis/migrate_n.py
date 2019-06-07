#!/usr/bin/env python 

"ipyrad.analysis wrapper for migrate-n"

from __future__ import print_function
from builtins import range

import os
import numpy as np
from ipyrad.assemble.utils import IPyradError


MISSING_IMPORTS = """
To use the ipa.structure module you must install two additional 
libraries which can be done with the following conda command. 

...
"""


class Migrate(object):
    """
    Analysis tool for creating migrate-n input and params files.
    """
    def __init__(
        self, 
        data, 
        name="test",
        workdir="analysis-migrate", 
        imap=None, 
        minmap=None,
        maxloci=None,
        minsnps=0,
        seed=None,
        ):

        # store attributes
        self.name = name
        self.data = os.path.realpath(os.path.expanduser(data))
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.imap = imap
        self.minmap = minmap
        self.seqfile = os.path.join(self.workdir, self.name + ".migrate-seq")
        self.paramfile = os.path.join(self.workdir, self.name + ".migrate-param")
        self.minsnps = minsnps
        self.maxloci = maxloci
        np.random.seed(seed)

        # require imap
        if not self.imap:
            raise IPyradError("imap dictionary required to group samples to pops")

        # check minmap
        if not self.minmap:
            self.minmap = {i: 1 for i in self.imap}
        else:
            if not all([i in self.minmap for i in self.imap]):
                raise IPyradError("keys of minmap and imap do not match.")

        # create workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # check imap pop names
        # if len(i) > 10 for i in         



        
    def write_seqfile(self):
        """
        Write a migrate-n formatted sequence file with samples from imap
        filtered by minmap and limited to maxloci and minsnps. The maxloci
        param randomly samples loci using the argument 'seed' from init.

        Names MUST be 10 characters and MUST NOT start with an int.

        Creates migrate-n format from the .loci file which has already been
        filtered by the .populations minsample args

        npops nloci
        loclen loclen loclen loclen loclen
        ninds ninds ninds ninds A
        ...
        ...
        ninds ninds ninds ninds B
        ...
        ...
        ninds ninds ninds ninds C
        """
        # iterate over loci lines
        indat = iter(open(self.data, 'r'))

        # step 1 load all loci and filter by imap, minmap, and minSNPs
        locdict = {}
        loci = []
        while 1:
            try:
                line = next(indat)
            except StopIteration:
                indat.close()
                break

            # end of locus
            if line.endswith("|\n"):

                # convert to an array
                arrdict = {}                
                names = []
                seqs = []
                for name, seq in locdict.items():
                    names.append(name)
                    seqs.append(list(seq))
                arr = np.array(seqs).astype(bytes).view(np.int8)

                # remove site that are all Ns
                drop = np.all(arr == 78, axis=0)
                arr = arr[:, ~drop]
                
                # count variants
                nvars = np.invert(np.all(arr == arr[0], axis=0)).sum()
                if nvars >= self.minsnps:

                    # check imap, minmap coverage
                    filtered = 0
                    for pop in self.imap:
                        minsamp = self.minmap[pop]
                        pnames = self.imap[pop]
                        nidxs = [names.index(i) for i in pnames if i in pnames]
                        arrdict[pop] = arr[nidxs, :]
                        if len(nidxs) <= minsamp:
                            filtered = 1

                    # passed filtering!
                    if not filtered:
                        loci.append(arrdict)  

                # clear the locus
                locdict = {}

            # just another line of a locus
            else:
                name, seq = line.split()
                locdict[name] = seq

        # step 3 filter by maxloci
        if self.maxloci:
            locidxs = sorted(np.random.choice(range(len(loci)), self.maxloci))
        else:
            locidxs = range(len(loci))

        # step 4 format to file
        with open(self.seqfile, 'w') as out:

            # write header
            nloci = min(len(loci), self.maxloci)
            out.write(
                "{} {} (npops nloci from {}.loci)\n"
                .format(len(self.imap), nloci, self.name)
            )

            # write locus lengths
            dummy = list(self.imap.keys())[0]
            loclens = [loci[i][dummy].shape[1] for i in locidxs]
            out.write(" ".join([str(i) for i in loclens]) + "\n")

            # write loci
            for pop in self.imap:

                # get all loci arrays for this pop
                locs = [loci[i][pop] for i in locidxs]

                # get nsamples in each loc
                nsamps = [loc.shape[0] for loc in locs]
                nsampline = [str(i) for i in nsamps] + [str(pop)]
                out.write(" ".join(nsampline) + "\n")

                # get seqarr for this pop for this locus
                lines = []
                for loc in locs:
                    for sidx, seq in enumerate(loc):
                        seq = b"".join(seq.view("S1")).decode()
                        line = "ind_{:<6}{}".format(sidx, seq)
                        lines.append(line)

                out.write("\n".join(lines) + "\n")



        # ## read in data to sample names
        # loci  = infile.read().strip().split("|")[:-1]
        # for loc in loci:
        #     samps = [i.split()[0].replace(">","") for i in loc.split("\n") if ">" in i]
        #     ## filter for coverage
        #     GG = []
        #     for group,mins in MINS:
        #         GG.append( sum([i in samps for i in taxa[group]]) >= int(mins) )
        #     if all(GG):
        #         keep.append(loc)

        # ## print data to file
        # print >>outfile, len(taxa), len(keep), "( npops nloci from {}.loci".format(self.data.name)

        # ## print all data for each population at a time
        # done = 0
        # for group in taxa:
        #     ## print a list of lengths of each locus
        #     if not done:
        #         loclens = [len(loc.split("\n")[1].split()[-1].replace("x","n").replace("n","")) for loc in keep]
        #         print >>outfile, " ".join(map(str,loclens))
        #         done += 1

        #     ## print a list of number of individuals in each locus
        #     indslist = []
        #     for loc in keep:
        #         samps = [i.split()[0].replace(">","") for i in loc.split("\n") if ">" in i]
        #         inds = sum([i in samps for i in taxa[group]])
        #         indslist.append(inds)
        #     print >>outfile, " ".join(map(str,indslist)), group

        #     ## print sample id, spaces, and sequence data
        #     #for loc in range(len(keep)):
        #     for loc in range(len(keep)):
        #         seqs = [i.split()[-1] for i in keep[loc].split("\n") if \
        #                 i.split()[0].replace(">","") in taxa[group]]
        #         for i in range(len(seqs)):
        #             print >>outfile, group[0:8]+"_"+str(i)+\
        #                   (" "*(10-len(group[0:8]+"_"+str(i))))+seqs[i].replace("x","n").replace("n","")
                
        # outfile.close()


    def array_to_migration_matrix(self, arr):
        rows = ["".join(i) for i in arr]
        print("{" + " ".join(rows) + "}")


    def _write_paramsfile(self):
        pass




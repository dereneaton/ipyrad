#!/usr/bin/env python

"""
Plot the pairwise locus sharing and pairwise missingness.
"""

import tempfile
import h5py
import numpy as np
import pandas as pd
from .locus_extracter import LocusExtracter

# import tested at init
try:
    import toyplot
except ImportError:
    pass
_TOYPLOT_IMPORT = """
This ipyrad analysis tool requires the following software
that you can install with conda using this command:

   conda install toyplot -c conda-forge
"""



class Sharing:
    def __init__(
        self, 
        data,
        imap=None, 
        minmap=None, 
        mincov=0, 
        minsnps=0,
        maxmissing=1.0,
        exclude=None,
        ):
        """
        ...

        Parameters
        ----------
        data: The .seqs.hdf5 database file from ipyrad.
        """
        # check imports
        if not sys.modules.get("toyplot"):
            raise ImportError(_TOYPLOT_IMPORT)

        # load args
        self.data = data
        self.imap = imap
        self.minmap = minmap
        self.mincov = mincov
        self.minsnps = minsnps
        self.maxmissing = maxmissing
        self.exclude = exclude

        # to be filled
        self.phymap = None
        self.slens = None
        self.snames = None
        self.faidict = None

        # run functions
        self._load_phymap()
        self._apply_filters()


    def _load_phymap(self):
        "load the phymap to get information for selecting scaffolds by idx"
        with h5py.File(self.data, 'r') as io5:

            # load the phymap with spans of each RAD tag pre-filtering
            self.phymap = io5["phymap"][:]

            # parse names and lengths from db
            scafnames = [i.decode() for i in io5["scaffold_names"][:]]
            scaflens = io5["scaffold_lengths"][:]

            # organize as a DF
            self.scaffold_table = pd.DataFrame(
                data={
                    "scaffold_name": scafnames,
                    "scaffold_length": scaflens,
                }, 
                columns=["scaffold_name", "scaffold_length"],
            )


            # mask for min scafflen
            # if self.scaffold_minlen:
                # self.scaffold_table = (
                    # self.scaffold_table[self.scaffold_table.scaffold_length > self.scaffold_minlen]
                # )
                # self.scaffold_table.reset_index(inplace=True)
                #self.scaffold_table.reset_index(inplace=True, drop=True)            
            # if self.scaffold_minlen:
            #     self.mask_minlen = np.array(scaflens) > self.scaffold_minlen
            #     scafnames = np.array(scafnames)[self.mask_minlen]
            #     scaflens = np.array(scaflens)[self.mask_minlen]            



    def _apply_filters(self):
        """
        ...
        """

        # filter all loci
        ext = LocusExtracter(
            name="coverage-plot",
            workdir=tempfile.gettempdir(),
            data=self.data,
            imap=self.imap,
            minmap=self.minmap,
            minsnps=self.minsnps,
            maxmissing=self.maxmissing,
            consensus_reduce=self.consensus_reduce,
            mincov=self.mincov,
            exclude=self.exclude,
            quiet=True,
        )
        ext.ipcluster['cores'] = 4
        ext.run(auto=True)

        # subset phymap to the scaffolds in scaffold_idx
        self._phymap = ext.phymap.copy()
        self.phymap = ext.phymap[~ext.phymap.filtered].reset_index()




    def draw(self, scaffold_idxs=None, width=None, height=None, **kwargs):

        # check for .apply_filters()


        # get phymap for scaffs
        scaffilter = np.sum(
            [self.phymap.chroms == i + 1 for i in scaffold_idxs],
            axis=0).astype(bool)
        phymap = self.phymap[scaffilter]

        # report number of loci passing filter on scaff idxs
        print(
            "{} loci passed filtering on selected scaffolds."
            .format(phymap.shape[0]))

        # subset scaff table
        scaffold_table = self.scaffold_table.iloc[scaffold_idxs, :]

        # set up
        width = (width if width else 800)
        height = (height if height else (100 + scaffold_table.shape[0] * 20))
        canvas = toyplot.Canvas(width, height)
        axes = canvas.cartesian(**kwargs)

        # adjusters
        nudge = scaffold_table.scaffold_length.max() * 0.01
        baseline = 0

        # iterate over scaffolds
        for scaff in scaffold_table.index:

            # draw the length of the scaffold
            axes.rectangle(
                0, scaffold_table.scaffold_length[scaff], 
                baseline, baseline - 1,
                color=toyplot.color.Palette()[0],
                opacity=0.7,
            )

            # add markers on the scaffolds
            scafmarkers = phymap[phymap.chroms == scaff + 1]
            axes.scatterplot(
                scafmarkers.pos0,
                np.repeat(baseline - 0.5, scafmarkers.pos0.size),
                color='black',
                marker="|",
            )

            axes.text(
                scaffold_table.scaffold_length[scaff] + nudge,
                baseline - 0.5,
                scaffold_table.scaffold_name[scaff],
                style={"text-anchor": "start"},
                color="black",
            )

            baseline -= 1.5
        axes.y.show = False
        axes.x.label.text = "genomic position (bp)"
        return canvas, axes

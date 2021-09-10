#!/usr/bin/env python


"""
Generates a toyplot drawing of results from abba-baba tests in ipa.
"""

from dataclasses import dataclass, field
import toyplot
import toytree


class Drawing:
    def __init__(self, res, tax, tree, width=500, height=500, sort=False, prune=False, fade=False, zscoreTH=2.5):
        
        self.tests = tax
        self.res = res
        self.ntests = res.shape[0]
        self.zscoreTH = zscoreTH
        self.fade = fade
      
        # if prune tree
        if prune:
            intree = set([])
            for cell in self.tests.values.flatten():
                for tax_ in cell.split(","):
                    intree.add(tax_)
            tree = tree.drop_tips(
                [i for i in tree.get_tip_labels() if i not in intree]
            )
        
        
        # define tree, original tree or prunned tree
        self.tree = tree
        
        
        if sort:
            # split to make cell into a list
            sindex = (
                self.tests
                .applymap(lambda x: x.split(","))
                .applymap(self.tree.get_mrca_idx_from_tip_labels)
                .sort_values(by=["p4", "p3", "p2", "p1"])
            ).index

            # rearrange tables by sindex
            self.tests = self.tests.loc[sindex]
            self.res = self.res.loc[sindex]
            self.tests.reset_index(drop=True, inplace=True)
            self.res.reset_index(drop=True, inplace=True)

        # canvas and axes components
        self.canvas = toyplot.Canvas(width, height)
        self.add_tree_to_canvas()
        self.add_zscores_to_canvas()
        self.add_histos_to_canvas()
        self.add_test_idxs_to_canvas()
        self.add_tip_names_to_canvas()
        self.add_tests_to_canvas()


    def add_tree_to_canvas(self):
        ax0 = self.canvas.cartesian(bounds=("50%", "90%", "5%", "19%"), show=False)
        self.tree.draw(
            axes=ax0, 
            ts='n', 
            layout='d', 
            tip_labels=False, 
            tip_labels_align=True, 
            xbaseline=0.5,
        )
        ax0.rectangle(
            0, self.tree.ntips, 
            0, self.tree.treenode.height, 
            style={"fill": "none"},
        )


    def add_test_idxs_to_canvas(self):
        # test names
        ax4 = self.canvas.cartesian(bounds=("91%", "95%", "21%", "80%"), show=False)
        ax4.rectangle(
            0, 1, 
            0, self.ntests + 1, 
            style={"fill": "none"})
        ax4.text(
            np.repeat(0, self.ntests),
            np.arange(self.ntests) + 1, 
            [str(i) for i in range(self.ntests)][::-1],
            style={"fill": "black", "text-anchor": "start"}
        )


    def add_tip_names_to_canvas(self):
        # tip names
        ax5 = self.canvas.cartesian(bounds=("50%", "90%", "80%", "97%"), show=False)
        ax5.rectangle(0, self.tree.ntips, 0, 1, style={"fill": "none"})
        ax5.text(
            np.arange(self.tree.ntips) + 0.5,
            np.repeat(0.9, self.tree.ntips),
            self.tree.get_tip_labels(),
            angle=-90,
            style={"fill": "black", "text-anchor": "start"},
            annotation=True,
        )


    def add_tests_to_canvas(self):
        # add tests bars to axes
        ax1 = self.canvas.cartesian(
            bounds=("50%", "90%", "21%", "80%"), 
            show=False,
            padding=0,
        )

        # spacer rect
        ax1.rectangle(
            0, self.tree.ntips, 
            0, self.ntests + 1, 
            style={
                "fill": "grey", 
                "fill-opacity": 0.1, 
            },
        )

        # coloring
        COLORS = toyplot.color.Palette()
        colors = [COLORS[0], COLORS[1], toyplot.color.black, COLORS[7]]
        opacities = [1, 1, 1, 1]
        TIPS = self.tree.get_tip_labels()

        # draw blocks
        for idx in range(self.ntests):

            # line tracing
            hidx = self.ntests - idx
            ax1.hlines(hidx, color=toyplot.color.black, style={"stroke-dasharray": "2,4"})
            
            
            #if fade option is true, make half transparent non significant blocks
            if self.fade:
                # check if Z is significant and set opacities for every block
                if self.res.Z[idx] < self.zscoreTH:
                    opacities = [0.6, 0.6, 1, 1] #make both P1 and P2 transparent
                else:
                    if self.res.D[idx] > 0:
                        opacities = [0.6, 1, 1, 1] #make P1 transparent
                    else:
                        opacities = [1, 0.6, 1, 1] #make P2 transparent
        
            

            # get test [name1, name2, name3]
            for cidx, pop in enumerate(["p1", "p2", "p3", "p4"]):
                test = self.tests.iloc[idx][pop]

                
                # get name indices [0, 2, 3]
                tidxs = sorted([TIPS.index(i) for i in test.split(",")])
                                
                # draw blocks connecting index to next until no more.
                blocks = []
                
                # declare a block as [names, initial tip, last tip]
                block = [test.replace(",","\n"), tidxs[0], tidxs[0]]
                for i in range(1, len(tidxs)):
                    if tidxs[i] - tidxs[i - 1] == 1:
                        block[-1] = tidxs[i]
                    else:
                        blocks.append(block)
                        block = [test, tidxs[i], tidxs[i]]

                blocks.append(block)
                blocks[-1][-1] = tidxs[-1]

                
                # draw them (left, right, top, bottom)
                for block in blocks:
                    ax1.rectangle(
                        a=block[1] + 0.25,
                        b=block[2] + 0.75,
                        c=hidx + 0.25, 
                        d=hidx - 0.25,
                        title=block[0],
                        style={
                            "fill": colors[cidx],
                            "stroke": toyplot.color.black,
                            "opacity": opacities[cidx],
                            "stroke-width": 0.5,
                        },
                    )
        ax1.hlines(
            [0, self.ntests + 1], 
            style={"stroke": toyplot.color.black, "stroke-width": 1.5}
        )
        ax1.vlines(
            [0, self.tree.ntips], 
            style={"stroke": toyplot.color.black, "stroke-width": 1.5},
        )        


    def add_zscores_to_canvas(self):
        # add zscores bars to axes
        ax2 = self.canvas.cartesian(
            bounds=("25%", "47%", "21%", "80%"), 
            yshow=False,
            padding=0,
        )

        # the longest bar space
        maxz = max(self.res.Z) + (max(self.res.Z) * .10)

        # spacer rect
        ax2.rectangle(
            -maxz, 0,
            0, self.ntests + 1, 
            style={
                "fill": "grey", 
                "fill-opacity": 0.1,
            },
        )

        # add data bars
        for idx in range(self.ntests):
            hidx = self.ntests - idx
            ax2.hlines(hidx, color='black', style={"stroke-dasharray": "2,4"})
            ax2.rectangle(
                0, -self.res.Z[idx],
                hidx - 0.25, hidx + 0.25, 
                color=toyplot.color.black,
                title="Z-score: " + str(round(-self.res.Z[idx], 2))
            )


        # stylring 
        ax2.x.spine.show = False
        ax2.x.label.text = "Z-score"
        ax2.x.ticks.locator = toyplot.locator.Extended(5, only_inside=True)
        ax2.vlines(
            [ax2.x.domain.min, ax2.x.domain.max, 0, -maxz], 
            style={"stroke": toyplot.color.black, "stroke-width": 1.5},
        )
        ax2.hlines(
            [0, self.ntests + 1],
            style={"stroke": toyplot.color.black, "stroke-width": 1.5},            
        ) 
        
        #zscore threshold
        if -maxz < -self.zscoreTH:
            ax2.vlines(
                -self.zscoreTH, 
                style={
                    "stroke": "grey", 
                    "stroke-dasharray": "2,4", 
                    "stroke-width": 1,
                })


    def add_histos_to_canvas(self):
        # add histograms to axes
        ax3 = self.canvas.cartesian(
            bounds=("5%", "22%", "21%", "80%"), 
            yshow=False, 
            padding=0,
        )

        zmin = min(self.res.D - 3.25 * self.res.bootstd[0])
        zmax = max(self.res.D + 3.25 * self.res.bootstd[0])

        # draw outline and fill
        ax3.rectangle(
            zmin, zmax,
            0, self.ntests + 1, 
            style={
                "fill": "grey", 
                "fill-opacity": 0.1, 
            },
        )

        # iterate over tests to add histos
        for idx in range(self.ntests):
            hidx = self.ntests - idx

            # get fill color
            if self.res.Z[idx] < self.zscoreTH:
                fill = toyplot.color.Palette()[7]
            else:
                if self.res.D[idx] > 0:
                    fill = toyplot.color.Palette()[1]
                else:
                    fill = toyplot.color.Palette()[0]

            # histogram fill
            points = np.linspace(zmin, zmax, 30)
            density = sc.norm.pdf(
                points, loc=self.res.D[idx], scale=self.res.bootstd[idx],
            )
            ax3.fill(
                points, density / density.max() * 0.7,
                baseline=np.repeat(hidx - 0.25, len(points)),
                style={
                    "stroke": 'black', 
                    "stroke-width": 0.5, 
                    "fill": fill},
                title="D-statistic: " + str(round(self.res.D[idx], 2))
            )

        # Z=0 indicator    
        ax3.vlines(
            0, 
            style={
                "stroke": "grey", 
                "stroke-dasharray": "2,4", 
                "stroke-width": 1,
            })

        ax3.vlines(
            [zmin, zmax],
            style={"stroke": "black", "stroke-width": 1.5},
        )
        ax3.hlines(
            [0, self.ntests + 1],
            style={"stroke": "black", "stroke-width": 1.5},
        )        

        # style axes
        ax3.x.label.text = "D-statistic"
        ax3.x.spine.show = False
        ax3.x.ticks.locator = toyplot.locator.Explicit(
            [zmin, 0.0, zmax],
            ["{:.1f}".format(i) for i in [zmin, 0.0, zmax]],
        )

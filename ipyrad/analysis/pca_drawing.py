#!/usr/bin/env python

"""
Generates PCA plots. This class is not for general use, and is used
internally from the ipa.pca tool in its .draw() function to return
a toyplot drawing.
"""

from typing import Optional, Dict, List
import itertools
from dataclasses import dataclass, field
from copy import deepcopy

import numpy as np
import toyplot
import toyplot.svg
import toyplot.pdf
try:
    from sklearn.linear_model import LinearRegression
    from sklearn.neighbors import NearestCentroid
except ImportError:
    pass
from ipyrad.assemble.utils import IPyradError


@dataclass
class Drawing:
    tool: 'ipa.pca'
    ax0: int=0
    ax1: int=1
    width: int=300
    height: int=300
    cycle: int=8
    colors: List[str]=None
    opacity: Optional[float]=None
    shapes: List[str]=None
    size: int=10
    legend: bool=True
    legend_width: int=100
    label: Optional[str]=None
    outfile: Optional[str]=None
    imap: Optional[Dict[str,List[str]]]=None
    axes: Optional['toyplot.coordinates.Cartesian']=None

    # non param attributes to be filled
    nreplicates: int = field(default=None, init=False)
    mean_variances: np.ndarray = field(default=None, init=False)
    canvas: 'toyplot.Canvas' = field(default=None, init=False)
    imap_rev: Dict[str,str] = field(default=None, init=False)
    pstyles: Dict = field(default_factory=dict, init=False)
    rstyles: Dict = field(default_factory=dict, init=False)

    # non-param attributes
    def __post_init__(self):
        """Check some parameter entries."""
        # use a copy of the tool
        self.tool = deepcopy(self.tool)

        # require an imap
        self.imap = self.imap if self.imap else self.tool.imap

        # require run to have been called.
        if self.tool._variances is None:
            raise IPyradError(
                "You must first call run() before calling draw().")

        # check that requested axes exist
        max_axis = max(self.ax0, self.ax1)
        naxes = self.tool._variances[0].size
        assert max_axis < naxes, (
            f"Cannot draw axis {max_axis}, data set only has {naxes} axes.")

    def run(self):
        """Run subfunctions to generate drawing."""
        self._get_mean_run_variances()
        self._get_regressed_run_loadings()
        self._get_canvas_and_axes()
        self._get_marker_styles()
        marks = self._draw_markers()
        self._add_legend()
        self._save_figure()
        return self.canvas, self.axes, marks

    def _get_mean_run_variances(self):
        """Get mean mean variance on each axis across replicate runs."""
        self.nreplicates = len(self.tool._variances)
        self.mean_variances = np.array(
            list(self.tool._variances.values())
        ).mean(axis=0)

    def _get_regressed_run_loadings(self):
        """Runs a CLUMPP-like test to align replicate runs.

        Replicate runs may have their axes flipped randomly so we run
        a simple linear regression on the loadings for each axis for
        replicate runs to flip axes so they align the same as rep 0.
        """
        model = LinearRegression()
        for idx in range(1, self.nreplicates):
            for axis in [self.ax0, self.ax1]:
                orig = self.tool._loadings[0][:, axis].reshape(-1, 1)
                new = self.tool._loadings[idx][:, axis].reshape(-1, 1)
                swap = (self.tool._loadings[idx][:, axis] * -1).reshape(-1, 1)

                # get r^2 for both model fits
                model.fit(orig, new)
                coeff0 = model.coef_[0][0]
                model.fit(orig, swap)
                coeff1 = model.coef_[0][0]

                # if swapped fit is better make this the data
                if coeff1 > coeff0:
                    self.tool._loadings[idx][:, axis] = (
                        self.tool._loadings[idx][:, axis] * -1)

    def _get_canvas_and_axes(self):
        """Setup and style the Canvas and Cartesian axes styles.

        The axes can be further styled by the user after the plot is
        returned if they want, default style here turns tick marks on,
        increases font size, and adds labels that show the variance
        explained for PC models.

        If a legend is requested then the requested width of the
        Canvas will be increased by 100px on the right side.
        """
        # get axis labels for PCA or TSNE plot
        if self.mean_variances[self.ax0] >= 0.0:
            exp0 = self.mean_variances[self.ax0] * 100
            exp1 = self.mean_variances[self.ax1] * 100
            xlab = f"PC{self.ax0} ({exp0:.1f}% explained)"
            ylab = f"PC{self.ax1} ({exp1:.1f}% explained)"
        else:
            xlab = f"{self.tool._model} component 1"
            ylab = f"{self.tool._model} component 2"

        # add room for legend
        margin_right = 60
        if self.legend:
            self.width += self.legend_width
            margin_right += self.legend_width

        # create new Canvas and axes with room for a legend.
        if not self.axes:
            self.canvas = toyplot.Canvas(self.width, self.height)
            self.axes = self.canvas.cartesian(
                xlabel=xlab,
                ylabel=ylab,
                bounds=(60, -margin_right, 60, -60),
                padding=20,
            )
        else:
            self.axes.x.label.text = xlab
            self.axes.y.label.text = ylab

        # style axes
        self.axes.x.spine.style["stroke-width"] = 1.5
        self.axes.y.spine.style["stroke-width"] = 1.5
        self.axes.x.ticks.labels.style["font-size"] = "12px"
        self.axes.y.ticks.labels.style["font-size"] = "12px"
        self.axes.x.label.style['font-size'] = "14px"
        self.axes.y.label.style['font-size'] = "14px"
        self.axes.x.ticks.show = True
        self.axes.y.ticks.show = True
        self.axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
        self.axes.y.ticks.locator = toyplot.locator.Extended(only_inside=True)

        if self.label:
            self.axes.label.text = self.label
            self.axes.label.style['font-size'] = "18px"
            self.axes.label.offset = 25

    def _get_marker_styles(self):
        """Build marker styles (color, size, shape).

        Builds for individual or replicate (low opacity) markers,
        and able to cycle over few or many categories of IMAP.
        """
        # make reverse imap dictionary
        self.imap_rev = {}
        for pop, vals in self.imap.items():
            for val in vals:
                self.imap_rev[val] = pop

        # the max number of pops until color cycle repeats
        # If the passed in number of colors is big enough to cover
        # the number of pops then set cycle to len(colors)
        # If colors == None this first `if` falls through (lazy evaluation)
        if (self.colors is not None) and len(self.colors) >= len(self.imap):
            self.cycle = len(self.colors)
        else:
            self.cycle = min(self.cycle, len(self.imap))

        # get color list repeating in cycles of cycle
        if not self.colors:
            self.colors = itertools.cycle(
                toyplot.color.broadcast(
                    toyplot.color.brewer.map("Spectral"), shape=self.cycle,
                )
            )
        else:
            self.colors = itertools.cycle(self.colors)
            # assert len(colors) == len(imap), "len colors must match len imap"

        # get shapes list repeating in cycles of cycle up to 5 * cycle
        if not self.shapes:
            self.shapes = itertools.cycle(np.concatenate([
                np.tile("o", self.cycle),
                np.tile("s", self.cycle),
                np.tile("^", self.cycle),
                np.tile("d", self.cycle),
                np.tile("v", self.cycle),
                np.tile("<", self.cycle),
                np.tile("x", self.cycle),
            ]))
        else:
            self.shapes = itertools.cycle(self.shapes)
        # else:
            # assert len(shapes) == len(imap), "len colors must match len imap"

        # assign styles to populations and to legend markers (no replicates)
        for pop in self.imap:
            icolor = next(self.colors)
            ishape = next(self.shapes)

            try:
                color = toyplot.color.to_css(icolor)
            except Exception:
                color = icolor

            self.pstyles[pop] = toyplot.marker.create(
                size=self.size,
                shape=ishape,
                mstyle={
                    "fill": color,
                    "stroke": "#262626",
                    "stroke-opacity": 1.0,
                    "stroke-width": 1.5,
                    "fill-opacity": (self.opacity if self.opacity else 0.75),
                },
            )

            self.rstyles[pop] = toyplot.marker.create(
                size=self.size,
                shape=ishape,
                mstyle={
                    "fill": color,
                    "stroke": "none",
                    "fill-opacity": (
                        self.opacity / self.nreplicates if self.opacity
                        else 0.9 / self.nreplicates
                    ),
                },
            )

    def _draw_markers(self):
        """Assigns pstyles to samples and draws as scatterplot marks.
        """
        # assign styled markers to data points
        pmarks = []
        rmarks = []
        for name in self.tool.names:
            pop = self.imap_rev[name]
            pmarks.append(self.pstyles[pop])
            rmarks.append(self.rstyles[pop])

        # if not replicates then just plot the points
        if self.nreplicates < 2:
            mark = self.axes.scatterplot(
                self.tool._loadings[0][:, self.ax0],
                self.tool._loadings[0][:, self.ax1],
                marker=pmarks,
                title=self.tool.names,
            )
            return mark

        # add the replicates cloud points
        for rep in range(self.nreplicates):

            # get transformed coordinates and variances
            _ = self.axes.scatterplot(
                self.tool._loadings[rep][:, self.ax0],
                self.tool._loadings[rep][:, self.ax1],
                marker=rmarks,
            )

        # compute centroids
        xarr = np.concatenate([
            np.array([
                self.tool._loadings[idx][:, self.ax0], 
                self.tool._loadings[idx][:, self.ax1],
            ]).T
            for idx in range(self.nreplicates)
        ])
        yarr = np.tile(np.arange(len(self.tool.names)), self.nreplicates)
        clf = NearestCentroid()
        clf.fit(xarr, yarr)

        # draw centroids
        mark = self.axes.scatterplot(
            clf.centroids_[:, 0],
            clf.centroids_[:, 1],
            title=self.tool.names,
            marker=pmarks,
        )
        return mark

    def _add_legend(self):
        """Optionally adds a legend to the right side of figure.

        Getting this to fit nicely with long names can be tricky, 
        which is why we added the legend_width option. However, it
        may be useful to play with more tweaking in the future.
        """
        if self.legend and (self.canvas is not None):
            self.canvas.legend(
                list(self.pstyles.items()),
                bounds=(-self.legend_width - 60, -60, 60, -60)
            )

    def _save_figure(self):
        """Write figure to pdf/svg."""
        if self.outfile and (self.canvas is not None):
            if self.outfile.endswith(".pdf"):
                toyplot.pdf.render(self.canvas, self.outfile)
            elif self.outfile.endswith(".svg"):
                toyplot.svg.render(self.canvas, self.outfile)
            else:
                raise IPyradError(
                    "outfile arg must end with .svg or .pdf to select "
                    "a supported format.")

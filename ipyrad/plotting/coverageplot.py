#!/usr/bin/env python2

""" plot coverage depths for samples """

from __future__ import print_function
import toyplot


def depthgram(data, samples=[], dims=(0,0), type="total"):
	""" plots histogram of coverages across clusters"""

	## get canvas dimensions based on n-samples
	if any(dims):
		## user-supplied dimensions (...)
		print(dims)
	else:
		if sum(dims) <= 4:
			## set dimension to N samples 
			dims[1] = nsamples
		else:
			dims[0] = 4
			dims[1] = nsamples/4

	## create canvas
	canvas = toyplot.Canvas(width=250*dims[0], height=250*dims[1])

	## select samples to be plotted
	if samples:
		subsamples = [i for i in data.samples if i.name in samples]
	else:
		subsamples = data.samples

	## fill in each panel of canvas with a sample
	for panel, sample in enumerate(subsamples):
		axes = canvas.axes(grid=(dims[0], dims[1], panel))
		mark = axes.bars(np.histogram(sample, 20))


# if __name__ == "__main__":
# 	import ipyrad as ip
# 	import numpy as np
# 	TEST = ip.Sample()
# 	TEST.depths["total"] = np.random.normal(20, 3, 2000)

# 	#DATA = ip.Assemble()
# 	#DATA.link_samples(TEST)

# 	coverageplot(TEST)



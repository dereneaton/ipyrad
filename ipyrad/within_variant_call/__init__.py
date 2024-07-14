#!/usr/bin/env python

"""Step 4: variant calls within samples.

Substeps:
1. denovo: infer H,E from clusters.
2. call site and allele variants.
3. store haploid or diploid consensus alleles.
4. apply filters to exclude bad clusters.
5. apply trimming to clean edges of good clusters.
6. write stats for H-hat, E-hat, H, filters, trimming.
"""

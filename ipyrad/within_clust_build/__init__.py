#!/usr/bin/env python

"""Step 3 of denovo assembly: build clusters and align.

Substeps:
1. Cluster are built from (_derep, _matches) files.
2. Muscle alignment is run on each cluster.
3. Update Sample objects and save json.

"""

#!/usr/bin/env python

try:
    from tetrad import Tetrad as external_tetrad
    Tetrad = external_tetrad
except ImportError:
    pass


class TetradMessage(object):
    """
    You must first install tetrad: 'conda install tetrad -c eaton-lab'
    """
    def __init__(self,
        name=None, 
        data=None,
        workdir="analysis-tetrad",
        nquartets=0, 
        nboots=0,
        load=False,
        save_invariants=False,
        *args, 
        **kwargs):    
        raise ImportError(MESSAGE)


Tetrad = TetradMessage


MESSAGE = """
To run tetrad you must first intall the tetrad package. 
Please run: 

  conda install tetrad -c eaton-lab
"""

#!/usr/bin/env python

class IPyradError(Exception):
    """Raise a custom exception that will report with traceback.

    This is used to catch and report internal errors in the code, 
    and the traceback will include the source error and error type
    for debugging.
    """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

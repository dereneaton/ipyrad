#!/usr/bin/env python2

""" extra class objects """

from __future__ import print_function
import multiprocessing

 
class Worker(multiprocessing.Process):
    """ multiprocessing object """
 
    def __init__(self, work_queue, result_queue, func):
 
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.func = func
 
    def run(self):
        while not self.kill_received:
            # get a task
            if self.work_queue.empty():
                break
            else:
                #job = self.work_queue.get_nowait()
                job = self.work_queue.get()
 
            # the actual processing
            res = self.func(*job)
            # store the result
            self.result_queue.put(res)

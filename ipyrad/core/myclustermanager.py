

import os
from traitlets.config.configurable import LoggingConfigurable
from traitlets import Dict, Instance, Float

from ipyparallel.apps.ipclusterapp import IPClusterStart

from ipyparallel.apps.ipclusterapp import ipclusterapp



klass = ipclusterapp.find_launcher_class("MPI")

clusterapp = ipclusterapp.IPClusterApp()


class MyClusterManager(LoggingConfigurable):
    """ Launch ipcluster if not already started for this pid """

    ## we're gonna use the u'default' profiles
    profiles = Dict()

    ## traitlets objects.
    delay = Float(1., config=True, help="delay (in s) between starting"\
                                      +"the controller and the engines")
    loop = Instance('zmq.eventloop.ioloop.IOLoop')


    def build_launcher(self, clsname, kind=None):
        """import and instantiate a Launcher based on importstring"""
        try:
            klass = find_launcher_class(clsname, kind)
        except (ImportError, KeyError):
            self.log.fatal("Could not import launcher class: %r"%clsname)
            self.exit(1)

        launcher = klass(
            work_dir=u'.', parent=self, log=self.log,
            profile_dir=self.profile_dir.location, cluster_id=self.cluster_id,
        )
        return launcher


    def build_launchers(self, profile_dir, cluster_id, controller, ):
        #from ipyparallel.apps.ipclusterapp import IPClusterStart
        try:
            c_klass = clusterapp.find_launcher_class(controller, "Controller")
            e_klass = clusterapp.find_launcher_class(controller, "EngineSet")
        except (ImportError, KeyError):
            LOGGER.fatal("could not launch type %s", controller)

        launcher = klass(
            work_dir=u'.', parent=self, log=self.log,
            profile_dir=self.profile_dir.location, cluster_id=self.cluster_id)


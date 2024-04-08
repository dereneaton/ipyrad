#!/usr/bin/env python

"""A function to report messages from ipyclient to the logger.

Call `track_remote_jobs` to track jobs running on an ipyclient and
report messages to the logger. Messages are sent to the logger from
remote functions by calling `print("@@LOGLEVEL: message")`.

Example Usage
-------------
>>> rasyncs = {}
>>> for sname, sample in data.samples.items():
>>>     rasyncs[sname] = lbview.apply(func, sample)
>>> results = track_remote_jobs(rasyncs, self.ipyclient)

"""

from loguru import logger

logger = logger.bind(name="ipyrad")


def progress(remote_messages: str) -> None:
    for msg in remote_messages.split("@@")[1:]:
        log_level, log_msg = msg.split(":", 1)
        logger.log(log_level.strip(), log_msg.strip())


def track_remote_jobs(rasyncs, ipyclient):
    # add a callback to log stdout from engine when a job finishes
    for rasync in rasyncs.values():
        rasync.add_done_callback(lambda x: progress(x.stdout))

    # wait, catch results and/or interrupts, and return results.
    results = {}
    try:
        ipyclient.wait()
        for job, rasync in rasyncs.items():
            results[job] = rasync.result()

    except KeyboardInterrupt:
        logger.error("KeyboardInterrupt by user.")
        ipyclient.cluster.signal_engines_sync(signum=2)
        raise KeyboardInterrupt("KeyboardInterrupt by user.")
    return results


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    with ip.Cluster(cores=4) as ipyclient:
        lbview = ipyclient.load_balanced_view()

        # jobs report to DEBUG
        rasyncs = {}
        for i in range(4):
            rasyncs[i] = lbview.apply(lambda x: print(f"@@DEBUG: {x}"), i)
        track_remote_jobs(rasyncs, ipyclient)

        # jobs report to WARNING
        rasyncs = {}
        for i in range(4):
            rasyncs[i] = lbview.apply(lambda x: print(f"@@WARNING: {x}"), i)
        track_remote_jobs(rasyncs, ipyclient)

        # jobs report an ERROR. Code must still raise the exception itself
        # by examining the rasyncs.
        rasyncs = {}
        for i in range(4):
            rasyncs[i] = lbview.apply(lambda x: print(f"@@ERROR: {x}"), i)
        track_remote_jobs(rasyncs, ipyclient)

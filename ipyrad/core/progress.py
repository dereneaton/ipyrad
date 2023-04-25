#!/usr/bin/env python

from loguru import logger

logger = logger.bind(name="ipyrad")


def progress(remote_messages: str) -> None:
    for msg in remote_messages.split("@@")[1:]:
        log_level, log_msg = msg.split(":", 1)
        logger.log(log_level.strip(), log_msg.strip())


# not yet used.
def track_remote_jobs(rasyncs, ipyclient):

    # add a callback to log progress from stdout on engines
    for rasync in rasyncs.values():
        rasync.add_done_callback(lambda x: progress(x.stdout))

    # wait, catch results and/or interrupts, and return results.
    results = {}
    try:
        ipyclient.wait()
        for job, rasync in rasyncs.items():
            results[job] = rasync.result()

    except KeyboardInterrupt:
        ipyclient.cluster.signal_engines_sync(signum=2)
        raise KeyboardInterrupt("KeyboardInterrupt by user.")
    return results

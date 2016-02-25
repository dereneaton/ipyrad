#!/usr/bin/env python2.7

import re
import sys
import fileinput
import subprocess

# Fetch version from git tags, and write to version.py.
# Also, when git is not available (PyPi package), use stored version.py.
initfile = "ipyrad/__init__.py"

version_git = sys.argv[1]
print("Setting new version to - {}".format(version_git))

# Write version to ipyrad/__init__.py
# Also set default __loglevel__ to ERROR so we don't check in
# DEBUG by accident.
for line in fileinput.input(initfile, inplace=1):
    if line.strip().startswith("__version__"):
        line = "__version__ = \""+version_git+"\""

    if line.strip().startswith("__loglevel__"):
        line = "__loglevel__ = \"ERROR\""

    print(line.strip("\n"))

try:
    subprocess.call(["git", "add", initfile])
    subprocess.call(["git", "commit", "-m \"Updating ipyrad/__init__.py to "+\
                    "version - {}".format(version_git)])
    subprocess.call(["git", "push"])
    subprocess.call(["git", "tag", "-a", version_git, "-m", "Updating ipyrad to "+\
                    "version - {}".format(version_git)])
    subprocess.call(["git", "push", "origin", version_git])
except Exception as e:
    print("Something broke - {}".format(e))

print("Push new version of conda installer")

try:
    subprocess.call(["conda", "build", "conda.recipe"])
except Exception as e:
    print("something broke - {}".format(e))

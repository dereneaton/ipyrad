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
for line in fileinput.input(initfile, inplace=1):
        if "__version__" in line:
            line = "__version__ = \""+version_git+"\""
        print(line.strip("\n"))
        elif: "__loglevel__" in line:
            line = "__loglevel__ = \"ERROR\""
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

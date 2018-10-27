#!/usr/bin/env python

from __future__ import print_function
import re
import sys
import time
import fileinput
import subprocess

# Fetch version from git tags, and write to version.py.
# Also, when git is not available (PyPi package), use stored version.py.
initfile = "ipyrad/__init__.py"

version_git = sys.argv[1]
print("Setting new version to - {}".format(version_git))

# Make sure we you have the current changes in master w/o any
# merge conflicts
try:
    subprocess.check_call(["git", "pull"])
except Exception as inst:
    sys.exit("You have unmerged conflicts from master. Fix these then "
             + "rerun versioner.py\n{}".format(inst))

# Make sure to track previous release version
prev_version = ""
with open(initfile, 'r') as infile:
    for line in infile:
        if line.strip().startswith("__version__"):
            prev_version = str(line.split("=")[1].strip().strip("\""))
print("last version - {}".format(prev_version))

# Update docs/releasenotes.rst to include commit history for this version
release_file = "docs/releasenotes.rst"

# Get all commits since last tag
# git log --pretty=oneline tagA..
cmd = "git log --pretty=oneline {}..".format(prev_version)
new_commits = subprocess.check_output(cmd.split())

# Split off just the first element, cuz we don't care about the
# commit tag
commit_lines = [x.split(" ", 1) for x in new_commits.split("\n")]

checkfor = "Merge branch 'master' of https://github.com/dereneaton/ipyrad"
# Write updates to releasenotes.rst
for line in fileinput.input(release_file, inplace=1):
    if line.strip().startswith("=========="):
        line += "\n"
        line += version_git + "\n"
        line += "-" * len(version_git) + "\n"
        line += "**" + time.strftime("%B %d, %Y") + "**\n\n"
        for commit in commit_lines:
            try:
                # Squash merge commits from the releasenotes cuz it annoying
                # Also any cosmetic commits
                if commit[1] == checkfor:
                    continue
                if "cosmetic" in commit[1]:
                    continue
                line += "- " + commit[1] + "\n"
            except:
                pass
    print(line.strip("\n"))


# Write version to ipyrad/__init__.py
for line in fileinput.input(initfile, inplace=1):
    if line.strip().startswith("__version__"):
        line = "__version__ = \""+version_git+"\""
    print(line.strip("\n"))


try:
    subprocess.call(["git", "add", release_file])
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
    subprocess.call(["conda", "build", "conda.recipe/ipyrad"])
except Exception as e:
    print("something broke - {}".format(e))

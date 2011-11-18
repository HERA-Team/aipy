#! /usr/bin/env python
"""
Helps you figure out which version of AIPY you have installed.
Lists:
    version
    install location
    last commit
    branch
"""
import aipy
print "AIPY base version in use:",
try:
    print aipy.__version__
except(AttributeError):
    print "Version number not Found"
print "Install location:",aipy.__file__
print "Last Git commit log:"
try:
    print aipy.__gitlog__
except(AttributeError):
    print "     Gitlog not found"
print "Branch:",
try:
    print aipy.__branch__
except(AttributeError):
    print "Branch Unknown"


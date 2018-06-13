#!/usr/bin/env python

"""
Helps you figure out which version of AIPY you have installed.
Lists:
    version
    install location
    last commit
    branch
"""

from __future__ import print_function, division, absolute_import

import aipy
print("AIPY base version in use: ", end='')
try:
    print(aipy.__version__)
except(AttributeError):
    print("Version number not Found")
print("Install location:", aipy.__file__)
print("Last Git commit log:")
try:
    print(aipy.__gitlog__)
except(AttributeError):
    print("not found")
print("Branch: ", end='')
try:
    print(aipy.__branch__)
except(AttributeError):
    print("unknown")

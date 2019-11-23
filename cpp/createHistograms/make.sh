#!/bin/bash

g++ `root-config --cflags --glibs` -g createHistograms.cxx makeHistograms.cxx -o createHistograms
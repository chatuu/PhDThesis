#!/bin/bash

g++ `root-config --cflags --glibs` -g createHistograms.cxx makeHistograms.cxx -o createHistograms

g++ `root-config --cflags --glibs` -g tracksVsProngs.cxx makeHistograms.cxx -o tracksVsProngs

#!/bin/bash
cafe -b -q --numuccinc -l $1 createSpectra.C+ > ./logFiles/createSpectra_log_$(date +\%Y-\%m-\%d_\%H\%M).log

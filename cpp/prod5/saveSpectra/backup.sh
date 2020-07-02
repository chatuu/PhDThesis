#!/bin/bash

tar -zcf Backup_SaveSpectra_$(date +\%Y-\%m-\%d_\%H\%M).tar --exclude='*.root' --exclude='*.png' --exclude='*.out' --exclude='*.tar' --exclude-vcs --exclude='*.o' --exclude='*.d' --exclude='*.pcm' --exclude='*.so' --exclude='*.swp' --exclude='gridJob/*' --exclude='tmp/*' --exclude='*debug' --exclude='*.tar.bz2' *

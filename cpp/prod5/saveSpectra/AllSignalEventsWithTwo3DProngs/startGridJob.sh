#!/bin/bash

declare -a fileList=("createSpectra.C"
	"vars.h"
	"headers.h"
	"structs.h"
	"cuts.h")
echo -e "\n"
echo "copying the files:"

for i in "${fileList[@]}"; do
	echo "$i"
	cp $i ../gridJob2
done

echo -e "Accessing gridJob2 folder:\n"
cd ../gridJob2/
ls -lrt
tar -zc -f XSec_Testing.tar --exclude='*.root' --exclude='*.png' --exclude='*.out' --exclude='*.tar' --exclude-vcs --exclude='*.o' --exclude='*.d' --exclude='tmp/*' --exclude='*debug' --exclude='*.tar.bz2' *

echo "done..!"

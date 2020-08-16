#!/bin/bash

location=$(pwd)
echo "#!/bin/bash
submit_cafana.py -n \$1 -ss -r S20-05-04 -o /pnfs/nova/scratch/users/ckuruppu --user_tarball XSec_Testing.tar ${location}/createSpectra.C \$2" >submitJob.sh

declare -a fileList=("createSpectra.C"
	"vars.h"
	"headers.h"
	"structs.h"
	"cuts.h"
	"submitJob.sh")
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

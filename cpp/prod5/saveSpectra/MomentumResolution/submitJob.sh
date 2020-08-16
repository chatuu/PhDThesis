#!/bin/bash
submit_cafana.py -n $1 -ss -r S20-05-04 -o /pnfs/nova/scratch/users/ckuruppu --user_tarball XSec_Testing.tar /nova/app/users/ckuruppu/workingThesis/cpp/prod5/saveSpectra/MomentumResolution/createSpectra.C $2

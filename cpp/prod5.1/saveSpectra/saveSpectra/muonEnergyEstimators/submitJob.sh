#!/bin/bash
submit_cafana.py -n $1 -ss -r development -o /pnfs/nova/scratch/users/ckuruppu --user_tarball XSec_Testing.tar /nova/app/users/ckuruppu/workingThesis/cpp/prod5.1/saveSpectra/saveSpectra/muonEnergyEstimators/createSpectra.C $2

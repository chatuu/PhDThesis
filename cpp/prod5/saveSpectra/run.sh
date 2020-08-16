#!/bin/bash
cafe -l $1 -b -q createSpectraInt.C+
cp createSpectraInt.C createSpectra.C
sed -i 's/SpectrumLoader lNDMC("def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");/\/\/SpectrumLoader lNDMC("def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");/g' createSpectra.C

sed -i 's/\/\/SpectrumLoader lNDMC("defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");/SpectrumLoader lNDMC("defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");/g' createSpectra.C



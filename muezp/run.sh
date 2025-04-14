#!/bin/bash

if [ ! "$#" = 2 ]; then
    >&2 echo "usage: $(basename "$0") <incoming_energy> <id>"
    exit 1
fi
INCOMING_ENERGY="$1"
ID=$["$2"]
MUON_ENERGY_ID=$[${ID} / 100]
ZP_MASS_ID=$[${ID} % 100]
MUON_ENERGY="$(head -$[${MUON_ENERGY_ID}+1] mup_energy_"${INCOMING_ENERGY}"GeV_pb_median.txt | tail -1)"
ZP_MASS="$(head -$[${ZP_MASS_ID}+1] mup_zp_mass_"${INCOMING_ENERGY}"GeV_pb_median.txt | tail -1)"

cat run.dat \
    | sed 's@\$SUB:workdir@/tmp/MueZp_'"${MUON_ENERGY}"'GeV_Zp_'"${ZP_MASS}"'GeV@g' \
    | sed 's@\$SUB:zp_mass@'"${ZP_MASS}"'@g' \
    | sed 's@\$SUB:nevent@100000@g' \
    | sed 's@\$SUB:muon_energy@'"${MUON_ENERGY}"'@g' \
    | sed 's@\$SUB:electron_energy@0.511e-3@g' \
    | /publicfs/cms/user/gaoleyun/MG5_aMC_v3_6_2/bin/mg5_aMC

mkdir -p mg5root
cp '/tmp/MueZp_'"${MUON_ENERGY}"'GeV_Zp_'"${ZP_MASS}"'GeV'/Events/run_01/*.root mg5root/'MueZp_'"${MUON_ENERGY}"'GeV_Zp_'"${ZP_MASS}"'GeV'.root

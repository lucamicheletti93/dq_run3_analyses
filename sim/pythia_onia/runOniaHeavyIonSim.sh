NJOBS=20
NCORES=5
NEVENTSPERJOB=1000
MODE=kHeavyIon
PROCESS=kSoftQCD
SYSTEM=kPbPb
SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))
#parallel -j $NCORES "root -l -q -b 'simulateOnia.cc($NEVENTSPERJOB, $MODE, $PROCESS, false, 5360, $SYSTEM, {}, \"sim_pythia_onia_OO_5.36TeV_LHC25i4/pythia8_onia_${MODE}_${PROCESS}_seed{}.root\")' > log_onia_${MODE}_seed{}.txt" ::: $(seq $SEEDSTART $SEEDEND)
parallel -j $NCORES "root -l -q -b 'simulateOnia.cc($NEVENTSPERJOB, $MODE, $PROCESS, false, 5360, $SYSTEM, {}, \"sim_pythia_HI_PbPb_5.36TeV/pythia8_${MODE}_${PROCESS}_seed{}.root\")' > log_onia_${MODE}_seed{}.txt" ::: $(seq $SEEDSTART $SEEDEND)
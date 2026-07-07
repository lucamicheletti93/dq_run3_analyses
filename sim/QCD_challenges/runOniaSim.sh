NJOBS=50
NCORES=5
NEVENTSPERJOB=500000
MODE=kMonash
PROCESS=kSoftQCD
SYSTEM=kPP
SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))
parallel -j $NCORES "root -l -q -b 'simulateOnia.cc($NEVENTSPERJOB, $MODE, $PROCESS, false, 5360, $SYSTEM, {}, \"sim_pythia_onia_5.36TeV/pythia8_onia_${MODE}_${PROCESS}_seed{}.root\")' > log_onia_${MODE}_seed{}.txt" ::: $(seq $SEEDSTART $SEEDEND)
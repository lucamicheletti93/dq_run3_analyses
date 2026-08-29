NJOBS=100
NCORES=20
NEVENTSPERJOB=1000000
MODE=kMonash
PROCESS=kSoftQCD
SYSTEM=kPP
SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))
parallel -j $NCORES "root -l -q -b 'simulateOnia.cc($NEVENTSPERJOB, $MODE, $PROCESS, false, 13600, $SYSTEM, {}, \"sim_pythia_onia_13.6TeV/pythia8_onia_${MODE}_${PROCESS}_seed{}.root\")' > log_onia_${MODE}_seed{}.txt" ::: $(seq $SEEDSTART $SEEDEND)
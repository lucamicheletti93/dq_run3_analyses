NJOBS=50
NCORES=5
NEVENTSPERJOB=500000
MODE=kMonash
PROCESS=kSoftQCD
SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))
parallel -j $NCORES "root -l -q -b 'simulateOnia.cc($NEVENTSPERJOB, $MODE, $PROCESS, false, 13600, {}, \"output/pythia8_onia_${MODE}_${PROCESS}_seed{}.root\")' > log_onia_${MODE}_seed{}.txt" ::: $(seq $SEEDSTART $SEEDEND)
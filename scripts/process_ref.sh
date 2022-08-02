#!/bin/sh

REF1="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/acq$1"
REF="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/ref$1" 

./Judith -c convert -i $REF1.bin -o $REF1.root -n 5000
### Alignment ###
./Judith -c coarseAlign -i $REF.root -r configs/sps-kartel.cfg -t configs/sps-tb2.cfg -n 5000
#./Judith -c fineAlign -i $REF.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg
### end ##
./Judith -c process -i $REF.root -o $REF-p.root -r configs/sps-kartel.cfg -t configs/sps-tb2.cfg -R $REF-results.root -n 5000


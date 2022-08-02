#!/bin/sh

RAW="/afs/cern.ch/work/b/bojan/Judith/data/may2018/raw/acq$1"
CONV="/afs/cern.ch/work/b/bojan/Judith/data/may2018/conv/ref$1" 
ANA="/afs/cern.ch/work/b/bojan/Judith/data/may2018/ana/ref$1" 

./Judith -c convert -i $RAW.bin -o $CONV.root -n 3000000
### Alignment ###
#./Judith -c coarseAlign -i $CONV.root -r configs/sps-kartel-june2017.cfg -t configs/sps-tb-may2018.cfg
#./Judith -c fineAlign -i $CONV.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg
### end ##
#./Judith -c process -i $CONV.root -o $ANA-p.root -r configs/sps-kartel-june2017.cfg -t configs/sps-tb-may2018.cfg -R $ANA-results.root


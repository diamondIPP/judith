#!/bin/sh

RAW="/afs/cern.ch/work/b/bojan/Judith/data/oct2017/raw/acq$1"
CONV="/afs/cern.ch/work/b/bojan/Judith/data/oct2017/conv/ref$1" 
ANA="/afs/cern.ch/work/b/bojan/Judith/data/oct2017/ana/ref$1" 

./Judith -c convert -i $RAW.bin -o $CONV.root -n 900000
### Alignment ###
#./Judith -c coarseAlign -i $CONV.root -r configs/sps-kartel-oct2017.cfg -t configs/sps-tb-oct2017.cfg
./Judith -c coarseAlign -i $CONV.root -r configs/sps-kartel-june2017.cfg -t configs/sps-tb-june2017.cfg
#./Judith -c fineAlign -i $CONV.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg
### end ##
#./Judith -c process -i $CONV.root -o $ANA-p.root -r configs/sps-kartel-oct2017.cfg -t configs/sps-tb-oct2017.cfg -R $ANA-results.root
./Judith -c process -i $CONV.root -o $ANA-p.root -r configs/sps-kartel-june2017.cfg -t configs/sps-tb-june2017.cfg -R $ANA-results.root

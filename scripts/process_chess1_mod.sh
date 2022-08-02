#!/bin/sh

RAW="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/raw/acq$1"
CONV="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/conv/ref$1" 
ANA="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/ana/ref$1" 
DUT="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/conv/run$1"
DUTANA="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/ana/run$1"
#DUTANA="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/run$1"


#./Judith -c convert -i $RAW.bin -o $CONV.root 
./Judith -c synchronize -i $CONV.root -I $DUT.root -o $ANA-s.root -O $DUTANA-s.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb2.cfg
### Alignment ###
#./Judith -c coarseAlign -i $ANA-s.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg -n 10000
#./Judith -c fineAlign -i $ANA-s.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg -n 100000
#./Judith -c coarseAlignDUT -i $ANA-s.root -I $DUTANA-s.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb-dut.cfg
#./Judith -c fineAlignDUT -i $ANA-s.root -I $DUTANA-s.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb-dut.cfg
### end ##
#./Judith -c process -i $ANA-s.root -o $ANA-sp.root -r configs/sps-kartel.cfg -t configs/sps-tb2.cfg -R $ANA-results.root
#./Judith -c process -i $DUTANA-s.root -o $DUTANA-sp.root -r configs/sps-dut-FEI4-C.cfg -t configs/sps-tb-dut.cfg -R $DUTANA-results.root
#./Judith -c analysisDUT -i $REFANA-sp.root -I $DUTANA-sp.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb.cfg -R $DUTANA-analysis.root


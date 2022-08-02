#!/bin/sh

REF1="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/acq$1"
REF="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/ref$1" 
DUT="/afs/cern.ch/work/b/bojan/Judith/data/june_2016/run$1"

./Judith -c convert -i $REF1.bin -o $REF1.root 
#./Judith -c synchronize -i $REF.root -I $DUT.root -o $REF-s.root -O $DUT-s.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb2.cfg
### Alignment ###
#./Judith -c coarseAlign -i $REF-s.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg -n 10000
#./Judith -c fineAlign -i $REF-s.root -r configs/sps-kartel.cfg -t configs/sps-tb.cfg -n 100000
#./Judith -c coarseAlignDUT -i $REF-s.root -I $DUT-s.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb-dut.cfg
#./Judith -c fineAlignDUT -i $REF-s.root -I $DUT-s.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb-dut.cfg
### end ##
#./Judith -c process -i $REF-s.root -o $REF-sp.root -r configs/sps-kartel.cfg -t configs/sps-tb2.cfg -R $REF-results.root
#./Judith -c process -i $DUT-s.root -o $DUT-sp.root -r configs/sps-dut-FEI4-C.cfg -t configs/sps-tb-dut.cfg -R $DUT-results.root
#./Judith -c analysisDUT -i $REF-sp.root -I $DUT-sp.root -r configs/sps-kartel.cfg -d configs/sps-dut-FEI4-C.cfg -t configs/sps-tb.cfg -R $DUT-analysis.root


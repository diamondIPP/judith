CC = g++
CFLAGS = `root-config --cflags` -O2 -Wall
LFLAGS = `root-config --ldflags --glibs` -O1
BUILD_DIR := build
SRCPATH = src
EXECUTABLE = Judith
OBJECTS = $(BUILD_DIR)/configparser.o $(BUILD_DIR)/inputargs.o $(BUILD_DIR)/main.o $(BUILD_DIR)/examplelooper.o $(BUILD_DIR)/configloopers.o $(BUILD_DIR)/synchronize.o $(BUILD_DIR)/coarsealign.o $(BUILD_DIR)/processevents.o $(BUILD_DIR)/finealigndut.o $(BUILD_DIR)/looper.o $(BUILD_DIR)/coarsealigndut.o $(BUILD_DIR)/noisescan.o $(BUILD_DIR)/analysis.o $(BUILD_DIR)/finealign.o $(BUILD_DIR)/applymask.o $(BUILD_DIR)/analysisdut.o $(BUILD_DIR)/kartelconvert.o $(BUILD_DIR)/storageio.o $(BUILD_DIR)/event.o $(BUILD_DIR)/hit.o $(BUILD_DIR)/cluster.o $(BUILD_DIR)/track.o $(BUILD_DIR)/plane.o $(BUILD_DIR)/occupancy.o $(BUILD_DIR)/singleanalyzer.o $(BUILD_DIR)/hitinfo.o $(BUILD_DIR)/correlation.o $(BUILD_DIR)/dualanalyzer.o $(BUILD_DIR)/dutcorrelation.o $(BUILD_DIR)/depiction.o $(BUILD_DIR)/configanalyzers.o $(BUILD_DIR)/matching.o $(BUILD_DIR)/examplesingleanalyzer.o $(BUILD_DIR)/clusterinfo.o $(BUILD_DIR)/efficiency.o $(BUILD_DIR)/trackinfo.o $(BUILD_DIR)/dutdepiction.o $(BUILD_DIR)/dutresiduals.o $(BUILD_DIR)/residuals.o $(BUILD_DIR)/eventinfo.o $(BUILD_DIR)/syncfluctuation.o $(BUILD_DIR)/exampledualanalyzer.o $(BUILD_DIR)/sensor.o $(BUILD_DIR)/device.o $(BUILD_DIR)/alignment.o $(BUILD_DIR)/noisemask.o $(BUILD_DIR)/configmechanics.o $(BUILD_DIR)/trackmatcher.o $(BUILD_DIR)/trackmaker.o $(BUILD_DIR)/clustermaker.o $(BUILD_DIR)/eventdepictor.o $(BUILD_DIR)/configprocessors.o $(BUILD_DIR)/synchronizer.o $(BUILD_DIR)/processors.o $(BUILD_DIR)/largesynchronizer.o 

$(info *** ROOTSYS:** $(ROOTSYS) ***)

all: Judith

Judith: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o $(EXECUTABLE)

$(BUILD_DIR)/configparser.o: $(SRCPATH)/configparser.cpp
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $(SRCPATH)/configparser.cpp -o $(BUILD_DIR)/configparser.o

$(BUILD_DIR)/inputargs.o: $(SRCPATH)/inputargs.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/inputargs.cpp -o $(BUILD_DIR)/inputargs.o

$(BUILD_DIR)/main.o: $(SRCPATH)/main.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/main.cpp -o $(BUILD_DIR)/main.o

$(BUILD_DIR)/examplelooper.o: $(SRCPATH)/loopers/examplelooper.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/examplelooper.cpp -o $(BUILD_DIR)/examplelooper.o

$(BUILD_DIR)/configloopers.o: $(SRCPATH)/loopers/configloopers.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/configloopers.cpp -o $(BUILD_DIR)/configloopers.o

$(BUILD_DIR)/synchronize.o: $(SRCPATH)/loopers/synchronize.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/synchronize.cpp -o $(BUILD_DIR)/synchronize.o

$(BUILD_DIR)/coarsealign.o: $(SRCPATH)/loopers/coarsealign.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/coarsealign.cpp -o $(BUILD_DIR)/coarsealign.o

$(BUILD_DIR)/processevents.o: $(SRCPATH)/loopers/processevents.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/processevents.cpp -o $(BUILD_DIR)/processevents.o

$(BUILD_DIR)/finealigndut.o: $(SRCPATH)/loopers/finealigndut.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/finealigndut.cpp -o $(BUILD_DIR)/finealigndut.o

$(BUILD_DIR)/looper.o: $(SRCPATH)/loopers/looper.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/looper.cpp -o $(BUILD_DIR)/looper.o

$(BUILD_DIR)/coarsealigndut.o: $(SRCPATH)/loopers/coarsealigndut.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/coarsealigndut.cpp -o $(BUILD_DIR)/coarsealigndut.o

$(BUILD_DIR)/noisescan.o: $(SRCPATH)/loopers/noisescan.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/noisescan.cpp -o $(BUILD_DIR)/noisescan.o

$(BUILD_DIR)/analysis.o: $(SRCPATH)/loopers/analysis.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/analysis.cpp -o $(BUILD_DIR)/analysis.o

$(BUILD_DIR)/finealign.o: $(SRCPATH)/loopers/finealign.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/finealign.cpp -o $(BUILD_DIR)/finealign.o

$(BUILD_DIR)/applymask.o: $(SRCPATH)/loopers/applymask.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/applymask.cpp -o $(BUILD_DIR)/applymask.o

$(BUILD_DIR)/analysisdut.o: $(SRCPATH)/loopers/analysisdut.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/loopers/analysisdut.cpp -o $(BUILD_DIR)/analysisdut.o

$(BUILD_DIR)/kartelconvert.o: $(SRCPATH)/converters/kartelconvert.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/converters/kartelconvert.cpp -o $(BUILD_DIR)/kartelconvert.o

$(BUILD_DIR)/storageio.o: $(SRCPATH)/storage/storageio.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/storage/storageio.cpp -o $(BUILD_DIR)/storageio.o

$(BUILD_DIR)/event.o: $(SRCPATH)/storage/event.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/storage/event.cpp -o $(BUILD_DIR)/event.o

$(BUILD_DIR)/hit.o: $(SRCPATH)/storage/hit.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/storage/hit.cpp -o $(BUILD_DIR)/hit.o

$(BUILD_DIR)/cluster.o: $(SRCPATH)/storage/cluster.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/storage/cluster.cpp -o $(BUILD_DIR)/cluster.o

$(BUILD_DIR)/track.o: $(SRCPATH)/storage/track.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/storage/track.cpp -o $(BUILD_DIR)/track.o

$(BUILD_DIR)/plane.o: $(SRCPATH)/storage/plane.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/storage/plane.cpp -o $(BUILD_DIR)/plane.o

$(BUILD_DIR)/occupancy.o: $(SRCPATH)/analyzers/occupancy.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/occupancy.cpp -o $(BUILD_DIR)/occupancy.o

$(BUILD_DIR)/singleanalyzer.o: $(SRCPATH)/analyzers/singleanalyzer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/singleanalyzer.cpp -o $(BUILD_DIR)/singleanalyzer.o

$(BUILD_DIR)/hitinfo.o: $(SRCPATH)/analyzers/hitinfo.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/hitinfo.cpp -o $(BUILD_DIR)/hitinfo.o

$(BUILD_DIR)/correlation.o: $(SRCPATH)/analyzers/correlation.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/correlation.cpp -o $(BUILD_DIR)/correlation.o

$(BUILD_DIR)/dualanalyzer.o: $(SRCPATH)/analyzers/dualanalyzer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/dualanalyzer.cpp -o $(BUILD_DIR)/dualanalyzer.o

$(BUILD_DIR)/dutcorrelation.o: $(SRCPATH)/analyzers/dutcorrelation.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/dutcorrelation.cpp -o $(BUILD_DIR)/dutcorrelation.o

$(BUILD_DIR)/depiction.o: $(SRCPATH)/analyzers/depiction.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/depiction.cpp -o $(BUILD_DIR)/depiction.o

$(BUILD_DIR)/configanalyzers.o: $(SRCPATH)/analyzers/configanalyzers.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/configanalyzers.cpp -o $(BUILD_DIR)/configanalyzers.o

$(BUILD_DIR)/matching.o: $(SRCPATH)/analyzers/matching.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/matching.cpp -o $(BUILD_DIR)/matching.o

$(BUILD_DIR)/examplesingleanalyzer.o: $(SRCPATH)/analyzers/examplesingleanalyzer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/examplesingleanalyzer.cpp -o $(BUILD_DIR)/examplesingleanalyzer.o

$(BUILD_DIR)/clusterinfo.o: $(SRCPATH)/analyzers/clusterinfo.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/clusterinfo.cpp -o $(BUILD_DIR)/clusterinfo.o

$(BUILD_DIR)/efficiency.o: $(SRCPATH)/analyzers/efficiency.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/efficiency.cpp -o $(BUILD_DIR)/efficiency.o

$(BUILD_DIR)/trackinfo.o: $(SRCPATH)/analyzers/trackinfo.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/trackinfo.cpp -o $(BUILD_DIR)/trackinfo.o

$(BUILD_DIR)/dutdepiction.o: $(SRCPATH)/analyzers/dutdepiction.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/dutdepiction.cpp -o $(BUILD_DIR)/dutdepiction.o

$(BUILD_DIR)/dutresiduals.o: $(SRCPATH)/analyzers/dutresiduals.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/dutresiduals.cpp -o $(BUILD_DIR)/dutresiduals.o

$(BUILD_DIR)/residuals.o: $(SRCPATH)/analyzers/residuals.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/residuals.cpp -o $(BUILD_DIR)/residuals.o

$(BUILD_DIR)/eventinfo.o: $(SRCPATH)/analyzers/eventinfo.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/eventinfo.cpp -o $(BUILD_DIR)/eventinfo.o

$(BUILD_DIR)/syncfluctuation.o: $(SRCPATH)/analyzers/syncfluctuation.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/syncfluctuation.cpp -o $(BUILD_DIR)/syncfluctuation.o

$(BUILD_DIR)/exampledualanalyzer.o: $(SRCPATH)/analyzers/exampledualanalyzer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/analyzers/exampledualanalyzer.cpp -o $(BUILD_DIR)/exampledualanalyzer.o

$(BUILD_DIR)/sensor.o: $(SRCPATH)/mechanics/sensor.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/mechanics/sensor.cpp -o $(BUILD_DIR)/sensor.o

$(BUILD_DIR)/device.o: $(SRCPATH)/mechanics/device.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/mechanics/device.cpp -o $(BUILD_DIR)/device.o

$(BUILD_DIR)/alignment.o: $(SRCPATH)/mechanics/alignment.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/mechanics/alignment.cpp -o $(BUILD_DIR)/alignment.o

$(BUILD_DIR)/noisemask.o: $(SRCPATH)/mechanics/noisemask.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/mechanics/noisemask.cpp -o $(BUILD_DIR)/noisemask.o

$(BUILD_DIR)/configmechanics.o: $(SRCPATH)/mechanics/configmechanics.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/mechanics/configmechanics.cpp -o $(BUILD_DIR)/configmechanics.o

$(BUILD_DIR)/trackmatcher.o: $(SRCPATH)/processors/trackmatcher.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/trackmatcher.cpp -o $(BUILD_DIR)/trackmatcher.o

$(BUILD_DIR)/trackmaker.o: $(SRCPATH)/processors/trackmaker.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/trackmaker.cpp -o $(BUILD_DIR)/trackmaker.o

$(BUILD_DIR)/clustermaker.o: $(SRCPATH)/processors/clustermaker.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/clustermaker.cpp -o $(BUILD_DIR)/clustermaker.o

$(BUILD_DIR)/eventdepictor.o: $(SRCPATH)/processors/eventdepictor.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/eventdepictor.cpp -o $(BUILD_DIR)/eventdepictor.o

$(BUILD_DIR)/configprocessors.o: $(SRCPATH)/processors/configprocessors.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/configprocessors.cpp -o $(BUILD_DIR)/configprocessors.o

$(BUILD_DIR)/synchronizer.o: $(SRCPATH)/processors/synchronizer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/synchronizer.cpp -o $(BUILD_DIR)/synchronizer.o

$(BUILD_DIR)/processors.o: $(SRCPATH)/processors/processors.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/processors.cpp -o $(BUILD_DIR)/processors.o

$(BUILD_DIR)/largesynchronizer.o: $(SRCPATH)/processors/largesynchronizer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)/processors/largesynchronizer.cpp -o $(BUILD_DIR)/largesynchronizer.o

clean:
	 rm $(BUILD_DIR)/*.o Judith
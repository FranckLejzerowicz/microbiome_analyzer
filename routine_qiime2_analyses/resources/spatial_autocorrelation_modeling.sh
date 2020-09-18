#!/bin/bash

java -cp SPATIALAUTOCORRELATIONPATH/bin/Autocorrelation.jar \
  edu.ucsf.SpatialAutocorrelation.SpatialAutocorrelationLauncher \
	--sBIOMPath=BIOM \
	--sOutputPath=OUTPUTPATH \
	--iMCMCChains=10 \
	--iMCMCIterations=10000 \
	--rgsDistanceNeighborhoods='0-150' \
	--bOutputData=true \
	--iPrevalenceMinimum=0 \
	--bNormalize=false \
	--iInitialScreeningIterations=10

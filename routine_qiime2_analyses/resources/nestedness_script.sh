#!/bin/bash

sIODir=$HOME/softs/SpatialAutocorrelation
sJavaDir=$sIODir/bin
sBiomPath=BIOM
sOutputDir=OUTPUTDIR
sTaxonRank=TAXONRANK
# e.g. otu
sNullModel=NULLMODEL
# e.g. equiprobablefixed
sAxis=AXIS
# e.g. sample
sComparisonMode=COMPARISONMODE
# e.g. betweeneachpairoftypes
sMetadataField=METADATAFIELD
# e.g. #SampleID
iPrevalenceMinimum=1
bSimulatefalse=false

#making nestedness graph
#java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Grapher.GrapherLauncher --help > $sIODir/doc/Autocorrelation.edu.ucsf.Nestedness.Grapher.GrapherLauncher.txt
java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Grapher.GrapherLauncher \
--sBIOMPath=$sBiomPath \
--bNormalize=false \
--sTaxonRank=$sTaxonRank \
--sOutputPath=$sOutputDir/graphs.csv \
--rgsSampleMetadataFields=area,province

#loading comparisons
#java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher --help > $sIODir/doc/Autocorrelation.edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher.txt
java -Xmx5g -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher \
--sBIOMPath=$sBiomPath \
--sOutputPath=$sOutputDir/comparisons.csv \
--bNormalize=false \
--sTaxonRank=$sTaxonRank \
--sMetadataField=$sMetadataField \
--iRandomSeed=1234 \
--sComparisonMode=$sComparisonMode \
--iNestednessPairs=250 \
--sNestednessAxis=$sAxis \
--iPrevalenceMinimum=1

#running statistics
#java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Calculator.CalculatorLauncher --help > $sIODir/doc/Autocorrelation.edu.ucsf.Nestedness.Calculator.CalculatorLauncher.txt
java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Calculator.CalculatorLauncher \
--sBIOMPath=$sBiomPath \
--sOutputPath=$sOutputDir/statistics.csv \
--bNormalize=false \
--sTaxonRank=$sTaxonRank \
--sComparisonsPath=$sOutputDir/comparisons.csv \
--iNullModelIterations=10000 \
--bOrderedNODF=false \
--sNestednessAxis=$sAxis \
--sNestednessNullModel=$sNullModel \
--iPrevalenceMinimum=$iPrevalenceMinimum \
--bSimulate=$bSimulate

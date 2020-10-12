#!/bin/bash

sJavaPath=SOFT
sBiomPath=BIOMPATH
sOutputDir=OUTPUTDIR
sNullModel=NULLMODEL
sAxis=AXIS
sComparisonMode=COMPARISONMODE
rgsSampleMetadataFields=METAS
iPrevalenceMinimum=1
bSimulatefalse=false

#making nestedness graph
java -cp $sJavaPath edu.ucsf.Nestedness.Grapher.GrapherLauncher --sBIOMPath=$sBiomPath --bNormalize=false --sTaxonRank=otu --sOutputPath=$sOutputDir/graphs.csv --rgsSampleMetadataFields=$rgsSampleMetadataFields
#loading comparisons
java -Xmx5g -cp $sJavaPath edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher --sBIOMPath=$sBiomPath --sOutputPath=$sOutputDir/comparisons.csv --bNormalize=false --sTaxonRank=otu --sMetadataField=METADATAFIELD --iRandomSeed=1234 --sComparisonMode=$sComparisonMode --iNestednessPairs=250 --sNestednessAxis=$sAxis --iPrevalenceMinimum=1
#running statistics
java -cp $sJavaPath edu.ucsf.Nestedness.Calculator.CalculatorLauncher --sBIOMPath=$sBiomPath --sOutputPath=$sOutputDir/statistics.csv --bNormalize=false --sTaxonRank=otu --sComparisonsPath=$sOutputDir/comparisons.csv --iNullModelIterations=10000 --bOrderedNODF=false --sNestednessAxis=$sAxis --sNestednessNullModel=$sNullModel --iPrevalenceMinimum=$iPrevalenceMinimum --bSimulate=$bSimulate

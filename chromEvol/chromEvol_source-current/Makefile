#! /usr/local/bin/gmake
#  $Id: Makefile,v 1.2 2006/07/19 15:22:21 eyalprivman Exp $

.PHONY: all lib clean exec

all: exec

CXX= g++
CC=$(CXX)

CXXFLAGS += -O3 -Wno-deprecated -DDOUBLEREP
Libsources= AddLog.cpp aaJC.cpp alphabet.cpp amino.cpp			\
   bestAlpha.cpp bestHKYparam.cpp bootstrap.cpp			\
   bblEM.cpp bblEMProportionalEB.cpp bblLSProportionalEB.cpp chebyshevAccelerator.cpp clustalFormat.cpp codon.cpp codonJC.cpp			\
   computeCounts.cpp computeDownAlg.cpp computeMarginalAlg.cpp				\
   computePijComponent.cpp computeUpAlg.cpp computeUpAlgFactors.cpp			\
   ConversionUtils.cpp countTableComponent.cpp datMatrixHolder.cpp distanceTable.cpp	\
   distribution.cpp errorMsg.cpp evaluateCharacterFreq.cpp				\
   fastStartTree.cpp fastaFormat.cpp fromCountTableComponentToDistance.cpp   fromCountTableComponentToDistancePropEB.cpp fromQtoPt.cpp				\
   generalGammaDistributionFixedCategories.cpp gammaDistribution.cpp gammaUtilities.cpp	\
   generalGammaDistribution.cpp getRandomWeights.cpp goldmanYangModel.cpp		\
   granthamChemicalDistances.cpp hky.cpp likeDist.cpp  \
   likeDistProp.cpp likelihoodComputation.cpp						\
   likelihoodComputationFactors.cpp logFile.cpp maseFormat.cpp				\
   molphyFormat.cpp nexusFormat.cpp nj.cpp njConstrain.cpp				\
   nucJC.cpp nucleotide.cpp numRec.cpp Parameters.cpp phylipFormat.cpp			\
   pijAccelerator.cpp readDatMatrix.cpp readTree.cpp recognizeFormat.cpp		\
   replacementModel.cpp seqContainerTreeMap.cpp			\
   sequence.cpp sequenceContainer.cpp simulateTree.cpp					\
   siteSpecificRate.cpp someUtil.cpp split.cpp splitMap.cpp				\
   splitTreeUtil.cpp stochasticProcess.cpp suffStatComponent.cpp			\
   talRandom.cpp tree.cpp treeIt.cpp treeUtil.cpp uniDistribution.cpp			\
   uniformDistribution.cpp givenRatesMLDistance.cpp distanceBasedSeqs2Tree.cpp					\
   posteriorDistance.cpp pairwiseGammaDistance.cpp doubleRep.cpp			\
   logRep.cpp indel.cpp indelModel.cpp mulAlphabet.cpp					\
   matrixUtils.cpp geneticCodeHolder.cpp				\
   samplingSequences.cpp normalDist.cpp				\
   tamura92.cpp bestTamura92param.cpp phylipSequentialFormat.cpp			\
   multipleStochasticProcess.cpp distributionPlusInvariant.cpp\
   gammaDistributionFixedCategories.cpp  generalGammaDistributionPlusInvariant.cpp  gammaDistributionPlusInvariant.cpp   \
   distributionPlusCategory.cpp simulateJumps.cpp seqeuncesFilter.cpp 					\
   gainLossAlphabet.cpp wYangModel.cpp codonUtils.cpp likelihoodComputation2Codon.cpp likeDist2Codon.cpp unObservableData.cpp likelihoodComputationGL.cpp	\
	chrNumberOptions.cpp chrNumOptimizer.cpp \
	evalParamChrNum.cpp chrNumberMng.cpp chrNumModel.cpp dataUtils.cpp finiteIntAlphabet.cpp \
	chrCountFormat.cpp simulateChangesAlongTree.cpp computeUpAlgChrNum.cpp simulateJumpsLargeAlphabet.cpp generalChrNumModel.cpp evalPPDistMng.cpp \
	simulatePPInPopulation.cpp


LIBNAME = chromEvol

lib: lib$(LIBNAME).a

lib$(LIBNAME).a: $(Libsources:.cpp=.o)
	ar rv $@ $?
#	ranlib $@

# LibCsources= cmdline.c
# LibCsources += getopt.c getopt1.c 

EXEC = chromEvol

$(EXEC): lib$(LIBNAME).a

exec: $(EXEC)

clean:
	rm *.o lib$(LIBNAME).a

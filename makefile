C		= gcc
CC		= g++
LD		= g++
DEBUG		= -g
OPTIMIZE	= -O
BOOM		= BOOM
GSL		= BOOM/GSL
CFLAGS		= $(OPTIMIZE) -fpermissive -w -D_THREAD_SAFE -D_REENTRANT -D__STL_PTHREADS -D_MT -D_PTHREADS
LDFLAGS		= $(OPTIMIZE) -lpthread
OBJ		= obj
LIBS		= -L$(BOOM) -lBOOM
AR		= ar -s

$(OBJ):
	mkdir $(OBJ)

clean:
	rm obj/*.o

$(BOOM)/libBOOM.a:
	make $(BOOM)/libBOOM.a

#---------------------------------------------------------
$(OBJ)/GslLabeledMatrixLoader.o:\
		$(GSL)/LabeledMatrixLoader.H \
		$(GSL)/LabeledMatrixLoader.C
	$(CC) $(CFLAGS) -o $(OBJ)/GslLabeledMatrixLoader.o -c \
		$(GSL)/LabeledMatrixLoader.C
#---------------------------------------------------------
$(OBJ)/GslMatrix.o:\
		$(GSL)/Matrix.H \
		$(GSL)/Matrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/GslMatrix.o -c \
		$(GSL)/Matrix.C
#---------------------------------------------------------
$(OBJ)/GslOptimizer.o:\
		$(GSL)/Optimizer.H \
		$(GSL)/Optimizer.C
	$(CC) $(CFLAGS) -o $(OBJ)/GslOptimizer.o -c \
		$(GSL)/Optimizer.C
#---------------------------------------------------------
$(OBJ)/GslVector.o:\
		$(GSL)/Vector.H \
		$(GSL)/Vector.C
	$(CC) $(CFLAGS) -o $(OBJ)/GslVector.o -c \
		$(GSL)/Vector.C
#---------------------------------------------------------
$(OBJ)/mud.o:\
		mud.C
	$(CC) $(CFLAGS) -o $(OBJ)/mud.o -c \
		mud.C
#---------------------------------------------------------
mud: \
		$(OBJ)/CommandLine.o \
		$(OBJ)/String.o \
		$(OBJ)/StrTokenizer.o \
		$(OBJ)/mud.o
	$(CC) $(LDFLAGS) -o mud \
		$(OBJ)/mud.o $(LIBS)
#---------------------------------------------------------
$(OBJ)/Message.o:\
		Message.H \
		Message.C
	$(CC) $(CFLAGS) -o $(OBJ)/Message.o -c \
		Message.C
#---------------------------------------------------------
$(OBJ)/FitchParsimony.o:\
		FitchParsimony.H \
		FitchParsimony.C
	$(CC) $(CFLAGS) -o $(OBJ)/FitchParsimony.o -c \
		FitchParsimony.C
#---------------------------------------------------------
$(OBJ)/NmerRateMatrix.o:\
		NmerRateMatrix.H \
		NmerRateMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/NmerRateMatrix.o -c \
		NmerRateMatrix.C
#---------------------------------------------------------
$(OBJ)/NmerSubstMatrix.o:\
		NmerSubstMatrix.H \
		NmerSubstMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/NmerSubstMatrix.o -c \
		NmerSubstMatrix.C
#---------------------------------------------------------
$(OBJ)/NthOrdRateMatrix.o:\
		NthOrdRateMatrix.H \
		NthOrdRateMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/NthOrdRateMatrix.o -c \
		NthOrdRateMatrix.C
#---------------------------------------------------------
$(OBJ)/NthOrdSubstMatrix.o:\
		NthOrdSubstMatrix.H \
		NthOrdSubstMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/NthOrdSubstMatrix.o -c \
		NthOrdSubstMatrix.C
#---------------------------------------------------------
$(OBJ)/RateMatrix.o:\
		RateMatrix.H \
		RateMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/RateMatrix.o -c \
		RateMatrix.C
#---------------------------------------------------------
$(OBJ)/AlignmentNmerTable.o:\
		AlignmentNmerTable.H \
		AlignmentNmerTable.C
	$(CC) $(CFLAGS) -o $(OBJ)/AlignmentNmerTable.o -c \
		AlignmentNmerTable.C
#---------------------------------------------------------
$(OBJ)/GapRateMatrix.o:\
		GapRateMatrix.H \
		GapRateMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/GapRateMatrix.o -c \
		GapRateMatrix.C
#---------------------------------------------------------
$(OBJ)/RateMatrixType.o:\
		RateMatrixType.H \
		RateMatrixType.C
	$(CC) $(CFLAGS) -o $(OBJ)/RateMatrixType.o -c \
		RateMatrixType.C
#---------------------------------------------------------
$(OBJ)/SubstitutionMatrix.o:\
		SubstitutionMatrix.H \
		SubstitutionMatrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/SubstitutionMatrix.o -c \
		SubstitutionMatrix.C
#---------------------------------------------------------
$(OBJ)/eigen-decompose.o:\
		eigen-decompose.C
	$(CC) $(CFLAGS) -o $(OBJ)/eigen-decompose.o -c \
		eigen-decompose.C
#---------------------------------------------------------
eigen-decompose: \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/eigen-decompose.o
	$(CC) $(LDFLAGS) -o eigen-decompose \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/eigen-decompose.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/JukesCantor.o:\
		JukesCantor.C\
		JukesCantor.H
	$(CC) $(CFLAGS) -o $(OBJ)/JukesCantor.o -c \
		JukesCantor.C
#--------------------------------------------------------
$(OBJ)/Kimura2P.o:\
		Kimura2P.C\
		Kimura2P.H
	$(CC) $(CFLAGS) -o $(OBJ)/Kimura2P.o -c \
		Kimura2P.C
#--------------------------------------------------------
$(OBJ)/Kimura2Param.o:\
		Kimura2Param.C\
		Kimura2Param.H
	$(CC) $(CFLAGS) -o $(OBJ)/Kimura2Param.o -c \
		Kimura2Param.C
#--------------------------------------------------------
$(OBJ)/FEL.o:\
		FEL.C\
		FEL.H
	$(CC) $(CFLAGS) -o $(OBJ)/FEL.o -c \
		FEL.C
#--------------------------------------------------------
$(OBJ)/HalpernBruno.o:\
		HalpernBruno.C\
		HalpernBruno.H
	$(CC) $(CFLAGS) -o $(OBJ)/HalpernBruno.o -c \
		HalpernBruno.C
#--------------------------------------------------------
$(OBJ)/HKY.o:\
		HKY.C\
		HKY.H
	$(CC) $(CFLAGS) -o $(OBJ)/HKY.o -c \
		HKY.C
#--------------------------------------------------------
$(OBJ)/REV.o:\
		REV.C\
		REV.H
	$(CC) $(CFLAGS) -o $(OBJ)/REV.o -c \
		REV.C
#--------------------------------------------------------
$(OBJ)/Phylogeny.o:\
		Phylogeny.C\
		Phylogeny.H
	$(CC) $(CFLAGS) -o $(OBJ)/Phylogeny.o -c \
		Phylogeny.C
#---------------------------------------------------------
$(OBJ)/UPGMA.o:\
		UPGMA.C
	$(CC) $(CFLAGS) -o $(OBJ)/UPGMA.o -c \
		UPGMA.C
#---------------------------------------------------------
UPGMA: \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/UPGMA.o
	$(CC) $(LDFLAGS) -o UPGMA \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/UPGMA.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
$(OBJ)/neighbor-joining.o:\
		neighbor-joining.C
	$(CC) $(CFLAGS) -o $(OBJ)/neighbor-joining.o -c \
		neighbor-joining.C
#---------------------------------------------------------
neighbor-joining: \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/neighbor-joining.o
	$(CC) $(LDFLAGS) -o neighbor-joining \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/neighbor-joining.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
$(OBJ)/train-parallel-H.o:\
		train-parallel-H.C
	mpiCC $(CFLAGS) -o $(OBJ)/train-parallel-H.o -c \
		train-parallel-H.C
#---------------------------------------------------------                      
train-parallel-H: \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/HalpernBruno.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/Message.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/train-parallel-H.o
	mpiCC $(LDFLAGS) -o train-parallel-H \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/Message.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/HalpernBruno.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/train-parallel-H.o \
	-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/train-phylohmm.o:\
		train-phylohmm.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-phylohmm.o -c \
		train-phylohmm.C
#---------------------------------------------------------
train-phylohmm: \
		$(BOOM)/libBOOM.a \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslOptimizer.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NthOrdFelsenstein.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/train-phylohmm.o
	$(LD) $(LDFLAGS) -o train-phylohmm \
		$(OBJ)/ContextType.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/NthOrdFelsenstein.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslOptimizer.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/train-phylohmm.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/FelsensteinsAlgorithm.o:\
		FelsensteinsAlgorithm.C\
		FelsensteinsAlgorithm.H
	$(CC) $(CFLAGS) -o $(OBJ)/FelsensteinsAlgorithm.o -c \
		FelsensteinsAlgorithm.C
#--------------------------------------------------------
$(OBJ)/RCO_Felsenstein.o:\
		RCO_Felsenstein.C \
		RCO_Felsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/RCO_Felsenstein.o -c \
		RCO_Felsenstein.C
#--------------------------------------------------------
$(OBJ)/LCO_Felsenstein.o:\
		LCO_Felsenstein.C \
		LCO_Felsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/LCO_Felsenstein.o -c \
		LCO_Felsenstein.C
#--------------------------------------------------------
$(OBJ)/NmerFelsenstein.o:\
		NmerFelsenstein.C \
		NmerFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/NmerFelsenstein.o -c \
		NmerFelsenstein.C
#--------------------------------------------------------
$(OBJ)/ACO_Felsenstein.o:\
		ACO_Felsenstein.C \
		ACO_Felsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/ACO_Felsenstein.o -c \
		ACO_Felsenstein.C
#--------------------------------------------------------
$(OBJ)/HOG_Felsenstein.o:\
		HOG_Felsenstein.C \
		HOG_Felsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/HOG_Felsenstein.o -c \
		HOG_Felsenstein.C
#--------------------------------------------------------
$(OBJ)/TRCO_Felsenstein.o:\
		TRCO_Felsenstein.C \
		TRCO_Felsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/TRCO_Felsenstein.o -c \
		TRCO_Felsenstein.C
#--------------------------------------------------------
$(OBJ)/FitchFelsenstein.o:\
		FitchFelsenstein.C \
		FitchFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/FitchFelsenstein.o -c \
		FitchFelsenstein.C
#---------------------------------------------------------
$(OBJ)/test-felsenstein.o:\
		test-felsenstein.C
	$(CC) $(CFLAGS) -o $(OBJ)/test-felsenstein.o -c \
		test-felsenstein.C
#---------------------------------------------------------
test-felsenstein: \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/FelsensteinsAlgorithm.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/test-felsenstein.o
	$(CC) $(LDFLAGS) -o test-felsenstein \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/FelsensteinsAlgorithm.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/test-felsenstein.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/reroot.o:\
		reroot.C
	$(CC) $(CFLAGS) -o $(OBJ)/reroot.o -c \
		reroot.C
#---------------------------------------------------------
reroot: \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/reroot.o
	$(CC) $(LDFLAGS) -o reroot \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/reroot.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/alignment-likelihood.o:\
		alignment-likelihood.C
	$(CC) $(CFLAGS) -o $(OBJ)/alignment-likelihood.o -c \
		alignment-likelihood.C
#---------------------------------------------------------
alignment-likelihood: \
		$(OBJ)/AlignmentNmerTable.o \
		$(OBJ)/GapRateMatrix.o \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/NmerFelsenstein.o \
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/alignment-likelihood.o
	$(CC) $(LDFLAGS) -o alignment-likelihood \
		$(OBJ)/AlignmentNmerTable.o \
		$(OBJ)/GapRateMatrix.o \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/NmerFelsenstein.o \
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/alignment-likelihood.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/build-PSA.o:\
		build-PSA.C
	$(CC) $(CFLAGS) -o $(OBJ)/build-PSA.o -c \
		build-PSA.C
#---------------------------------------------------------
build-PSA: \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/NmerFelsenstein.o \
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/build-PSA.o
	$(CC) $(LDFLAGS) -o build-PSA \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/NmerFelsenstein.o \
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/build-PSA.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------
#--------------------------------------------------------
$(OBJ)/IndelHistory.o:\
		IndelHistory.C\
		IndelHistory.H
	$(CC) $(CFLAGS) -o $(OBJ)/IndelHistory.o -c \
		IndelHistory.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/GapPattern.o:\
		GapPattern.C\
		GapPattern.H
	$(CC) $(CFLAGS) -o $(OBJ)/GapPattern.o -c \
		GapPattern.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/UsePattern.o:\
		UsePattern.C\
		UsePattern.H
	$(CC) $(CFLAGS) -o $(OBJ)/UsePattern.o -c \
		UsePattern.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/reduce-matrix-order.o:\
		reduce-matrix-order.C
	$(CC) $(CFLAGS) -o $(OBJ)/reduce-matrix-order.o -c \
		reduce-matrix-order.C
#---------------------------------------------------------
reduce-matrix-order: \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/reduce-matrix-order.o
	$(CC) $(LDFLAGS) -o reduce-matrix-order \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/reduce-matrix-order.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------


#--------------------------------------------------------
$(OBJ)/infer-ancestral-states.o:\
		infer-ancestral-states.C
	$(CC) $(CFLAGS) -o $(OBJ)/infer-ancestral-states.o -c \
		infer-ancestral-states.C
#---------------------------------------------------------
infer-ancestral-states: \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/infer-ancestral-states.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/FitchParsimony.o
	$(CC) $(LDFLAGS) -o infer-ancestral-states \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/infer-ancestral-states.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/FitchParsimony.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/compare-alignments.o:\
		compare-alignments.C
	$(CC) $(CFLAGS) -o $(OBJ)/compare-alignments.o -c \
		compare-alignments.C
#---------------------------------------------------------
compare-alignments: \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/compare-alignments.o
	$(CC) $(LDFLAGS) -o compare-alignments \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/compare-alignments.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------




#--------------------------------------------------------
$(OBJ)/ContextType.o:\
		ContextType.C\
		ContextType.H
	$(CC) $(CFLAGS) -o $(OBJ)/ContextType.o -c \
		ContextType.C
#---------------------------------------------------------


#--------------------------------------------------------
$(OBJ)/ACO-from-dual.o:\
		ACO-from-dual.C
	$(CC) $(CFLAGS) -o $(OBJ)/ACO-from-dual.o -c \
		ACO-from-dual.C
#---------------------------------------------------------
ACO-from-dual: \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/ACO-from-dual.o
	$(CC) $(LDFLAGS) -o ACO-from-dual \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/ACO-from-dual.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/NthOrdFelsenstein.o:\
		NthOrdFelsenstein.C\
		NthOrdFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/NthOrdFelsenstein.o -c \
		NthOrdFelsenstein.C
#---------------------------------------------------------
#--------------------------------------------------------
$(OBJ)/split-alignment.o:\
		split-alignment.C
	$(CC) $(CFLAGS) -o $(OBJ)/split-alignment.o -c \
		split-alignment.C
#---------------------------------------------------------
split-alignment: \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/split-alignment.o
	$(CC) $(LDFLAGS) -o split-alignment \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/split-alignment.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------
#---------------------------------------------------------
$(OBJ)/fix-maf-track-names.o:\
		fix-maf-track-names.C
	$(CC) $(CFLAGS) -o $(OBJ)/fix-maf-track-names.o -c \
		fix-maf-track-names.C
#---------------------------------------------------------
fix-maf-track-names: \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/fix-maf-track-names.o
	$(CC) $(LDFLAGS) -o fix-maf-track-names \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/fix-maf-track-names.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
$(OBJ)/first-n-alignments.o:\
		first-n-alignments.C
	$(CC) $(CFLAGS) -o $(OBJ)/first-n-alignments.o -c \
		first-n-alignments.C
#---------------------------------------------------------
first-n-alignments: \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/first-n-alignments.o
	$(CC) $(LDFLAGS) -o first-n-alignments \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/first-n-alignments.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------
$(OBJ)/alignment-stats.o:\
		alignment-stats.C
	$(CC) $(CFLAGS) -o $(OBJ)/alignment-stats.o -c \
		alignment-stats.C
#---------------------------------------------------------
alignment-stats: \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/alignment-stats.o
	$(CC) $(LDFLAGS) -o alignment-stats \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/alignment-stats.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/compare-to-ACO.o:\
		compare-to-ACO.C
	$(CC) $(CFLAGS) -o $(OBJ)/compare-to-ACO.o -c \
		compare-to-ACO.C
#---------------------------------------------------------
compare-to-ACO: \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/compare-to-ACO.o
	$(CC) $(LDFLAGS) -o compare-to-ACO \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/compare-to-ACO.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/PeriodicRateMatrix.o:\
		PeriodicRateMatrix.C\
		PeriodicRateMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/PeriodicRateMatrix.o -c \
		PeriodicRateMatrix.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/PeriodicSubstMatrix.o:\
		PeriodicSubstMatrix.C\
		PeriodicSubstMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/PeriodicSubstMatrix.o -c \
		PeriodicSubstMatrix.C
#---------------------------------------------------------


#--------------------------------------------------------
$(OBJ)/DegenerateDnaMatch.o:\
		DegenerateDnaMatch.C\
		DegenerateDnaMatch.H
	$(CC) $(CFLAGS) -o $(OBJ)/DegenerateDnaMatch.o -c \
		DegenerateDnaMatch.C
#---------------------------------------------------------


#---------------------------------------------------------                      
$(OBJ)/instantiate-matrix.o:\
		instantiate-matrix.C
	mpiCC $(CFLAGS) -o $(OBJ)/instantiate-matrix.o -c \
		instantiate-matrix.C
instantiate-matrix: \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/GslOptimizer.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/GslMatrix.o \
	$(OBJ)/GslVector.o \
	$(OBJ)/GslLabeledMatrixLoader.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/Message.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/instantiate-matrix.o
	mpiCC $(LDFLAGS) -o instantiate-matrix \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/Message.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/GslOptimizer.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/GslMatrix.o \
	$(OBJ)/GslVector.o \
	$(OBJ)/GslLabeledMatrixLoader.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/instantiate-matrix.o \
	-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/PSA.o:\
		PSA.C\
		PSA.H
	$(CC) $(CFLAGS) -o $(OBJ)/PSA.o -c \
		PSA.C
#---------------------------------------------------------


#---------------------------------------------------------
$(OBJ)/sample-n-alignments.o:\
		sample-n-alignments.C
	$(CC) $(CFLAGS) -o $(OBJ)/sample-n-alignments.o -c \
		sample-n-alignments.C
#---------------------------------------------------------
sample-n-alignments: \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/sample-n-alignments.o
	$(CC) $(LDFLAGS) -o sample-n-alignments \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslLabeledMatrixLoader.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/sample-n-alignments.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/GapRewritingModel.o:\
		GapRewritingModel.C\
		GapRewritingModel.H
	$(CC) $(CFLAGS) -o $(OBJ)/GapRewritingModel.o -c \
		GapRewritingModel.C
#---------------------------------------------------------

#---------------------------------------------------------
$(OBJ)/train-parallel-wobble.o:\
		train-parallel-wobble.C
	mpiCC $(CFLAGS) -o $(OBJ)/train-parallel-wobble.o -c \
		train-parallel-wobble.C
train-parallel-wobble: \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/Message.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/train-parallel-wobble.o
	mpiCC $(LDFLAGS) -o train-parallel-wobble \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/Message.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/train-parallel-wobble.o \
	-lgsl -lm -lgslcblas $(LIBS)



#---------------------------------------------------------
$(OBJ)/train-parallel.o:\
		train-parallel.C
	mpiCC $(CFLAGS) -o $(OBJ)/train-parallel.o -c \
		train-parallel.C
train-parallel: \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/Message.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/train-parallel.o
	mpiCC $(LDFLAGS) -o train-parallel \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/Message.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/train-parallel.o \
	-lgsl -lm -lgslcblas $(LIBS)

#--------------------------------------------------------
$(OBJ)/filter-alignment-columns.o:\
		filter-alignment-columns.C
	$(CC) $(CFLAGS) -o $(OBJ)/filter-alignment-columns.o -c \
		filter-alignment-columns.C
#---------------------------------------------------------
filter-alignment-columns: \
	$(OBJ)/filter-alignment-columns.o
	$(CC) $(LDFLAGS) -o filter-alignment-columns \
	$(OBJ)/filter-alignment-columns.o \
	-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
$(OBJ)/extract-target-sequence.o:\
		extract-target-sequence.C
	$(CC) $(CFLAGS) -o $(OBJ)/extract-target-sequence.o -c \
		extract-target-sequence.C
#---------------------------------------------------------
extract-target-sequence: \
		$(OBJ)/extract-target-sequence.o
	$(CC) $(LDFLAGS) -o extract-target-sequence \
		$(OBJ)/extract-target-sequence.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/predict-phase.o:\
		predict-phase.C
	$(CC) $(CFLAGS) -o $(OBJ)/predict-phase.o -c \
		predict-phase.C
#---------------------------------------------------------
predict-phase: \
		$(OBJ)/AlignmentNmerTable.o \
		$(OBJ)/GapRateMatrix.o \
		$(OBJ)/GapRewritingModel.o \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NmerFelsenstein.o\
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/NthOrdFelsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/Message.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o \
		$(OBJ)/predict-phase.o
	$(CC) $(LDFLAGS) -o predict-phase \
		$(OBJ)/predict-phase.o \
		$(OBJ)/AlignmentNmerTable.o \
		$(OBJ)/GapRateMatrix.o \
		$(OBJ)/GapRewritingModel.o \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NmerFelsenstein.o\
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/NthOrdFelsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/Message.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------
libphylib.a:	$(OBJ) \
		$(OBJ)/HalpernBruno.o \
		$(OBJ)/SingleGainParsimony.o \
		$(OBJ)/RateMatrixDecomposed.o \
		$(OBJ)/SubstitutionMixtureModel.o \
		$(OBJ)/AlignmentNmerTable.o \
		$(OBJ)/GapRateMatrix.o \
		$(OBJ)/GapRewritingModel.o \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NmerFelsenstein.o\
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/NthOrdFelsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/Message.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o
	$(AR) -r libphylib.a \
		$(OBJ)/HalpernBruno.o \
		$(OBJ)/SingleGainParsimony.o \
		$(OBJ)/RateMatrixDecomposed.o \
		$(OBJ)/SubstitutionMixtureModel.o \
		$(OBJ)/AlignmentNmerTable.o \
		$(OBJ)/GapRateMatrix.o \
		$(OBJ)/GapRewritingModel.o \
		$(OBJ)/DegenerateDnaMatch.o \
		$(OBJ)/Message.o \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/ContextType.o \
		$(OBJ)/RateMatrixType.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/NmerFelsenstein.o\
		$(OBJ)/NmerRateMatrix.o \
		$(OBJ)/NmerSubstMatrix.o \
		$(OBJ)/ACO_Felsenstein.o \
		$(OBJ)/TRCO_Felsenstein.o \
		$(OBJ)/RCO_Felsenstein.o \
		$(OBJ)/LCO_Felsenstein.o \
		$(OBJ)/HOG_Felsenstein.o \
		$(OBJ)/FitchFelsenstein.o \
		$(OBJ)/FitchParsimony.o \
		$(OBJ)/NthOrdFelsenstein.o \
		$(OBJ)/SubstitutionMatrix.o \
		$(OBJ)/JukesCantor.o \
		$(OBJ)/Kimura2Param.o \
		$(OBJ)/HKY.o \
		$(OBJ)/FEL.o \
		$(OBJ)/REV.o \
		$(OBJ)/RateMatrix.o \
		$(OBJ)/NthOrdRateMatrix.o \
		$(OBJ)/NthOrdSubstMatrix.o \
		$(OBJ)/PeriodicRateMatrix.o \
		$(OBJ)/PeriodicSubstMatrix.o
#--------------------------------------------------------
$(OBJ)/SingleGainParsimony.o:\
		SingleGainParsimony.C\
		SingleGainParsimony.H
	$(CC) $(CFLAGS) -o $(OBJ)/SingleGainParsimony.o -c \
		SingleGainParsimony.C
#--------------------------------------------------------
$(OBJ)/RateMatrixDecomposed.o:\
		RateMatrixDecomposed.C\
		RateMatrixDecomposed.H
	$(CC) $(CFLAGS) -o $(OBJ)/RateMatrixDecomposed.o -c \
		RateMatrixDecomposed.C
#---------------------------------------------------------
$(OBJ)/SubstitutionMixtureModel.o:\
		SubstitutionMixtureModel.C\
		SubstitutionMixtureModel.H
	$(CC) $(CFLAGS) -o $(OBJ)/SubstitutionMixtureModel.o -c \
		SubstitutionMixtureModel.C
#--------------------------------------------------------
$(OBJ)/delete-target-gaps.o:\
		delete-target-gaps.C
	$(CC) $(CFLAGS) -o $(OBJ)/delete-target-gaps.o -c \
		delete-target-gaps.C
#---------------------------------------------------------
delete-target-gaps: \
		$(OBJ)/delete-target-gaps.o
	$(CC) $(LDFLAGS) -o delete-target-gaps \
		$(OBJ)/delete-target-gaps.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/make-phase-0.o:\
		make-phase-0.C
	mpiCC $(CFLAGS) -o $(OBJ)/make-phase-0.o -c \
		make-phase-0.C
make-phase-0: \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/Message.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/make-phase-0.o
	mpiCC $(LDFLAGS) -o make-phase-0 \
	$(OBJ)/AlignmentNmerTable.o \
	$(OBJ)/GapRateMatrix.o \
	$(OBJ)/GapRewritingModel.o \
	$(OBJ)/DegenerateDnaMatch.o \
	$(OBJ)/Message.o \
	$(OBJ)/IndelHistory.o \
	$(OBJ)/ContextType.o \
	$(OBJ)/RateMatrixType.o \
	$(OBJ)/Phylogeny.o \
	$(OBJ)/NmerFelsenstein.o\
	$(OBJ)/NmerRateMatrix.o \
	$(OBJ)/NmerSubstMatrix.o \
	$(OBJ)/ACO_Felsenstein.o \
	$(OBJ)/TRCO_Felsenstein.o \
	$(OBJ)/RCO_Felsenstein.o \
	$(OBJ)/LCO_Felsenstein.o \
	$(OBJ)/HOG_Felsenstein.o \
	$(OBJ)/FitchFelsenstein.o \
	$(OBJ)/FitchParsimony.o \
	$(OBJ)/NthOrdFelsenstein.o \
	$(OBJ)/SubstitutionMatrix.o \
	$(OBJ)/JukesCantor.o \
	$(OBJ)/Kimura2Param.o \
	$(OBJ)/HKY.o \
	$(OBJ)/FEL.o \
	$(OBJ)/REV.o \
	$(OBJ)/RateMatrix.o \
	$(OBJ)/NthOrdRateMatrix.o \
	$(OBJ)/NthOrdSubstMatrix.o \
	$(OBJ)/PeriodicRateMatrix.o \
	$(OBJ)/PeriodicSubstMatrix.o \
	$(OBJ)/make-phase-0.o \
	-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
$(OBJ)/gamma.o:\
		gamma.C
	mpiCC $(CFLAGS) -o $(OBJ)/gamma.o -c \
		gamma.C
gamma: \
	$(OBJ)/gamma.o
	mpiCC $(LDFLAGS) -o gamma \
	$(OBJ)/gamma.o \
	-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------------------
$(OBJ)/subset-phylogeny.o:\
		subset-phylogeny.C
	$(CC) $(CFLAGS) -o $(OBJ)/subset-phylogeny.o -c \
		subset-phylogeny.C
#---------------------------------------------------------
subset-phylogeny: \
		$(OBJ)/subset-phylogeny.o
	$(CC) $(LDFLAGS) -o subset-phylogeny \
		$(OBJ)/subset-phylogeny.o \
		-L. -lphylib -lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/get-leaves.o:\
		get-leaves.C
	$(CC) $(CFLAGS) -o $(OBJ)/get-leaves.o -c \
		get-leaves.C
#---------------------------------------------------------
get-leaves: \
		$(OBJ)/get-leaves.o
	$(CC) $(LDFLAGS) -o get-leaves \
		$(OBJ)/get-leaves.o \
		-L. -lphylib -lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/get-phy-distance.o:\
		get-phy-distance.C
	$(CC) $(CFLAGS) -o $(OBJ)/get-phy-distance.o -c \
		get-phy-distance.C
#---------------------------------------------------------
get-phy-distance: \
		$(OBJ)/get-phy-distance.o
	$(CC) $(LDFLAGS) -o get-phy-distance \
		$(OBJ)/get-phy-distance.o \
		-L. -lphylib -lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/phy-to-graph.o:\
		phy-to-graph.C
	$(CC) $(CFLAGS) -o $(OBJ)/phy-to-graph.o -c \
		phy-to-graph.C
#---------------------------------------------------------
phy-to-graph: \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/phy-to-graph.o
	$(CC) $(LDFLAGS) -o phy-to-graph \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/phy-to-graph.o \
		-lgsl -lm -lgslcblas $(LIBS)
#--------------------------------------------------------
$(OBJ)/phy-to-newick.o:\
		phy-to-newick.C
	$(CC) $(CFLAGS) -o $(OBJ)/phy-to-newick.o -c \
		phy-to-newick.C
#---------------------------------------------------------
phy-to-newick: \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/phy-to-newick.o
	$(CC) $(LDFLAGS) -o phy-to-newick \
		$(OBJ)/IndelHistory.o \
		$(OBJ)/GslVector.o \
		$(OBJ)/Phylogeny.o \
		$(OBJ)/GslMatrix.o \
		$(OBJ)/phy-to-newick.o \
		-lgsl -lm -lgslcblas $(LIBS)
#---------------------------------------------

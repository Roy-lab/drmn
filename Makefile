LFLAG = -lgsl -lgslcblas 
SRC = Framework.C SpeciesDistManager.C	SpeciesClusterManager.C	Expert.C Matrix.C GeneExpManager.C MappedOrthogroupReader.C     MappedOrthogroup.C GeneMap.C GeneTreeManager.C GammaManager.C Gamma.C GeneTree.C GeneNameMapper.C NewickReader.C SpeciesFeatureManager.C DRMNPotential.C common/Error.C	   common/EvidenceManager.C  common/Evidence.C  common/VariableManager.C common/Variable.C LeastDirty.C LeastCFGLasso.C LeastLasso.C LeastL21.C GenericLearner.C Task.C

LIBPATH = lib
INCLPATH1 = include
INCLPATH3 = common

CC=g++
CFLAGS = -g --std=c++0x

BIN = learnDRMN

$(BIN): $(SRC)
	$(CC) $(SRC) -I . -I $(INCLPATH1)  -I $(INCLPATH3)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o learnDRMN -g -Werror -Wuninitialized

clean:
	rm -f $(BIN) *~

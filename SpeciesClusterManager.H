#ifndef _SPECIES_CLUSTER_MNGR_
#define _SPECIES_CLUSTER_MNGR_
/*This class reads in a set of clusters for each of the species. Then it maps each gene
 * to its orthogroup. If there are multiple representations of a gene in an orthogroup, then
 * we consider only one of them. 
 * */
#include <map>
#include <string>
#include "CommonTypes.H"

#include "Task.H"
#include "GenericLearner.H"
#include "LeastLasso.H"
#include "LeastL21.H"
#include "LeastDirty.H"
#include "LeastCFGLasso.H"

typedef enum
{
	LEASTL21,
	LEASTDIRTY,
	LEASTFUSED,
	GREEDY
}LearnMode;


using namespace std;
class Expert;
class GeneExpManager;
class MappedOrthogroupReader;
class GammaManager;
class GeneNameMapper;
class SpeciesFeatureManager;
class SpeciesDistManager;
class DRMNPotential;
//class Evidence;
class EvidenceManager;

typedef map<int,Expert*> CLUSTERSET;
typedef map<int,Expert*>::iterator CLUSTERSET_ITER;

class SpeciesClusterManager
{
	public:
	SpeciesClusterManager();
	~SpeciesClusterManager();
	int setOrthogroupReader(MappedOrthogroupReader*);

	int setGammaManager(GammaManager*);
	int setSrcSpecies(const char*);
	int setTrainOGIDs(map<int,int>&);
	int getTestPerformance(const char*,map<int,int>&);
	int clear();
	int readSpeciesData(const char*);
	int readSpeciesData(const char*, vector<string>& species);

	int initExperts();
	int getExtantClusterAssignment(map<int,map<string,int>*>&);
	int getExtantClusterAssignment(map<int,map<string,int>*>&,map<string,int>&);
	
	map<string,CLUSTERSET*>& getExtantSpeciesClusters();

	//Do an EM for EMINT. Could do one step or could do until convergence
	int estimateExpertParameters(const char*);

	
	// DC added -- init clusters per cell type (true) or to source (false; default)?
	int setInitExpertsPerSpecies(bool setPerSpecies);

	int setMaxAssignments(); // Sets maximal assignments and transition matrices between EMINT and DRMN

	// Do an EM for DRMN. (One step or until convergence?)
	int estimateDRMN(const char*);


	// DRMN utilities

	// Lets us retrieve both penalized and unpenalized likelihoods.
	// You can pass in either training or testing data here as OGIDs.
	// Pass in empty unpenalized/penalized; the function will populate them.
	int getDRMNScore_test(map<int,int>& ogids, double& unpenalized, double& penalized); 

	// Dumps expression, feature data, lists of regulators and targets per cluster for a species
	int dumpFeatures(const char* outDir, vector<string>& speciesList);
	int dumpFeaturesPerCluster(const char* outDir,string& specName);
	double getScore();

	int showClusters(const char*);
	int showClusters_Extant(const char*);
	int showClusters_Ancestral(const char*);
	int showMeans(const char* outputDir);
	int showMeans(const char* outputDir,int);
	int dumpAllInferredClusterAssignments(const char* outputDir);
	int dumpAllInferredClusterAssignments(const char* outputDir,int);
	//int dumpAllInferredClusters_Srcwise(const char* outputDir,vector<string>& speciesList);
	int dumpAllInferredClusters_SrcwiseGrouped(const char* outputDir,vector<string>& speciesList);
	int dumpAllInferredClusters_LCA(const char* outputDir,vector<string>& speciesList, string& lcaName);
	int dumpAllInferredClusterGammas(const char* outputDir,vector<string>& speciesList);
	//int dumpAllInferredClusters_Srcwise_Dup(const char* outputDir,vector<string>& speciesList);
	map<string,int>* getGenesForSpecies(string&);
	int setRandom(bool);
	int setRandSeed(int);
	int setMaxClusterCnt(int);
	int generateData(const char*,string&,vector<string>&); //deprecated - for Arboretum - ont sure if works here
	int setConstCov(double val); // Sets constant cov value (DC Add)

	

	int setRestrictedList(const char*); // Reads regulator OGIDs (from MRTLE)

	/*
	* Prints current data (cluster assignments, conditionals, etc) to a directory.
	*/
	int printCurrentData(const char* outputDir); // DO NOT USE

	/*
 	* Prints the DRMN results for training data, for each cell type.
 	* - Network (cell_net.txt)
 	* - Predicted expression (as mean) (cell_pred.txt, cell_pred_maxclust_exprtab.txt)
 	* - Sampled expression (sample from gaussian instead of using mean) (cell_sample_maxclust_exprtab.txt)
 	*/
	int showDRMNResults(const char*);
	int setMode(LearnMode lm, double r1, double r2, double r3);
	
	private:

	
	int readClusters(string&, const char*);
	int maximizationStep();
	int expectationStep();
	int expectationStep_Species(string& specName, CLUSTERSET* expertSet);

	
	int expectationStep_DRMN(); // Expectation step for DRMN: update gammas
	int expectationStep_DRMN_Species(string& specName, CLUSTERSET* expertSet);
	int doInferenceForGene(string& geneName,string& speciesName);
	//int maximizationStep1_DRMN(); // Maximization step 1 for DRMN (Transition probabilities)
	//int maximizationStep2_DRMN(); // Maximization step 2 for DRMN (regulators)
	double getDRMNScore();
	double getDRMNScoreForModule(int);

	int precomputeEmptyGraphPrior();
	int estimateRegProgs();
	int estimateRegProgs_PerModule(int moduleID);
	double assessScoreImprovement(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential** potPtr);
	//double assessScoreImprovement(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential* potPtr);
	double getScoreForPot(map<string,int>& geneSet,DRMNPotential* apot, string& cellType);
	double getScoreForPot_Tracetrick(map<string,int>& geneSet,DRMNPotential* apot, DRMNPotential* ppot, string& cellType);

	/* Ali */
	double assessScoreImprovement_Bookkeeping_Ali(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential** potPtr, map<int,INTDBLMAP*> baseCovar, Matrix* X, Matrix* Y);
	/* DC */
	double assessScoreImprovement_Bookkeeping(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential** potPtr, map<int,INTDBLMAP*> baseCovar);
	int copyCovar(map<int,INTDBLMAP*>& source, map<int,INTDBLMAP*>& target);
	int printCovar(map<int,INTDBLMAP*>& gCovar);
	int computeBaseUnnormCovar(Expert* e,string& cellType, map<int,INTDBLMAP*>& baseCovar);


	int estimateMeanCov(Expert* e, string& specName, int clusterID);
	int displaySpeciesClusters(const char* outFName,string& specName);
	int displaySpeciesClusterAssignments(const char* outFName,string& specName);
	int displaySpeciesClusterAssignments_NonSrcNames(const char* outFName,string& specName);
	int assignGenesToExperts();
	int assignGenesToExperts_FromMap();
	int showDRMNModelForCellType(const char* pathName,map<int,Expert*>* expertSet);

	// For each training gene, prints its true expression, cluster assignment, 
	// and predicted expression from each module. creates the cell_pred.txt file
	int showDRMNPredictionsForCellType(const char* pathName,map<int,Expert*>* expertSet,string& cellType);

	// For each training gene, prints just the predicted expression from the best cluster assignment
	// using its regulatory feature data.
	// Creates the cell_pred_maxclust_exprtab.txt file
	int showBestDRMNPredictionsForCellType(const char* pathName,map<int,Expert*>* expertSet,string& cellType);
	
	// For each training gene, samples an expression value from the best cluster assignment
	// (using the gene's own input regulatory feature data)
	// Creates the cell_sample_maxclust_exprtab.txt file.
	int showBestDRMNSamplesForCellType(const char* pathName, map<int,Expert*>* expertSet, string& cellType);

	/*
 	* Prints OGIDs file to feed into MRTLE step. cid is cluster ID to disambiguate in output files.
	* Eventually we want to replace this with direct communication
 	*/
	int makeOGIDsForRegnet(const char* ofname, vector<string>& speciesList);
	
	/*
	* Chooses an OGID to use for a module -- TO DO
	*/
	int setUpModuleOGIDs();

	/*
 	* Prints regulator and target files for all species for one cluster
 	*/
	//int printRegTargetFiles(const char* ofPrefix, int cid, vector<string>& speciesList);
	
	map<string,CLUSTERSET*> speciesExpertSet;
	map<string,map<string,int>*> speciesClusterSet_Genewise;

	//Key is the ogid and the subsequent key is the name of the species.
	map<int,map<string,int>*> mappedClusterAssignment; // We use these to hold onto current cluster assignment

	//First key is the species name. Second level key is the ogid
	map<int,map<string,map<int,double>*>*> mappedClusterGamma;
	map<string,map<int,int>*> speciesClusterSet_OGwise;
	map<string,int> testGenes;

	map<string,GeneExpManager*> speciesExprSet; // Stores gene expression

	// Stores other input cis-regulatory data per species
	// like in the expression manager, targets are stored as strings
	map<string,SpeciesFeatureManager*> speciesFeatSet; 


	MappedOrthogroupReader* mor; // stores orthology info for genes
	map<int,int> inputRegulatorOGs; // input regulator OGIds
	int moduleOG; // special OGID to be used for the "target" 
	
	GammaManager* gammaMgr;
	int maxClusterCnt;
	bool randFlag;
	
	GeneNameMapper gnm;
	char srcSpecies[256];
	
	gsl_rng* r; // RNG 
	int rseed; 
	map<int,int> initRandAssign;

	double constCov; // Fixed value for covariance (DC)
	bool fixCov; // If true, uses fixed value (DC)
	map<string,int> allRegulatorSet;
	map<string,int> allRegulatorIndex;
	map<string,int> varNameIDMap;
	map<int,string> varIDNameMap;
	gsl_matrix* ludecomp;
	gsl_permutation* perm;
	//This will store the status of edge across cell types
	map<string,STRINTMAP*> edgeCelltypeMap;
	//This will store and update how the prior changes as we add edges into each cell type's modules
	map<string,double> regModulePairSetPrior;
	map<string,EvidenceManager*> evMgrSet;
	map<int,int> trainSetOGIDs;
	map<int,int> testSetOGIDs;  
	
	// We use this to keep track of any OGIDs that were specified in the cluster assignment for each species
	map<string, map<int,int>*> legalTestOGIDs; 

	int maxEMINTIter; // iterations of EMINt for initialization
	int maxDRMNIter; // iterations of the full DRMN loop
	int regsPerIter; // iterations of regulator adding

	// DC added -- init clusters per cell type (true) or to source (false; default)?
	bool initExpertsPerSpecies;

	int estimateCov(Expert* e, Matrix* X, Matrix* Y, map<int,INTDBLMAP*>& gCovar);
	int estimateCov_All(Expert* e, Matrix* X, Matrix* Y, map<int,INTDBLMAP*>& gCovar);
	int addCov(Expert* e, Matrix* X, Matrix* Y, map<int,INTDBLMAP*>& gCovar, string regName);

	int estimateRegProgs_PerModule_LASSO(int moduleID);
	double makeOnePotential(Expert* e,vector<string>& regNames,string& cellType,bool& regStatus,DRMNPotential** potPtr, Matrix* X, Matrix* Y);
	int learnLEAST(vector<Task_T*>* allt,vector<Matrix*>& allW);
	int countAllW(vector<Matrix*> allW);
	int clearAllW(vector<Matrix*>& allW);

	LearnMode learnMode;
	double rho1;
	double rho2;
	double rho3;
};

#endif

#ifndef _SPECIES_DATA_
#define _SPECIES_DATA_

#include <map>
#include <vector>
#include <string>
using namespace std;
class Matrix;
class SpeciesDistManager
{
	public:
	SpeciesDistManager();
	~SpeciesDistManager();
	int readSpeciesTree(const char*);
	int setMaxClusters(int k);
	double getAncestralClustering(map<string,int>& extantStatus,map<string,int>& ancestralStatus);
	double getExtantClustering(map<string,int>& ancestralClustering, map<string,int>& extantClustering);
	double getConditionalProb(string& child,int parentCluster,int childCluster);
	double getConditionalProbForEdge(string& child,int parentEdgeStatus,int childEdgeStatus,int moduleID);
	//double getProbGainLoss(string&);
	//double getProbLossGain(string&);
	//double getProbLossLoss(string&);
	//int enumerateAllAncestors(map<string,int>& edgeStatus);
	double scoreAssignment(map<string,int>& jointAssign);
	double getEdgeStatusProb(map<string,int>& edgeStatus,int);
	int showInferredConditionals(const char* outputDir);
	int showInferredConditionals_ML(const char* outputDir);
	
	// Show tree restricted to a cluster -- like required for MRTLE
	int showTreeForCluster(int clustID);	
	int printTreeForCluster(int clustID, const char* outFN);	
	

	struct Species;

	struct Species 
	{
		string name;
		Species* parent;
		vector<Species*> children;
		Matrix* conditional;
		Matrix* conditional_ml;
		inline Matrix* getParams()
		{
			return conditional;
		}
		map<int,vector<int>*> sortedClustIDs;
		inline vector<int>* getSortedClusterIDs(int rowID)
		{
			vector<int>* sortedIDs=NULL;
			if(sortedClustIDs.find(rowID)!=sortedClustIDs.end())
			{
				sortedIDs=sortedClustIDs[rowID];
			}
			return sortedIDs;
		}
		inline int setSortedClusterIDs(int rowId,vector<int>* sortedIDs)
		{
			sortedClustIDs[rowId]=sortedIDs;
		}
	};
	
	Species* getRoot();
	int assignLevel();
	int assignLevel(SpeciesDistManager::Species*,int);
	int getLevelFromRoot(const char* nodeName);
	int getSpeciesListPrefix(vector<string>&);
	map<string,Species*>& getAllSpecies();
	Species* getSpecies(string&);
	Matrix* getConditional(string& spName);
	int resetTransitionProbability();
	int resetTransitionProbability(Species*);
	int normalizeTransitionMatrix();
	int normalizeTransitionMatrix(Species*);

	int initTransitionMatrix_ML();
	int initTransitionMatrix_ML(Species*);
	Matrix* getTransitionMatrix_ML(string& name);
	int normalizeTransitionMatrix_ML();
	int normalizeTransitionMatrix_ML(Species*);
	
	int assignExtantClustering(map<string,int>& ancestralClusterAssign,Species* node, map<string,int>& extantClusterAssign);
	double getSubTree(int parentEdgeStatus,Species*,map<string,int>&, int moduleID);
	double maxSubTree(int parentCluster, Species* child, map<string,int>& extantCluster,map<string,int>& ancestralCluster);
	//int enumerateChild(Species*,map<string,int>&,map<string,int>&,vector<map<string,int>*>&);
	int getMaxClusterAssignForChild(int parentClustID,Species* child);
	int showConditionals(const char*,SpeciesDistManager::Species* species);
	int showConditionals_ML(const char*,SpeciesDistManager::Species* species);
	//Should be called on a per datapoint level, where a datapoint represents the joint assignment of all
	//variables. The input parameter, dataPtProb is the probability density of the data point from the
	//ith cluster, where i is the key of the map.
	//The output parameter is the normalization constant for each node.
	//int estimateNormalizationConstants(map<string,map<int,double>*>& dataPtProb, map<string,double>& normConstants);
	//int estimateNormalizationConstant_Node(map<string,map<int,double>*>& dataPtProb, map<string,double>& normConstants,Species*,int ancestralClustID);
	
	private:
	
	// Prints tree for cluster to stdout (helper function)
	int showTreeForCluster(Species* node, int clustID, int childID);
	// Prints tree for cluster to a file so we can read it in for MRTLE (helper function)
	int printTreeForCluster(Species* node, int clustID, int childID, ofstream& oFile);

	int getSpeciesListPrefix(Species*,vector<string>&);
	map<string,Species*> speciesSet;
	Species* root;
	int maxClusterCnt;
	map<string,int> levelFromRoot;
};
#endif

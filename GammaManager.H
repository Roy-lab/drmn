#ifndef _GAMMA_MGR
#define _GAMMA_MGR
#include <map>
#include <string>
#include <unistd.h>
#include <gsl/gsl_randist.h>
using namespace std;

class SpeciesDistManager;
class Gamma;
class GeneTreeManager;
class MappedOrthogroupReader;

class GammaManager
{
	public:
	GammaManager();
	~GammaManager();
	int setMaxClusterCnt(int k);
	int setOrthogroupReader(MappedOrthogroupReader* aPtr);
	int setSpeciesDistManager(SpeciesDistManager* aPtr);
	SpeciesDistManager* getSpeciesDistManager(); // DC add

	int initGamma(int ogid, string& geneName, string& specName,int clustID);
	Gamma* getGammaForOGID(int);
	int setObservedProbs(int,map<int,double>,string& geneName, string& specName);
	int estimateLeafGamma(int,map<int,double>&,string& geneName, string& specName);
	int estimateLeafGamma(Gamma::Node* m);
	int estimateNonLeafPosterior();
	int estimateNonLeafPosterior_DRMN();


	double getNonLeafScore();
	Matrix* getGamma(int ogid,string& geneName, string& specName);
	Matrix* getNormTerm(int ogid,string& geneName, string& specName);
	double getPrior(string&,int,int);
	double getAllNodeScore();
	map<int,double>* getLeafLikelihood_store(int ogid,string& geneName);
	double getNodeScore(Gamma::Node*);
	int estimateNonLeafPosterior(Gamma::Node* node);
	int estimateNonLeafPosterior_DRMN(Gamma::Node* node);

	// DC added this to get scores for choice nodes
	// key is OGID; value is just for participation
	double getAllNodeScore(map<int,int>& ogids);

	int getAllClusterAssignments(map<int,map<string,int>*>& allClusterAssignments,bool);
	int getAllClusterAssignments_Grouped(map<int,map<string,int>*>& allClusterAssignments);
	int getAllClusterAssignments_Conditional(map<int,map<string,int>*>& allClusterAssignments);
	int revisitLeafAssignments(Gamma::Node* anode, map<string,int>& clusterAssignment_prior, map<string,int>* clusterAssignment_posterior);
	int getAllClusterGammas(map<int,map<string,map<int,double>*>*>& allClusterGammas);
	int getMaximalAssignment(Gamma::Node* anode,map<string,int>* assignment);
	int getMaximalAssignment_Grouped(Gamma::Node* anode,map<string,int>* assignment,int dupID);

	//This remains 1 until we reach a duplication node which is 2
	int getMaximalAssignment(Gamma::Node* anode,map<string,int>* assignment,int dupID);
	int getMaximalAssignment(Gamma::Node*,map<string,int>& assignment);
	int getMaximalConditionalAssignment(Gamma::Node* anode,map<string,int>* assignment,int dupID,int parentClusterID);
	int getClusterProb(Gamma::Node* anode,map<string,map<int,double>*>* prob);
	int estimateTransitionProbability();

	int estimateTransitionProbability_DRMN(); // DC

	int reestimateTransitionProbability();
	int showClusterFlipCnts();
	const char* getDupAncestor(int);
	inline int getTotalUpdatedParentCnt()
	{
		return gtMgr.getTotalAdditionalParents();
	}
	int showTotalUpdatedParents();
	int sampleAllClusterAssignments(map<int,map<string,int>*>& allClusterAssignments);
	
	int clearGammaSet(); // DC added to clear out gamma set -- is this enough?

	private:


	int updateTransitionMatrix(Gamma::Node*);
	int updateTransitionMatrix_ML(Gamma::Node*,int);
	int updateTransitionMatrix_GammaFix(Gamma::Node*);
	bool hasDuplicate(Gamma::Node*);

	
	int drawAssignment(Gamma::Node* anode, map<string,int>* assignment,int dupID,int parentID);
	int sampleAncestralCluster(gsl_rng* r,SpeciesDistManager::Species* root);
	int sampleChildCluster(gsl_rng* r,SpeciesDistManager::Species* child,int parentID);
	int sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>* sortedClustIDs);
	int sortIndices(Matrix* params,int parentID,vector<int>* sortedClustIDs);

	SpeciesDistManager* spdistMgr;
	//We will store one instance of Gamma for each OGid
	map<int,Gamma*> gammaSet;
	int maxClusterCnt;
	GeneTreeManager gtMgr;
	MappedOrthogroupReader* mgrpReader;
	bool showGammas;
	map<int,string> ogDupAncMap;
	map<string,Matrix*> clusterFlipCnt;
	gsl_rng* r;
};

#endif

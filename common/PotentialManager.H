/* The goal of this class is to take all the MarkovBlankets and the factors and
 *  the data and create all potentials involving valid pairings MarkovBlankets 
 *  and factors. We create all putative potentials upfront because we want to
 *  iterate over the data once and update counts of configurations in one go*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include <vector>
#include <map>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_permutation.h>

using namespace std;

class Potential;
class EvidenceManager;
class SlimFactor;

//This stores all the potentials
typedef map<int,Potential*> POTFUNC_MAP;
typedef map<int,Potential*>::iterator POTFUNC_MAP_ITER;


class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();
		int setEvidenceManager(EvidenceManager*);
		int setRestrictedNeighborSet(map<int,Variable*>&);
		int setOutputDir(const char*);
		int estimateAllMeanCov(bool,INTDBLMAP&,map<int,INTDBLMAP*>&,INTINTMAP&,const char*,const char*,int leaveOutData=-1);
		int estimateAllMeanCov(bool,INTDBLMAP&,map<int,INTDBLMAP*>&,INTINTMAP&);
		int estimateCovariance(bool random,INTDBLMAP* vcov, int vId, int uId);
		int readAllMeanCov(const char*, const char*);
		int initPooled();
		int init();
		int init(int);
		int initValidationSet(int);
		int initValidationSet(int,int);
		int reset();
		int initRandom();
		
		Error::ErrorCode populatePotentialsSlimFactors(map<int,SlimFactor*>&,VSET&);
		Error::ErrorCode estimateRandomInfo(map<int,SlimFactor*>&,VSET&,vector<double>&,int fSize);
		int estimateMarginalEntropies(map<int,SlimFactor*>&,VSET&,bool);
		int populateFactor(map<int,SlimFactor*>&,VSET&,SlimFactor*,bool random);
		int populateFactor_Buffer(map<int,SlimFactor*>&,VSET&,SlimFactor*,bool random);
		//We are going to store the set of potentials per factor. That is each 
		//POTVECT object stores the set of potentials for a particular factor
		int populatePotential(Potential*,bool random);
		double getPseudoLikelihood(SlimFactor*,VSET&,bool);
		double getPseudoLikelihood(SlimFactor* sFactor,VSET& varSet, bool train,int& status);
		double getPseudoLikelihood(SlimFactor* sFactor,VSET& varSet, bool train,int& status, Potential*);
		double getGaussianLikelihood(map<int,SlimFactor*>&,VSET&,bool);

		double getLikelihood(SlimFactor*,VSET&);
		double getLikelihood(SlimFactor*,VSET&,map<int,int>&);
		int estimateConditionalPotential(SlimFactor*,VSET&,Potential**,STRDBLMAP&);
		int populatePotential(Potential*,STRDBLMAP&);
		int estimateCanonicalPotential(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		int estimateCanonicalPotential_Abbeel(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		int estimateCanonicalPotential_Approximate(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		int resetPotFuncs();
		int estimateCanonicalPotential_Joint(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		Potential* getPotential(int);
		double getConditionalEntropy(int,INTINTMAP&,VSET&);

		double getSampleLikelihood(map<int,SlimFactor*>&, VSET&, INTINTMAP*);
		int getVariableSample(INTINTMAP&,VSET&,int,SlimFactor*,gsl_rng*);

		int clearJointEntropies();

	private:
		int getInstantiation(const string& confString, INTINTMAP& confMap);
		double estimateCanonicalValue(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* condPot);
		double estimateCanonicalValue_Joint(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* condPot);
		POTFUNC_MAP potFuncs;
		//This stores the joint potential for a joint configuration and then can
		//be reused over and over again to generate the conditional
		map<string,Potential*> canonicalPotFunc;
		EvidenceManager* evMgr;
		map<string,double> summaryStat;
		map<string,double> jointEntropies;
		INTDBLMAP globalMean;
		INTDBLMAP globalMean_Rand;
		map<int,INTDBLMAP*> globalCovar;
		map<int,INTDBLMAP*> globalCovar_Rand;
		map<int,Variable*> restrictedNeighborSet;
		char outputDir[1024];
		gsl_matrix* ludecomp;
		gsl_permutation* perm;

		map<int,Potential*> potBuffer;
};
#endif

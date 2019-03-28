#ifndef _DRMN_POTENTIAL
#define _DRMN_POTENTIAL
//A potential over continuous variables needs only store the mean and variance parameters.
//However to compute the pdf, we need to store the inverse of the covariance matrix and the
//determinant.
#include <fstream>
#define SMALL_PROB 1e-60
#define PI 3.1472
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#include "Matrix.H"
//class Evidence;
class DRMNPotential
{
	public:
		DRMNPotential();
		~DRMNPotential();
		
		int resetVarSet();
		int initMatrix(int);
		//This enum is used to decide what role the Variable
		//plays in the potential. A Variable can be either from a Markov Blanket
		//or from the factor. This is required while converting the joint into a conditional
		//of the factor variables give Markov blanket variables
		
		typedef enum
		{
			FACTOR,
			MARKOV_BNKT
		}VariableRole;
		
		int setAssocVariable(int,VariableRole);
		map<int,int>& getAssocVariables();
		int updateFactorVariable(int);

		//Initialize potential 
		int potZeroInit();
		int potZeroInit_MeanOnly();
		int potCustomInit(double,double);
		double getJointPotValueFor(INTDBLMAP&);
		//double getCondPotValueFor(int vId, INTDBLMAP&);
		double getJointPotValueFor(STRINTMAP&);
		double getJointPotValueForConf(string&);
		
		//The result of marginalization is put in to the second
		//parameter
		//Needed when we are updating the joint potential
		int updateMean(int,double);
		int updateCovariance(int,int,double);
		int dumpPotential(ostream&);
		int showVariables(ostream&);
		int makeValidJPD();
		int makeValidJPD(gsl_matrix*,gsl_permutation*);

		//Since both the joint and the conditional are in the same potential 
		//we can easily compute the entropy
		int calculateEntropy();
		int calculateJointEntropy();
		double getEntropy();
		double getJointEntropy();
		// not using this one
		double generateSample(INTDBLMAP&, int vId,gsl_rng*);

		int initMBCovMean();
		double getCondPotValueFor(INTDBLMAP&);
		//double getCondPotValueFor(map<int,Evidence*>*);
		double getCondPotValueFor(map<int,double>*);
		double computeLL_Tracetrick(int sampleSize);
		double getCondVariance();
		double getCondBias();
		INTDBLMAP& getCondWeight();
		int setParameters(double mean, double var, double coeffs);

		// predicts sample as mean of multivariate gaussian given markov blanket vals
		//double predictSample(map<int,Evidence*>* evidenceSet);
		double predictSample(map<int,double>* evidenceSet);
		// draws sample from multivariate gaussian given markov blanket vals
		//double generateSample(map<int,Evidence*>* evidenceSet, gsl_rng* r);
		double generateSample(map<int,double>* evidenceSet, gsl_rng* r);


	protected:
		INTINTMAP varSet;
		INTINTMAP factorVariables;
		INTINTMAP markovBlnktVariables;
		INTDBLMAP meanWrap;
		Matrix* covariance;
		Matrix* inverse;
		Matrix* mean;
		double determinant;
		double normFactor;
		double potEntropy;
		double jointEntropy;
		//Each variable is mapped to an index in the covariance or mean matrix
		//This map stores that information
		int matrixInd;
		INTINTMAP vIDMatIndMap;
		INTINTMAP matIndvIDMap;

		double mbcondVar;
		double mbcondMean_Part;
		INTDBLMAP mbcondMean_Vect;
};
#endif

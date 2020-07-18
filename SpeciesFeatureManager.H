#ifndef _SPECFEATMGR_
#define _SPECFEATMGR_

#include <map>
#include <string>
#include "CommonTypes.H"
#include <unistd.h>

using namespace std;

/*
* Storing the cis-regulatory feature data for a species (cell type).
* DC Jan 2017
*/


class SpeciesFeatureManager
{
	public:
	SpeciesFeatureManager();
	~SpeciesFeatureManager();
	int readFeatures(const char*);
	int readFeatures_Efficient(const char*);
	map<string,map<string,double>*>& getFeatures();

	//SR: This will return the number of targets that have an overlap with the input gene set. Let the SpeciesFeatureManager take care of this.
	int getTargetHitCntForRegulator(string&,map<string,int>&);
	
	// checks for presence of regulator-target value
	bool hasValue(string regulator, string target);

	// checks for presence of regulator
	bool hasRegulator(string regulator);
	// checks for presence of target
	bool hasTarget(string target);

	// Gets the value associated with a regulator and target
	double getValue(string regulator, string target);

	// Given empty regulator map, retrieves regulators in this species
	int getRegulators(map<string,int>& regulators);

	// Default value for reg-target edge if not present in data.
	double getDefaultValue();

	// Get expected fraction of targets
	double getExpectedTargetCnt(string&);
	//SR: Get the targets associated with a regulator
	map<string,double>* getTargetsForRegulator(string&);
	map<string,double>* getRegulatorsForTarget(string&);
	double* getRegulatorsForTarget_Efficient(string&);
	map<string,int>& getFeatureNames();
	private:
	map<string,map<string,double>*> motifNetwork;
	map<string,map<string,double>*> tgtTFNetwork;
	map<string,int> allTargets;
	map<string,int> regulatorNameIDMap; 
	
	map<string,int> tfNameIDMap;
	map<string,int> tgtNameIDMap;
	double DEFAULT_VALUE;
	double** motifNetwork_TFcentric;
	double** motifNetwork_Tgtcentric;

};
#endif

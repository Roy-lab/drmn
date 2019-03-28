#include "SpeciesFeatureManager.H"
#include <string.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>

/*
* Storing the cis-regulatory feature Feature for a species (cell type).
* DC Jan 2017
*/
SpeciesFeatureManager::SpeciesFeatureManager()
{
	DEFAULT_VALUE=0.0;
}

SpeciesFeatureManager::~SpeciesFeatureManager()
{
	allTargets.clear();
	regulatorNameIDMap.clear();
	for(map<string,map<string,double>*>::iterator ritr=motifNetwork.begin(); ritr!=motifNetwork.end(); ritr++)
	{
		map<string,double>* tmap = ritr->second;
		tmap->clear();
		delete tmap;
	}
	motifNetwork.clear();
}

/*
* Reads network from 3-column file:
* cis_reg_element \t target \t score
*/
int
SpeciesFeatureManager::readFeatures(const char* aPtr)
{
	ifstream inFile(aPtr);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1024);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string tfName;
		string tgtName;
		double score=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfName.append(tok);
			}
			else if(tokCnt==1)
			{
				tgtName.append(tok);
			}	
			else if(tokCnt==2)
			{
				score=atof(tok);				
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
	
		map<string,double>* tgtSet=NULL;
		if(motifNetwork.find(tfName)==motifNetwork.end())
		{
			tgtSet=new map<string,double>;
			motifNetwork[tfName]=tgtSet;
		}
		else
		{
			tgtSet=motifNetwork[tfName];
		}
		(*tgtSet)[tgtName]=score;
		allTargets[tgtName]=1; // Store target in all-targets
	}
	inFile.close();
	return 0;
}

/*
* Returns the feature data (motif network)
*/
map<string,map<string,double>*>& 
SpeciesFeatureManager::getFeatures() 
{
	return motifNetwork;
}


int 
SpeciesFeatureManager::getTargetHitCntForRegulator(string& regulatorName,map<string,int>& ageneSet)
{
	if (!hasRegulator(regulatorName))
	{
		return 0;
	}

	map<string,double>* targetSet=motifNetwork[regulatorName];
	int hitCnt=0;
	for(map<string,int>::iterator gIter=ageneSet.begin();gIter!=ageneSet.end();gIter++)
	{
		if(targetSet->find(gIter->first)!=targetSet->end())
		{
			hitCnt++;		
		}
	}
	return hitCnt;
}

/*
* Populates the provided regulator map (assumed empty...) with the regulators in this network.
* I guess we could also use this to get the union of regulators over multiple motif nets, so we don't check.
*/
int 
SpeciesFeatureManager::getRegulators(map<string,int>& regulators)
{
	for(map<string,map<string,double>*>::iterator gIter=motifNetwork.begin();gIter!=motifNetwork.end();gIter++)
	{
		regulators[gIter->first]=1;
	}
	return 0;
}

/*
* Checks for presence of regulator, target 
*/
bool
SpeciesFeatureManager::hasValue(string regulator, string target)
{
	if (!hasRegulator(regulator))
	{
		return false;
	}		
	map<string,double>* tgtSet=motifNetwork[regulator];
	return (tgtSet->find(target)!=tgtSet->end());
}

bool
SpeciesFeatureManager::hasRegulator(string regulator)
{
	return (motifNetwork.find(regulator)!=motifNetwork.end());
}

bool
SpeciesFeatureManager::hasTarget(string target)
{
	return (allTargets.find(target)!=allTargets.end());
}

/*
* Gets the value associated with a regulator and target.
*/
double 
SpeciesFeatureManager::getValue(string regulator, string target)
{
	if (!hasRegulator(regulator))
	{
		return false;
	}		
	map<string,double>* tgtSet=motifNetwork[regulator];
	return (*tgtSet)[target];
}

/*
* Gets default value for edge not present in data (eg 0).
*/
double
SpeciesFeatureManager::getDefaultValue()
{
	return DEFAULT_VALUE;
}	


double 
SpeciesFeatureManager::getExpectedTargetCnt(string& regName)
{
	if (!hasRegulator(regName))
	{
		return 0;
	}
	map<string,double>* tgtSet=motifNetwork[regName];
	double frac=((double)tgtSet->size())/((double)allTargets.size());
	return frac;
}

map<string,double>*
SpeciesFeatureManager::getTargetsForRegulator(string& regName)
{
	if (!hasRegulator(regName))
	{
		return NULL;
	}
	map<string,double>* tgtSet=motifNetwork[regName];
	return tgtSet;
}




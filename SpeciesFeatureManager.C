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
		map<string,double>* tfSet=NULL;
		if(tgtTFNetwork.find(tgtName)==tgtTFNetwork.end())
		{
			tfSet=new map<string,double>;
			tgtTFNetwork[tgtName]=tfSet;
		}
		else
		{
			tfSet=tgtTFNetwork[tgtName];
		}
		(*tfSet)[tfName]=score;
		allTargets[tgtName]=1; // Store target in all-targets
	}
	inFile.close();
	return 0;
}


int
SpeciesFeatureManager::readFeatures_Efficient(const char* aPtr)
{
	ifstream inFile(aPtr);
	char buffer[1024];
	int lineCnt=0;
	double score=0;
	int tfCnt=0;
	int tgtCnt=0;
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
		if(lineCnt==0)
		{
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					tfCnt=atoi(tok);
				}
				else if(tokCnt==1)
				{
					tgtCnt=atoi(tok);
				}	
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
			motifNetwork_TFcentric=new double*[tfCnt];
			motifNetwork_Tgtcentric=new double* [tgtCnt];
		}
		else
		{
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
			int tfID=0;
			if(tfNameIDMap.find(tfName)==tfNameIDMap.end())
			{
				tfID=tfNameIDMap.size();
				tfNameIDMap[tfName]=tfID;
				if(tfID>=tfCnt)
				{
					cout<<"Oops! Too many TFs " << tfID << " than we requested memory for: " << tfCnt << endl;
					exit(0);
				}
				motifNetwork_TFcentric[tfID]=new double[tgtCnt];
				for (int ti=0;ti<tgtCnt;ti++)
				{
					motifNetwork_TFcentric[tfID][ti]=0;
				}
			}
			else
			{
				tfID=tfNameIDMap[tfName];
			}
			int tgtID=0;
			if(tgtNameIDMap.find(tgtName)==tgtNameIDMap.end())
			{
				tgtID=tgtNameIDMap.size();
				tgtNameIDMap[tgtName]=tgtID;
				if(tgtID>=tgtCnt)
				{
					cout<<"Oops! Too many Targetss " << tgtID << " than we requested memory for: " << tgtCnt << endl;
					exit(0);
				}
				motifNetwork_Tgtcentric[tgtID]=new double[tfCnt];
				for (int ti=0;ti<tfCnt;ti++)
				{
					motifNetwork_Tgtcentric[tgtID][ti]=0;
				}
			}	
			else
			{
				tgtID=tgtNameIDMap[tgtName];
			}
			motifNetwork_Tgtcentric[tgtID][tfID]=score;
			motifNetwork_TFcentric[tfID][tgtID]=score;
		}
		/*map<string,double>* tgtSet=NULL;
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
		map<string,double>* tfSet=NULL;
		if(tgtTFNetwork.find(tgtName)==tgtTFNetwork.end())
		{
			tfSet=new map<string,double>;
			tgtTFNetwork[tgtName]=tfSet;
		}
		else
		{
			tfSet=tgtTFNetwork[tgtName];
		}
		(*tfSet)[tfName]=score;
		allTargets[tgtName]=1; // Store target in all-targets*/
		
		lineCnt++;
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


map<string,double>* 
SpeciesFeatureManager::getRegulatorsForTarget(string& tgtName)
{
	if(tgtTFNetwork.find(tgtName)==tgtTFNetwork.end())
	{
		return NULL;
	}
	map<string,double>* tfSet=tgtTFNetwork[tgtName];
	return tfSet;
}
double* 
SpeciesFeatureManager::getRegulatorsForTarget_Efficient(string& tgtName)
{
	if(tgtNameIDMap.find(tgtName)==tgtNameIDMap.end())
	{
		return NULL;
	}
	int tgtID=tgtNameIDMap[tgtName];
	double* tfs=motifNetwork_Tgtcentric[tgtID];
	return tfs;
}


map<string,int>& 
SpeciesFeatureManager::getFeatureNames()
{
	return tfNameIDMap;
}

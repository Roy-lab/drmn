#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
//#include "Evidence.H"
#include "EvidenceManager.H"


EvidenceManager::EvidenceManager()
{
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=-1;
}

EvidenceManager::~EvidenceManager()
{
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

/*
Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Continuous(const char* inFName)
{
	ifstream inFile(inFName);
	char buffer[160000];
	int lineNo=0;
	while(inFile.good())
	{
		inFile.getline(buffer,160000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		else if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		if(lineNo>=500)
		{
			lineNo++;
		//	continue;
		}

		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value
		while(tok!=NULL)
		{
			Evidence* evid;
			if(populateEvidence_Continuous(&evid,tok)==-1)
			{
				cout <<"Error while populating evidence " << endl;
				return Error::DATAFILE_ERR;
			}
			(*evidMap)[evid->getAssocVariable()]=evid;
			tok=strtok(NULL,"\t");
		}
		//if(evidenceSet.size()<=20)
		//{
			evidenceSet.push_back(evidMap);
		//}
	}
	inFile.close();
	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;
	//normalize();
	return Error::SUCCESS;
}
*/

/*
int
EvidenceManager::normalize()
{	
	for(int e=0;e<evidenceSet.size();e++)
	{
		EMAP* evidMap=evidenceSet[e];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			Evidence* evid=eIter->second;
			double aVal=evid->getEvidVal();
			double maxVal=maxAbsValues[eIter->first];
			aVal=aVal/maxVal;
			evid->setEvidVal(aVal);
		}
	}
	return 0;
}
*/

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	VSET& variableSet=vMgr->getVariableSet();
	int* randInds=new int[evidenceSet.size()];
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,randEvidenceSet.size());	
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[randInds[i]];
			EMAP* randEvidMap=randEvidenceSet[i];
			//Evidence* evid=(*evidMap)[vIter->first];
			//(*randEvidMap)[vIter->first]=evid;
			double evid=(*evidMap)[vIter->first];
			(*randEvidMap)[vIter->first]=evid;
		}
	}
	return 0;
}

int 
EvidenceManager::getNumberOfEvidences()
{
	return evidenceSet.size();
}

EMAP* 
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}


EMAP* 
EvidenceManager::getEvidenceAt(string& name)
{
	if(sampleNameIDs.find(name)==sampleNameIDs.end())
	{
		return NULL;
	}
	int ID=sampleNameIDs[name];
	return evidenceSet[ID];
}


EMAP* 
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}


int
EvidenceManager::addEvidence(EMAP* evidSet)
{
	evidenceSet.push_back(evidSet);
	return 0;
}

//A better to handle this might have been to store the name as another attribute!
int
EvidenceManager::addEvidenceWithName(EMAP* evidSet,string& name)
{
	int id=evidenceSet.size();
	sampleIDNames[id]=name;
	sampleNameIDs[name]=id;
	evidenceSet.push_back(evidSet);
	return 0;
}

int 
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int
EvidenceManager::generateValidationSet(const char* vFName, int vSetSize,gsl_rng* r)
{
	ifstream inFile(vFName);
	if(inFile.good())
	{
		char buffer[256];
		while(inFile.good())
		{
			inFile.getline(buffer,255);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			int dId=atoi(buffer);
			validationIndex[dId]=0;
		}
		inFile.close();
	}
	else
	{
		populateRandIntegers(r,validationIndex,evidenceSet.size(),vSetSize);
		ofstream oFile(vFName);
		for(INTINTMAP_ITER vIter=validationIndex.begin();vIter!=validationIndex.end();vIter++)
		{
			oFile << vIter->first << endl;
		}
		oFile.close();
	}
	return 0;
}

int 
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int 
EvidenceManager::setPreRandomizeSplitSeed(int seed)
{
	randseed=seed;
	return 0;
}

/*
int 
EvidenceManager::splitData(int s)
{
	int testSetSize=(evidenceSet.size()-validationIndex.size())/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	for(int i=0;i<evidenceSet.size();i++)
	{
		if(validationIndex.find(i)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
			testIndex[i]=0;
		}
		else
		{
			trainIndex[i]=0;
		}
		m++;
	}
	return 0;
}*/


int 
EvidenceManager::splitData(int s)
{
	int testSetSize=(evidenceSet.size()-validationIndex.size())/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	if(foldCnt==1)
	{
		//testStartIndex=-1;
		//testEndIndex=-1;
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		if(randseed<0)
		{
			randseed=getpid();
		}
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size());	
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if(validationIndex.find(eInd)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
			testIndex[eInd]=0;
		}
		else
		{
			trainIndex[eInd]=0;
		}
		m++;
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}


INTINTMAP& 
EvidenceManager::getTrainingSet()
{
	return trainIndex;
}

INTINTMAP& 
EvidenceManager::getTestSet()
{
	return testIndex;
}

INTINTMAP&
EvidenceManager::getValidationSet()
{	
	return validationIndex;
}

/*
int 
EvidenceManager::populateEvidence_Continuous(Evidence** evid,const char* evidStr)
{
	//first check for validity of evidStr
	if(strchr(evidStr,'=')==NULL)
	{
		return -1;
	}
	*evid=new Evidence;
	
	int currInd=0;
	int ttInd=0;
	int tokId=0;
	char tempTok[256];
	while(evidStr[currInd]!='\0')
	{
		if((evidStr[currInd]=='=') || 
		   (evidStr[currInd]==']') ||
		   (evidStr[currInd]==',')
		  )
		{
			tempTok[ttInd]='\0';
			ttInd=0;
			if(tokId==0)
			{
				//This is the variable
				int vId=atoi(tempTok);
				Variable* var=vMgr->getVariableAt(vId);
				(*evid)->assocVariable(vId);
			}
			else
			{
				char* pos=strchr(tempTok,'|');
				//Hard evidence
				if(pos==NULL)
				{
					double varVal=log(atof(tempTok));
					(*evid)->setEvidVal(varVal);
					int vId=(*evid)->getAssocVariable();
					if(maxAbsValues.find(vId)==maxAbsValues.end())
					{
						maxAbsValues[vId]=fabs(varVal);
					}
					else
					{
						double mVal=maxAbsValues[vId];
						if(mVal<fabs(varVal))
						{
							maxAbsValues[vId]=fabs(varVal);
						}
					}
				}
			}
			tokId++;
		}
		else if(evidStr[currInd]!='[')
		{
			tempTok[ttInd]=evidStr[currInd];
			ttInd++;
		}
		currInd++;
	}
	return 0;
}
*/

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, INTINTMAP& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(randInds.find(rind)!=randInds.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		randInds[rind]=0;
	}
	return 0;
}


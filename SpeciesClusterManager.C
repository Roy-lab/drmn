
/*
CMINT: An algorithm to cluster functional omics data from multiple celltypes
    Copyright (C) 2016 Sushmita Roy sushroy@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <ctime>
#include <algorithm>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "Matrix.H"
#include "GeneExpManager.H"
#include "SpeciesDistManager.H"
#include "GeneTree.H"
#include "GeneTreeManager.H"
#include "Gamma.H"
#include "GammaManager.H"
#include "GeneNameMapper.H"
#include "SpeciesFeatureManager.H"
#include "Error.H"
//#include "Evidence.H"
#include "EvidenceManager.H"
#include "DRMNPotential.H"
#include "Expert.H"
#include "SpeciesDistManager.H"
#include "SpeciesClusterManager.H"

SpeciesClusterManager::SpeciesClusterManager()
{
	fixCov=false;
	randFlag=false;
	r=NULL;
	gnm.readGeneNames();
	rseed=-1;

	regsPerIter=5;
	maxEMINTIter=10;
	maxDRMNIter=10; // original experiments set this to 10

	initExpertsPerSpecies=false; // default: init to source species
	learnMode = GREEDY;
	rho1 = 0;
	rho2 = 0;
	rho3 = 0;
}

SpeciesClusterManager::~SpeciesClusterManager()
{
}


int 
SpeciesClusterManager::setOrthogroupReader(MappedOrthogroupReader* aPtr)
{
	mor=aPtr;

	// Set OGID for module ID -- first unused integer starting with 0.
	// Set up 
	setUpModuleOGIDs();
	//cout << "SET UP MappedOrthogroups!!" << endl;
	return 0;
}


int
SpeciesClusterManager::setGammaManager(GammaManager* aPtr)
{
	gammaMgr=aPtr;
	return 0;
}

int 
SpeciesClusterManager::setSrcSpecies(const char* aName)
{
	strcpy(srcSpecies,aName);
	return 0;
}

int 
SpeciesClusterManager::setTrainOGIDs(map<int,int>& ogidSet)
{
	trainSetOGIDs.clear();
	for(map<int,int>::iterator oIter=ogidSet.begin();oIter!=ogidSet.end();oIter++)
	{
		trainSetOGIDs[oIter->first]=oIter->second;	
	}
	return 0;
}

int 
SpeciesClusterManager::getTestPerformance(const char* outputDir, map<int,int>& ogidSet)
{
	
	// Only test on the testsetOGIDS that actually have input cluster assignments
	// We figured this out in readClusters.

	// repopulate with new test set
	// do we ever use testSetOGIDs??
	testSetOGIDs.clear();
	for(map<int,int>::iterator oIter=ogidSet.begin();oIter!=ogidSet.end();oIter++)
	{
		
		testSetOGIDs[oIter->first]=oIter->second;
	}
	
	
	//We want to predict the test expression and we also want to get some kind of a score.
	// Forget about the score for now
	//First get the probabilities for each individual nodes. And then get the joint assignment across 
	//all nodes.
	map<int,map<string,int>*> testClusterAssignments;
	map<int,MappedOrthogroup*>& allOGs=mor->getMappedOrthogroups();
	for(map<int,int>::iterator oIter=ogidSet.begin();oIter!=ogidSet.end();oIter++)
	{
		if(allOGs.find(oIter->first)==allOGs.end())
		{
			cout <<"No og by ID " << oIter->first <<endl;
			exit(0);
		}
		
		int ogid=oIter->first;
		MappedOrthogroup* mgrp=allOGs[ogid];
		map<int,map<string,string>*>& aset=mgrp->getGeneSets();
		for(map<int,map<string,string>*>::iterator gIter=aset.begin();gIter!=aset.end();gIter++)
		{
			
			// need to make sure we have any genes in this ogid/species
			bool doOGID=false; // set to true once we find a match
		
		
			map<string,string>* speciesGeneMap=gIter->second;
			for(map<string,string>::iterator sIter=speciesGeneMap->begin();sIter!=speciesGeneMap->end();sIter++)
			{
				string& specName=(string&)sIter->first;
				string& genename=(string&)sIter->second;
				
				// don't try to do inference or anything for
				// genes that are in the ogids file but for which 
				// we have no species data.
				// this doesn't fix the problem...
				// don't try to do anything with genes that are in the OGIDS file but
				// (a) have no expression data
				// (b) have no initial cluster assignment
				if (legalTestOGIDs.find(specName) == legalTestOGIDs.end() )
				{
					//cout << "Illegal test OGID " << ogid << " for species " << specName << endl;
					continue;
				}
				map<int,int>* specLegal=legalTestOGIDs[specName];
				
				//if (legalTestOGIDs[specName]->find(ogid)=i=(*legalTestOGIDs[specName])->end())
				if (specLegal->find(ogid) == specLegal->end())
				{
					continue;
				}
				
				//init gamma needs a cluster assignment, but this can be arbitrary for a test gene
				gammaMgr->initGamma(ogid,genename,specName,0);
				doInferenceForGene(genename,specName);
				doOGID=true;
			}	
			//Now call gamma manager to do inference jointly across all cell types 
			
			// skip if ogid not done
			if (!doOGID)
			{
				continue;
			}
			Gamma* gamma=gammaMgr->getGammaForOGID(ogid);
			gammaMgr->estimateNonLeafPosterior(gamma->root);
			map<string,int>* maxAssignment=new map<string,int>;
			gammaMgr->getMaximalAssignment(gamma->root,maxAssignment,1);
			testClusterAssignments[ogid]=maxAssignment;
		}
	}
	map<string,ofstream*> filePtrs;
	map<string,ofstream*> mapFilePtrs; // pointers to files for just the predicted expression for the best module
	//Once the test genes are inferred we can write out the test predictions per species
	for(map<string,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		char aFName[1024];
		sprintf(aFName,"%s/%s_pred_test.txt",outputDir,eIter->first.c_str());
		ofstream* oFile=new ofstream(aFName);
		filePtrs[eIter->first]=oFile;

		// predicted, species-specific expression
		sprintf(aFName,"%s/%s_pred_test_maxclust_exprtab.txt",outputDir,eIter->first.c_str());
		ofstream* mapFile=new ofstream(aFName);
		mapFilePtrs[eIter->first]=mapFile;
		(*mapFile) << "Gene\t" << eIter->first << endl;
	}
	for(map<int,map<string,int>*>::iterator oIter=testClusterAssignments.begin();oIter!=testClusterAssignments.end();oIter++)
	{
		map<string,int>* geneSet=oIter->second;
		for(map<string,int>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			//Get the gene name and the cell type
			char tempBuff[1024];
			strcpy(tempBuff,gIter->first.c_str());
			char* tok=strchr(tempBuff,':');
			if(tok==NULL)
			{
				cout <<"Bad format in gene name " << tempBuff<<endl;
				exit(0);
			}
			*tok='\0';
			string geneName(tempBuff);
			string specName(tok+1);
			ofstream* oFileP=filePtrs[specName];
			ofstream* mapFileP=mapFilePtrs[specName];
			int clusterID=gIter->second;
			EvidenceManager* evMgr=evMgrSet[specName];
			EMAP* evidSet=evMgr->getEvidenceAt(geneName);
			int varID=varNameIDMap["Expression"];
			//Evidence* evid=(*evidSet)[varID];
			//double eval=evid->getEvidVal();
			double eval=(*evidSet)[varID];
			(*oFileP) <<geneName <<"\t"<< eval <<"\t" <<clusterID;
			(*mapFileP) << geneName;
			map<int,Expert*>* expertSet=speciesExpertSet[specName];
			for(map<int,Expert*>::iterator fIter=expertSet->begin();fIter!=expertSet->end();fIter++)
			{
				Expert* f=fIter->second;
				DRMNPotential* dPot=f->getDRMNPotential();
				double pred=dPot->predictSample(evidSet);
				(*oFileP) <<"\t"<<pred;

				if (fIter->first==clusterID)
				{
					(*mapFileP) << "\t" << pred << endl;
				}
			}
			(*oFileP) <<endl;
		}
	}
	for(map<string,ofstream*>::iterator fIter=filePtrs.begin();fIter!=filePtrs.end();fIter++)
	{
		fIter->second->close();
		delete fIter->second;
	}
	filePtrs.clear();
	for(map<string,ofstream*>::iterator fIter=mapFilePtrs.begin();fIter!=mapFilePtrs.end();fIter++)
	{
		fIter->second->close();
		delete fIter->second;
	}
	mapFilePtrs.clear();

		
	return 0;
}


int
SpeciesClusterManager::clear()
{
	//I think we need to essentially clear the experts and speciesClusterSet_Genewise and OGwise? 
	// We also need to clear the GammaManager.

	for(map<string,CLUSTERSET*>::iterator eIter=speciesExpertSet.begin();eIter!=speciesExpertSet.end();eIter++)
	{
		CLUSTERSET* speciesExpert=eIter->second;
		for(CLUSTERSET_ITER cIter=speciesExpert->begin();cIter!=speciesExpert->end();cIter++)
		{
			Expert* e=cIter->second;
			delete e;
		}
		speciesExpert->clear();
		delete speciesExpert;
	}
	for(map<string,map<string,int>*>::iterator eIter=speciesClusterSet_Genewise.begin();eIter!=speciesClusterSet_Genewise.end();eIter++)
	{
		map<string,int>* set=eIter->second;
		set->clear();
		delete set;
	}
	speciesClusterSet_Genewise.clear();

	trainSetOGIDs.clear();
	testSetOGIDs.clear();
	
	// clear the legal test IDs
	for(map<string, map<int,int>*>::iterator lIter=legalTestOGIDs.begin(); lIter!=legalTestOGIDs.end();lIter++)
	{
		map<int,int>* set=lIter->second;
		set->clear();
		delete set;
	}
	legalTestOGIDs.clear();
	
	
	// clear out the gammas	
	gammaMgr->clearGammaSet();
	
	
	return 0;
}

int
SpeciesClusterManager::setMaxClusterCnt(int k)
{
	maxClusterCnt=k;
	return 0;
}

int 
SpeciesClusterManager::setRandom(bool flag)
{
	randFlag=flag;
	return 0;
}

int
SpeciesClusterManager::setRandSeed(int aseed)
{
	rseed=aseed;
	return 0;
}

/*
 * Sets constant covariance value and flag that we will use it.
 */
int
SpeciesClusterManager::setConstCov(double val)
{
        constCov=val;
	fixCov=true;
        return 0;
}

/*
 * Sets whether to initialize cluster params within each cell type (true)
 * or to source species params (false)
 */
int
SpeciesClusterManager::setInitExpertsPerSpecies(bool setPerSpecies)
{
        initExpertsPerSpecies=setPerSpecies;
        return 0;
}


/**
* Read species data -- now updated to read species-specific network data.
* Config file tokens:
* Species name, initial clusters, expression data, network data
* SR: We will also keep a variable ID Name mapping here
* **This is where the random seed actually gets set if the user provides one.**
* If seed is -1, then uses PID.
*/
int 
SpeciesClusterManager::readSpeciesData(const char* clusterFName)
{
	vector<string> empty;
	return readSpeciesData(clusterFName, empty);
}

/**
* Version that allows us to only read data for a subset of species
*/
int 
SpeciesClusterManager::readSpeciesData(const char* clusterFName, vector<string>& specNames)
{
	clock_t t=clock(); // time this step

	ifstream inFile(clusterFName);
	char buffer[1024];
	r=gsl_rng_alloc(gsl_rng_default); // RNG initialized here
	if(rseed==-1)
	{
		rseed=getpid();
	}
	gsl_rng_set(r,rseed); // OK, random seed is finally set here.
	cout << rseed << endl;

	bool failed=false; // some data problem

	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string speciesName;
		string clusterFName;
		string expressionFName;
		string netFeatFName=""; // cisregulatory features
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				speciesName.append(tok);
			}
			else if(tokCnt==1)
			{
				clusterFName.append(tok);
			}
			else if(tokCnt==2)
			{
				expressionFName.append(tok);
			}
			else if (tokCnt==3)
			{
				netFeatFName.append(tok);	
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		// skip if we specified a list of species and that list does not contain this species
		if (specNames.size() !=0 && (find(specNames.begin(), specNames.end(), speciesName) == specNames.end()) )
		{
			cout << "Skipping config file line for species " << speciesName << endl; 
			continue; 
		}

		//populate ONLY if this does not exist
		if(speciesExprSet.find(speciesName)==speciesExprSet.end())
		{
			GeneExpManager* exprManager=new GeneExpManager;
			//exprManager->readExpression(expressionFName.c_str());
			exprManager->readExpression_Withheader(expressionFName.c_str());
			speciesExprSet[speciesName]=exprManager;
			int vID=0;
			//harcoding this now
			varNameIDMap["Expression"]=vID;
			varIDNameMap[vID]="Expression";
			vID++;
			// read cisreg features
			SpeciesFeatureManager* featureManager=new SpeciesFeatureManager;
			cout << netFeatFName << endl;
			if (netFeatFName.compare("") != 0)
			{
				//featureManager->readFeatures(netFeatFName.c_str());
				featureManager->readFeatures_Efficient(netFeatFName.c_str());
			}
			else
			{
				cout << "Missing cis-regulatory feature data for " << speciesName << endl;
				failed=true;
			}
			
			//speciesFeatSet[speciesName]=featureManager;
			//SR: We will maintain the names of all regulators in one place
			map<string,int>& regSet=featureManager->getFeatureNames();
			for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
			{
				allRegulatorSet[rIter->first]=0;
				if(varNameIDMap.find(rIter->first)==varNameIDMap.end())
				{
					varNameIDMap[rIter->first]=vID;
					varIDNameMap[vID]=rIter->first;
					vID++;
				}
			}	
			//Let's combine expression and features into evidenceManager like object.
			EvidenceManager* evManager=new EvidenceManager;
			evMgrSet[speciesName]=evManager;
			//Now for every gene, get it's expression and its species features. We will create an evidenceSet for each gene and an evidence 
			//expression or the regulatory features.
			map<string,vector<double>*>& geneSet=exprManager->getGeneSet();
			for(map<string,vector<double>*>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				vector<double>* exp=gIter->second;
				EMAP* evidMap=new EMAP;
				//If there are multiple measures for this gene, we are just going to use the mean
				double collapsedVal=0;
				for(int i=0;i<exp->size();i++)
				{
					collapsedVal=collapsedVal+(*exp)[i];
				}
				collapsedVal=collapsedVal/exp->size();
				//Now make an Evidence object
				//Evidence* evid=new Evidence;
				//int varID=varNameIDMap["Expression"];
				////evid->assocVariable(varID);
				//evid->setEvidVal(collapsedVal);
				//(*evidMap)[varID]=evid;
				int varID=varNameIDMap["Expression"];
				(*evidMap)[varID]=collapsedVal;
				//Now for this gene get all its regulatory features	
				//map<string,double>* tfs=featureManager->getRegulatorsForTarget((string&)gIter->first);
				double* tfs=featureManager->getRegulatorsForTarget_Efficient((string&)gIter->first);
				//for(map<string,map<string,double>*>::iterator rIter=predictiveFeatureSet.begin();rIter!=predictiveFeatureSet.end();rIter++)
				if(tfs!=NULL)
				{
					//for(map<string,double>::iterator rIter=tfs->begin();rIter!=tfs->end();rIter++)
					for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
					{
						//map<string,double>* tgts=rIter->second;
						//double regFeatVal=rIter->second;
						double regFeatVal=tfs[rIter->second];
						/*if(tgts->find(gIter->first)==tgts->end())
						{	
							regFeatVal=featureManager->getDefaultValue();
						}
						else
						{
							regFeatVal=(*tgts)[gIter->first];
						}*/
						//Evidence* evid=new Evidence;
						//int varID=varNameIDMap[rIter->first];
						////evid->assocVariable(varID);
						//evid->setEvidVal(regFeatVal);
						//(*evidMap)[varID]=evid;
						int varID=varNameIDMap[rIter->first];
						(*evidMap)[varID]=regFeatVal;
					}
					//Now we just add the evidMap to evMgr;
					//evManager->addEvidenceWithName(evidMap,(string&)gIter->first);
				}
				//else
				//{
				//	//We didn't have features, we skip this one
				//	delete evidMap;
				//}
				evManager->addEvidenceWithName(evidMap,(string&)gIter->first);
			}
			delete featureManager;
		}
		//Only this needs to be read each time
		int readOK = readClusters(speciesName,clusterFName.c_str());
		if (readOK != 0) 
		{
			cerr << "ClusterID exceeds max allowed." << endl;
			failed=true;
			break;
		}
	}
	inFile.close();

	int ri=0;
	for(map<string,int>::iterator rIter=allRegulatorSet.begin();rIter!=allRegulatorSet.end();rIter++)
	{
		allRegulatorIndex[rIter->first] = ri;
		ri++;
	}

	if (failed)
	{
		return 1;
	}
	//Do some other inits to learn the params
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);

	//cout <<"Read expression, clusterings, cis-regulatory feature data AND CREATING EVMGR for " << speciesClusterSet_Genewise.size() << " species" << endl;
	t=clock()-t;
	cout << "Time to read all input data for " << speciesClusterSet_Genewise.size() << " species (m) = " << ((long double)t/CLOCKS_PER_SEC)/60.0 << endl;
	return 0;
}

int
SpeciesClusterManager::initExperts()
{
	for(map<string,CLUSTERSET*>::iterator aIter=speciesExpertSet.begin();aIter!=speciesExpertSet.end();aIter++)
	{
		CLUSTERSET* speciesExpert=aIter->second;
		for(CLUSTERSET_ITER cIter=speciesExpert->begin();cIter!=speciesExpert->end();cIter++)
		{
			Expert* e=cIter->second;
			// init cluster params per species, by its own data?
			if (initExpertsPerSpecies)
			{
				estimateMeanCov(e,(string&)aIter->first,cIter->first);
			}
			else // Update to initialize expert params to source species?
			{
				std::string srcSpecStr(srcSpecies);
				estimateMeanCov(e,srcSpecStr,cIter->first);
			}
		}
	}
	return 0;
}

/*
int 
SpeciesClusterManager::getExtantClusterAssignment(map<int,map<string,int>*>& clusterAssignments)
{
	for(map<int,map<string,INTDBLMAP*>*>::iterator gIter=gammas.begin();gIter!=gammas.end();gIter++)
	{
		map<string,INTDBLMAP*>* gValSet=gIter->second;
		map<string,int>* clustermembers=NULL;
		if(clusterAssignments.find(gIter->first)==clusterAssignments.end())
		{
			clustermembers=new map<string,int>;
			clusterAssignments[gIter->first]=clustermembers;
		}
		else
		{
			clustermembers=clusterAssignments[gIter->first];
		}
		for(map<string,INTDBLMAP*>::iterator sIter=gValSet->begin();sIter!=gValSet->end();sIter++)
		{
			//For each species get the best cluster id for this og
			INTDBLMAP* gamma_i=sIter->second;
			double maxPval=-1;
			int maxexpertID=-1;
			for(INTDBLMAP_ITER eIter=gamma_i->begin();eIter!=gamma_i->end();eIter++)
			{
				if(eIter->second>maxPval)
				{
					maxPval=eIter->second;
					maxexpertID=eIter->first;
				}
			}
			(*clustermembers)[sIter->first]=maxexpertID;
		}
	}
	return 0;
}*/


map<string,CLUSTERSET*>& 
SpeciesClusterManager::getExtantSpeciesClusters()
{
	return speciesExpertSet;
}

/**
* This is EMINT. We use this to initialize everything.
*
*/
int
SpeciesClusterManager::estimateExpertParameters(const char* outputDir)
{
	double currScore=0;
	bool convergence=false;
	//bool convergence=true;
	int iter=0;
	//while((!convergence)&&(iter<100))
	//while((!convergence)&&(iter<10)) // 
	//while((!convergence)&&(iter<2))
	while((!convergence)&&(iter<maxEMINTIter)) 
	{
		cout << "EMINT expectation step: compute gammas and transition matrices." << endl;
		expectationStep();

		cout << "Maximization: Update module parameters for EMINT." << endl;
		maximizationStep();


		//dumpAllInferredClusterAssignments(outputDir,iter);
		double newScore=getScore();
		double diff=fabs(newScore-currScore);
		if((iter>0) && (diff<0.5))
		{
			convergence=true;
		}
		cout <<"Iter: " << iter << " score: " << newScore << " diff " << diff << endl;
		currScore=newScore;
		iter++;
	}
	if(convergence)
	{
		cout <<"EMINT Convergence at iteration:" << iter << endl;
	}
	return 0;
}

/**
* Run DRMN.
* Draft 1: Keep EMINT machinery mostly the same -- run expectationStep and maximizationStep. However, each iteration we now 
* make a hard assignment to modules, which we then pass to MRTLE.
* Yes, this is used!
*/
int
SpeciesClusterManager::estimateDRMN(const char* outputDir)
{
	int success=0;

	double currScore=0;
	bool convergence=false;
	int iter=0;
	//This will store the empty graph per cell type where for each module there will be no regulators.
	precomputeEmptyGraphPrior();
	//while((!convergence)&&(iter<100))
	//while((!convergence)&&(iter<10))  // SET THIS for number of iterations // MAX ITER
	while((!convergence)&&(iter<maxDRMNIter) && success==0)
	{
		// At this stage we have a set of modules, so we are going to attach regulators to each module
		// This is like doing the M step in DRMN
		estimateRegProgs(iter);
		// Hard module assignment and recompute transition probs.
		expectationStep_DRMN();
		//re-estimate transition probabilities. Not needed because we are estimating transition probs in the expectationStep_DRMN
		//reestimateTransitionProbabilities();
		// show output: Don't really need this unless we are debugging
		dumpAllInferredClusterAssignments(outputDir,iter);
		double newScore=getDRMNScore();
		double diff=fabs(newScore-currScore);
		if((iter>0) && (diff<0.5))
		{
			convergence=true;
		}
		cout <<"Iter: " << iter << " score: " << newScore << " diff " << diff << endl;
		currScore=newScore;
		//Do the hard assignment
		assignGenesToExperts_FromMap();

		// Make output directory -- and kill the program immediately if we can't!
		char intermediateOutDir[1024];
		sprintf(intermediateOutDir,"%s/iter%d",outputDir,iter);
		char command[1024];

		sprintf(command,"mkdir -p %s",intermediateOutDir);
		const int mkdir_err = system(command);
		if (mkdir_err != 0)
		{
			cerr << "Cannot create output directory " << intermediateOutDir << " !!" << endl;
			//success=1; 
		}
		success=abs(success)+abs(showDRMNResults(intermediateOutDir));
		iter++;
	}
	if(convergence)
	{
		cout <<"DRMN Convergence at iteration:" << iter << endl;
	} 
	else 
	{
		cout <<"DRMN hit max iters " << iter << endl;
	}

	// had a problem earlier
	if (success != 0) return success;


	success=showDRMNResults(outputDir);
	return success;
}


int
SpeciesClusterManager::showDRMNResults(const char* drmnOutputDir)
{
	int errors=0;
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		char outputPath[1024];
		sprintf(outputPath,"%s/%s_net.txt",drmnOutputDir,cIter->first.c_str());
		CLUSTERSET* expertSet=cIter->second;
		int errcode = showDRMNModelForCellType(outputPath,expertSet);
		errors=errors+abs(errcode);

		sprintf(outputPath,"%s/%s_pred.txt",drmnOutputDir,cIter->first.c_str());
		errcode=showDRMNPredictionsForCellType(outputPath,expertSet,(string&)cIter->first);
		errors=errors+abs(errcode);
		
		sprintf(outputPath,"%s/%s_pred_maxclust_exprtab.txt",drmnOutputDir,cIter->first.c_str());
		errcode=showBestDRMNPredictionsForCellType(outputPath,expertSet,(string&)cIter->first);
		errors=errors+abs(errcode);

		sprintf(outputPath,"%s/%s_sample_maxclust_exprtab.txt",drmnOutputDir,cIter->first.c_str());
		errcode=showBestDRMNSamplesForCellType(outputPath, expertSet, (string&)cIter->first);
		errors=errors+abs(errcode);
		
	}
	//SR added show the transition matrices
	SpeciesDistManager* sdMgr=gammaMgr->getSpeciesDistManager();
	sdMgr->showInferredConditionals(drmnOutputDir);
	if (errors!=0) return 1;
	else return 0;
}

int
SpeciesClusterManager::showDRMNModelForCellType(const char* pathName,map<int,Expert*>* expertSet)
{
	ofstream oFile(pathName);
	for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
	{
		Expert* e=eIter->second;
		DRMNPotential* dPot=e->getDRMNPotential();
		INTDBLMAP& potWts=dPot->getCondWeight();
		for(INTDBLMAP_ITER wIter=potWts.begin();wIter!=potWts.end();wIter++)
		{
			oFile <<"Module"<<eIter->first<<"\t"<< varIDNameMap[wIter->first] << "\t" << wIter->second << endl;
		}
	}
	oFile.close();
	if (oFile.bad())
	{
		cerr << "Error writing to " << pathName << endl;
		return 1; 
	}
	return 0;
}

int
SpeciesClusterManager::showDRMNPredictionsForCellType(const char* pathName,map<int,Expert*>* expertSet,string& cellType)
{
	ofstream oFile(pathName);
	//Let's print the gene's actual expression and the predicted expression from the experts
	oFile<<"Gene\tTrueExpr\tModuleID";
	for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
	{
		oFile<<"\tPred"<< eIter->first;
	}
	oFile<<endl;

	
	EvidenceManager* evMgr=evMgrSet[cellType];
	for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
	{
		Expert* e=eIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			EMAP* evidSet=evMgr->getEvidenceAt((string&)gIter->first);
			int varID=varNameIDMap["Expression"];
			//Evidence* evid=(*evidSet)[varID];
			//double eval=evid->getEvidVal();
			double eval=(*evidSet)[varID];
			oFile <<gIter->first<<"\t"<< eval <<"\t" <<eIter->first;
			for(map<int,Expert*>::iterator fIter=expertSet->begin();fIter!=expertSet->end();fIter++)
			{
				Expert* f=fIter->second;
				DRMNPotential* dPot=f->getDRMNPotential();
				double pred=dPot->predictSample(evidSet);
				oFile<<"\t"<<pred;
			}
			oFile <<endl;
		}
	}
	oFile.close();
	if (oFile.bad())
	{
		cerr << "Error writing to " << pathName << endl;
		return 1; 
	}

	return 0;
}

/**
 * Prints the predictions for the best module for each training gene.
 */
int
SpeciesClusterManager::showBestDRMNPredictionsForCellType(const char* pathName,map<int,Expert*>* expertSet,string& cellType)
{
	ofstream oFile(pathName);
	oFile<<"Gene\t" << cellType << endl;

	EvidenceManager* evMgr=evMgrSet[cellType];
	for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
	{
		Expert* e=eIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			EMAP* evidSet=evMgr->getEvidenceAt((string&)gIter->first);
			int varID=varNameIDMap["Expression"];
			//Evidence* evid=(*evidSet)[varID];
			//double eval=evid->getEvidVal();
			double eval=(*evidSet)[varID];
			DRMNPotential* dPot=e->getDRMNPotential();
			double pred=dPot->predictSample(evidSet);
			oFile <<gIter->first<<"\t"<< pred << endl;
		}
	}
	oFile.close();
	if (oFile.bad())
	{
		cerr << "Error writing to " << pathName << endl;
		return 1; 
	}

	return 0;
}

/**
 * Prints a sample for the best module for each training gene.
 */
int
SpeciesClusterManager::showBestDRMNSamplesForCellType(const char* pathName,map<int,Expert*>* expertSet,string& cellType)
{
	ofstream oFile(pathName);
	oFile<<"Gene\t" << cellType << endl;

	EvidenceManager* evMgr=evMgrSet[cellType];
	for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
	{
		Expert* e=eIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			EMAP* evidSet=evMgr->getEvidenceAt((string&)gIter->first);
			int varID=varNameIDMap["Expression"];
			//Evidence* evid=(*evidSet)[varID];
			//double eval=evid->getEvidVal();
			double eval=(*evidSet)[varID];
			DRMNPotential* dPot=e->getDRMNPotential();
			// sample from gaussian using rng r (owned by speciesclustermanager)
			double pred=dPot->generateSample(evidSet, r);
			oFile <<gIter->first<<"\t"<< pred << endl;
		}
	}
	oFile.close();
	if (oFile.bad())
	{
		cerr << "Error writing to " << pathName << endl;
		return 1; 
	}

	return 0;
}


/**
* Makes hard assignments to modules and updates transition probabilities based on the max assignments.
* DC ADD for DRMN
*/
int
SpeciesClusterManager::setMaxAssignments()
{
	// mappedClusterAssignment is a member variable -- we're updating it 
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,true);
	gammaMgr->showClusterFlipCnts();
	gammaMgr->reestimateTransitionProbability();
	
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,false);
	gammaMgr->reestimateTransitionProbability();

	// Update assignment of genes to experts
	assignGenesToExperts_FromMap();
	//Initialze the score for each expert now
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		EvidenceManager* evMgr=evMgrSet[cIter->first];
		CLUSTERSET* expertSet=cIter->second;
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			Expert* e=eIter->second;
			//get the genes assigned to this expert and get the score for these genes
			map<string,int>& geneSet=e->getGeneSet();
			//need to worry about the fact that the old expert's mean was multi-variate possibly?
			double mean=0;
			double cv=0;
			int varID=varNameIDMap["Expression"];
			for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				EMAP* evidMap=evMgr->getEvidenceAt((string&)gIter->first);
				//Only care about expression
				//Evidence* evid=(*evidMap)[varID];
				//double v=evid->getEvidVal();
				double v=(*evidMap)[varID];
				mean=mean+v;
			}	
			mean=mean/geneSet.size();
			for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				EMAP* evidMap=evMgr->getEvidenceAt((string&)gIter->first);
				//Only care about expression
				//Evidence* evid=(*evidMap)[varID];
				//double v=evid->getEvidVal();
				double v=(*evidMap)[varID];
				double diff=mean-v;
				cv=cv+(diff*diff);
			}	
			cv=cv/(geneSet.size()-1);
			DRMNPotential* pot=new DRMNPotential;	
			pot->setAssocVariable(varID,DRMNPotential::FACTOR);
			pot->potZeroInit();
			pot->updateMean(varID,mean);
			pot->updateCovariance(varID,varID,cv);
			pot->initMBCovMean();
			e->setDRMNPotential(pot);
			double pll=0;
			for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				EMAP* evidMap=evMgr->getEvidenceAt((string&)gIter->first);
				//Only care about expression
				//Evidence* evid=(*evidMap)[varID];	
				double ll=pot->getCondPotValueFor(evidMap);
				pll=pll+log(ll);
			}	
			double varCnt=pot->getAssocVariables().size();
			double paramCnt=2*varCnt;
			paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
			double p_score=pll-((paramCnt/2)*log(geneSet.size()));
			e->setLLScore(pll);
			//e->setLLScore(p_score);
		}
	}
	return 0;
}

//DRMN function needed to add the multi-cell type graph prior
//This is needed only for the regulator selection to make sure the regulators are varying gradually
//Since this is the same species, we can assume that the regulators/elements are named the same
int
SpeciesClusterManager::precomputeEmptyGraphPrior()
{
	SpeciesDistManager* sdMgr=gammaMgr->getSpeciesDistManager();
	
	STRINTMAP edgeStatus;
	for(map<string,GeneExpManager*>::iterator sIter=speciesExprSet.begin();sIter!=speciesExprSet.end();sIter++)
	{
		edgeStatus[sIter->first]=0;
	}
	//For each module and regulator we need to add an edge per cell type
	for(int k=0;k<maxClusterCnt;k++)
	{
		for(map<string,int>::iterator rIter=allRegulatorSet.begin();rIter!=allRegulatorSet.end();rIter++)
		{
			char regModulePair[1024];
			sprintf(regModulePair,"%s-Module%d",rIter->first.c_str(),k);
			string key(regModulePair); 
			//Let's compute the prior so that we are sensitive to the module we are using
			double prior=sdMgr->getEdgeStatusProb(edgeStatus,k);
			double logPrior=log(prior);
			regModulePairSetPrior[key]=logPrior;
			if(edgeCelltypeMap.find(key)==edgeCelltypeMap.end())
			{
				STRINTMAP* edgeSpecAssign=new STRINTMAP;
				for(map<string,GeneExpManager*>::iterator sIter=speciesExprSet.begin();sIter!=speciesExprSet.end();sIter++)
				{
					(*edgeSpecAssign)[sIter->first]=0;
				}
				edgeCelltypeMap[key]=edgeSpecAssign;
			}
		}
	}
	//Done with the creation of the empty graph. Now compute edge specific prior values		
	edgeStatus.clear();
	return 0;
}

int
SpeciesClusterManager::readClusters(string& specName, const char* fName)
{
	int success=0; // switch to 1 if any problem
	
	EvidenceManager* evMgr=evMgrSet[specName];

	ifstream inFile(fName);
	CLUSTERSET* cset=new CLUSTERSET;
	map<string,int>* geneset=new map<string,int>;
	speciesExpertSet[specName]=cset;
	speciesClusterSet_Genewise[specName]=geneset;

	int sansOGIDs=0; // keep track of how many genes don't have OGIDs
	int sansExpr=0; // keep track of how many genes don't have expression
	int sansFeat=0; // keep track of how many genes don't have features
	
	// we will store legal test OGIDS for this species
	// legal means they have a cluster assignment, are in the test set,
	// and aren't in the training set.
	map<int,int>* legal=new map<int,int>;
	legalTestOGIDs[specName]=legal;
	
	// keep track of the cluster IDs that we've found.
	map<int,int> foundClusterIDs;

	char buffer[1024];
	GeneExpManager* geMgr=speciesExprSet[specName];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int clustid=0;
		string genename;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				genename.append(tok);
			}	
			else if(tokCnt==1)
			{
				clustid=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}

		// is cluster ID legal?
		if (clustid >= maxClusterCnt || clustid < 0)
		{
			cerr << "ClusterID " << clustid << " for " << genename << ":" << specName << " is outside bounds of [0," << maxClusterCnt-1 << "]." << endl;
			success=1;
			break;
		}

		int ogid=mor->getMappedOrthogroupID(genename.c_str(),specName.c_str());
		if(ogid==-1)
		{
			//cout <<"OGID -1 for " << genename << ":" << specName << endl;
			sansOGIDs++;
			continue;
		}
			
		
		if(geMgr->getExp(genename)==NULL)
		{
			sansExpr++;
			continue;
		}
		EMAP* evidSet=evMgr->getEvidenceAt(genename);
		if (evidSet == NULL)
		{
			sansFeat++;
			continue;
		}

		// save legal OGID IF it has expression data AND cluster assignment
		(*legal)[ogid]=1;
		
		if(trainSetOGIDs.size()>0 && trainSetOGIDs.find(ogid)==trainSetOGIDs.end())
		{
			continue;
		}
		
		
		Expert* expert=NULL;
		if(randFlag) // randomize initial clusters
		{
			cout << "Now randomizing input clusters..." << endl;
			if(initRandAssign.find(ogid)==initRandAssign.end())
			{
				//don't use this clustid but some other one
				double pval=gsl_ran_flat(r,0,1);
				double step=1.0/(double)maxClusterCnt;
				int newclustid=(int)(floor(pval/step));
				if(newclustid>=maxClusterCnt)
				{
					newclustid=maxClusterCnt-1;
				}
				clustid=newclustid;
				//initRandAssign[ogid]=clustid;
			}
			else
			{
				clustid=initRandAssign[ogid];
			}
		}
		foundClusterIDs[clustid]++; // tally up genes in this cluster

		if(cset->find(clustid)==cset->end())
		{
			expert=new Expert;
			(*cset)[clustid]=expert;
		}
		else
		{
			expert=(*cset)[clustid];
		}
		(*geneset)[genename]=clustid;
		expert->assignGeneToExpert(genename.c_str());
		gammaMgr->initGamma(ogid,genename,specName,clustid);
		//Gamma* g=gammaMgr->getGammaForOGID(ogid);
		//gammaMgr->estimateNonLeafPosterior(g->root);
	}
	cout <<"Read " << geneset->size() << " genes in " << specName << endl;
	cout << "\tSkipped " << sansOGIDs << " genes that weren't in the OGIDs master list." << endl;
	cout << "\tSkipped " << sansExpr << " genes that didn't have expression data." << endl;
	cout << "\tSkipped " << sansFeat << " genes that didn't have feature data." << endl;
	inFile.close();

	// initialize empty experts for remaining IDs?
	// DC Dec 2017
	/*for (int k=0; k < maxClusterCnt; k++)
	{
		// no genes here?
		if (foundClusterIDs[k]==0)
		{
			cout << "Initializing empty cluster " << k << endl;
			Expert* expert=new Expert;
			(*cset)[k]=expert;
		}
	}*/


	return success;
}

int 
SpeciesClusterManager::maximizationStep()
{
	for(map<string,CLUSTERSET*>::iterator aIter=speciesExpertSet.begin();aIter!=speciesExpertSet.end();aIter++)
	{
		CLUSTERSET* speciesExpert=aIter->second;
		for(CLUSTERSET_ITER cIter=speciesExpert->begin();cIter!=speciesExpert->end();cIter++)
		{
			Expert* e=cIter->second;
			estimateMeanCov(e,(string&)aIter->first,cIter->first);
		}
	}
	return 0;
}


/**
* EMINT original expectation step.
* Computes gammas and transition probabilities.
*/
int 
SpeciesClusterManager::expectationStep()
{
	//First get the gammas for each leaf node
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* speciesClusters=sIter->second;
		expectationStep_Species((string&)sIter->first,speciesClusters);
	}
	//The use the speciesdist manager's conditional distributions to infer the rest of the gammas.
	gammaMgr->estimateNonLeafPosterior();
	gammaMgr->estimateTransitionProbability();
	return 0;
}

/**
* DRMN expectation step.
* Computes gammas 
*/
int
SpeciesClusterManager::expectationStep_DRMN()
{
	cout << "Doing DRMN expectation step." << endl;
	cout << "Getting observed probs from experts." << endl;
	//First get observed probs from each expert.
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* speciesClusters=sIter->second;
		expectationStep_DRMN_Species((string&)sIter->first,speciesClusters);
	}
	
	cout << "Computing gammas and normterms as conditionals" << endl;
	//This would probably be the same, I don't see why it might be different.
	//gammaMgr->estimateNonLeafPosterior_DRMN();
	gammaMgr->estimateNonLeafPosterior();

	//Then use the speciesdist manager's conditional distributions to compute the joint.
	cout << "Doing DRMN - update gammas, normterms to represent joint" << endl;
	gammaMgr->estimateTransitionProbability();

	return 0;
}


/**
* SR: changed the name from maximizationStep2_DRMN to estimateRegProgs. 
* The goal is to go over each module and find the best regulator that explains the expression of the module.
*/
int
SpeciesClusterManager::estimateRegProgs(int iter)
{
	cout << "Doing DRMN maximization step 2: updating regulatory programs." << endl;
	
	// Actually, we'll have one instance of MRTLE per CLUSTER, for all species.
	// First, loop over clusters and print out maintenance probs
	//cout << "Estimating new regulatory programs per cluster." << endl;

	// TODO: problem if no genes in this cluster! skip it.
	for (int k=0; k<maxClusterCnt; k++)
	{
		cout << "Estimating reg programs for all cell types, for state " << k << endl;
		// does this module have an expert initialized? if not, skip it. TODO
        //estimateRegProgs_PerModule(k); // SR: only need to send the cluster ids here. 
		if (learnMode == GREEDY)
		{
			cout << "doing GREEDY!" << endl;
        	estimateRegProgs_PerModule(k); // SR: only need to send the cluster ids here. 
		}
		else
		{
			cout << "doing LASSO!" << endl;
        	estimateRegProgs_PerModule_LASSO(k,iter); // SR: only need to send the cluster ids here. 
		}
	}
    return 0;
}



double 
SpeciesClusterManager::getScore()
{
	double totalComplexity=0;
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
		CLUSTERSET* expertSet=cIter->second;
		map<string,int>* speciesGenes=speciesClusterSet_Genewise[cIter->first];
		int paramCnt=(expertSet->size()*2*(speciesGenes->size()));
		double modelComplexity=(paramCnt/2)*log(speciesGenes->size());
		totalComplexity=totalComplexity+modelComplexity;
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetClip();
		}
		for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
		{
			int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),cIter->first.c_str());
			vector<double>* exprProf=speciesMgr->getExp(vIter->first);
			map<int,double>* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
			//double* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
			for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
			{
				Expert* e=eIter->second;
				double pdf=e->getOutputPDF(exprProf);
				if(std::isnan(pdf))
				{
					cout <<"PDF is nan for " << cIter->first << " " << vIter->first << " for expert " << eIter->first << endl;
				}	
				(*mixOutProbs)[eIter->first]=pdf;
			}
		}
		/*for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			cout <<cIter->first<<":" << eIter->first <<" " << eIter->second->getClip()<< " minpdf " << eIter->second->getMinPDF() << endl;;
		}*/
	}
	double netLL_unpen=gammaMgr->getAllNodeScore();
	double netLL=netLL_unpen+totalComplexity;
	cout <<"Unpenalized score= " << netLL_unpen << " Penalized score=" << netLL << endl;
	return netLL_unpen;
}

//DRMN score computation will be similar to CMINT, except now we will use the cell type specific EvidenceManager
double 
SpeciesClusterManager::getDRMNScore()
{
	double totalComplexity=0;
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
		EvidenceManager* evMgr=evMgrSet[cIter->first];
		CLUSTERSET* expertSet=cIter->second;
		map<string,int>* speciesGenes=speciesClusterSet_Genewise[cIter->first];
		int paramCnt=(expertSet->size()*2*(speciesGenes->size()));
		double modelComplexity=(paramCnt/2)*log(speciesGenes->size());
		totalComplexity=totalComplexity+modelComplexity;
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetClip();
		}
		for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
		{
			int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),cIter->first.c_str());
			EMAP* evidSet=evMgr->getEvidenceAt((string&)vIter->first);
			map<int,double>* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
			//double* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
			for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
			{
				Expert* e=eIter->second;
				double pdf=e->getDRMNProb(evidSet);
				if(std::isnan(pdf))
				{
					cout <<"PDF is nan for " << cIter->first << " " << vIter->first << " for expert " << eIter->first << endl;
				}	
				(*mixOutProbs)[eIter->first]=pdf;
			}
		}
		/*for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			cout <<cIter->first<<":" << eIter->first <<" " << eIter->second->getClip()<< " minpdf " << eIter->second->getMinPDF() << endl;;
		}*/
	}
	double netLL_unpen=gammaMgr->getAllNodeScore();
	double netLL=netLL_unpen+totalComplexity;
	cout <<"Unpenalized score= " << netLL_unpen << " Penalized score=" << netLL << endl;
	return netLL_unpen;
}

//This function simply goes over each cell type for a module and returns the score for the module

double 
SpeciesClusterManager::getDRMNScoreForModule(int moduleID)
{
	double moduleScore=0;
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
		EvidenceManager* evMgr=evMgrSet[cIter->first];
		CLUSTERSET* expertSet=cIter->second;
		Expert* e=(*expertSet)[moduleID];
		// don't bother if expert is null?? No -- this should never happen.
		if (e==NULL)
		{
			cerr << "NULL EXPERT in getDRMNScoreForModule " << moduleID << endl;
		}
		moduleScore=moduleScore+e->getLLScore();
		
	}
	return moduleScore;
}

// Get penalized and unpenalized data likelihood from DRMN model.
// DC added 12/2017
int 
SpeciesClusterManager::getDRMNScore_test(map<int,int>& testOGIDs, double& netLL_unpen, double& netLL)
{
	double totalComplexity=0;

	// iterate over species (cell types)
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		int ogidsForSpecies=0; // We will count up the number of test OGIDs in this species/cell type

		GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
		EvidenceManager* evMgr=evMgrSet[cIter->first];
		CLUSTERSET* expertSet=cIter->second;

		// compute complexity from training data size
		// we instead need to do it for testing data size, yes?
		map<string,int>* speciesGenes=speciesClusterSet_Genewise[cIter->first];
		int paramCnt=(expertSet->size()*2*(speciesGenes->size()));
		double modelComplexity=(paramCnt/2)*log(speciesGenes->size());
		totalComplexity=totalComplexity+modelComplexity;

		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetClip();
		}
		// speciesGenes was populated with training genes only.

		// training data (old)
		//for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
		for (map<int,int>::iterator testIter=testOGIDs.begin();testIter!=testOGIDs.end();testIter++)
		{
			//int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),cIter->first.c_str()); //train only
			//EMAP* evidSet=evMgr->getEvidenceAt((string&)vIter->first);
			//map<int,double>* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
			int ogid=testIter->first;

			// now we need to go through some contortions to get the genename!
			map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
			MappedOrthogroup* mgrp=allOgs[ogid];
			// get 'hits' for this species
			GeneMap* geneMap=mgrp->getSpeciesHits(cIter->first.c_str());
			// I think the keys here finally the gene names!
			map<string,map<string,STRINTMAP*>*>& geneSet=geneMap->getGeneSet();
			// count up for the number of genes here
			ogidsForSpecies+=geneSet.size();
			// now for each of those...
			for (map<string,map<string,STRINTMAP*>*>::iterator nameIter=geneSet.begin(); nameIter!=geneSet.end(); nameIter++)
			{
				string testgene=nameIter->first;
				EMAP* evidSet=evMgr->getEvidenceAt(testgene);
				if (evidSet == NULL)
				{
					continue;
				}
				map<int,double>* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,testgene);
				
				//double* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
				for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
				{
					Expert* e=eIter->second;
					double pdf=e->getDRMNProb(evidSet);
					if(std::isnan(pdf))
					{
						cout <<"PDF is nan for " << cIter->first << ":" << testgene << " for expert " << eIter->first << endl;
					}	
					(*mixOutProbs)[eIter->first]=pdf;
				}
			}
		}


		/*for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			cout <<cIter->first<<":" << eIter->first <<" " << eIter->second->getClip()<< " minpdf " << eIter->second->getMinPDF() << endl;;
		}*/
	}

	// now we take these as args
	//double netLL_unpen=gammaMgr->getAllNodeScore(); // no need to declare
	//double netLL=netLL_unpen+totalComplexity; // no need to declare
	cout << "Submitting " << testOGIDs.size() << " test OGIDs for scoring." << endl;
	netLL_unpen=gammaMgr->getAllNodeScore(testOGIDs);
	netLL=netLL_unpen+totalComplexity;

	//cout <<"Unpenalized DRMN test score= " << netLL_unpen << " Penalized DRMN test score=" << netLL << endl;
	return 0;
}



/*
* EMINT.
* Estimates mean and covariance for each state.
*/
int
SpeciesClusterManager::estimateMeanCov(Expert* e, string& specName, int clusterID)
{
	GeneExpManager* exprMgr=speciesExprSet[specName];
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	vector<double>* expr=exprMgr->getExp(speciesGenes->begin()->first);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	int dim=expr->size();
	Matrix* mean=e->getMean();
	if(mean==NULL)
	{
		mean=new Matrix(1,dim);
	}
	double totaldp=0;
	for(int d=0;d<dim;d++)
	{
		double sum=0;
		double meanval=0;
		for(map<string,int>::iterator dIter=speciesGenes->begin();dIter!=speciesGenes->end();dIter++)
		{
			//Assuming that each evidence is actually a joing assignment to RVs. What we
			//want is a mean and variance for each experiment
			vector<double>* geneExpr=exprMgr->getExp(dIter->first);
			if(geneExpr==NULL)
			{
				cout <<"No gene by name " << dIter->first << " in species " << specName << endl;
			}
			int ogid=mor->getMappedOrthogroupID(dIter->first.c_str(),specName.c_str());
			//This gamma matrix is a kXk matrix which stores the joint probability of being in cluster i
			//given its ancestor was in some other species. 
			Matrix* gamma_i_k_s=gammaMgr->getGamma(ogid,(string&)dIter->first,specName);
			//Need to sum over all possible ways in which we could end up in clusterID
			for(int r=0;r<gamma_i_k_s->getRowCnt();r++)
			{
				double g_i=gamma_i_k_s->getValue(r,clusterID);
				if(g_i>1)
				{
					cout <<"Weird gamma found for " << dIter->first << " at row " << r << endl;
				}
				sum=sum+g_i;
				/*if(sum>10000)
				{
					cout <<"Weird sum " << sum << " found at " << dIter->first << " at row " << r << endl;
				}*/

				double eval=(*geneExpr)[d];
				meanval=meanval+(g_i*eval);
			}
		}
		if(sum==0)
		{
			cout <<"No members in cluster " << clusterID << endl;
		}
		meanval=meanval/sum;
		mean->setValue(meanval,0,d);
		totaldp=sum;
	}
	Matrix* covariance=e->getCovariance();
	if(covariance==NULL)
	{
		covariance=new Matrix(dim,dim);
	}
	covariance->setAllValues(0);
	for(int i=0;i<dim;i++)
	{
		double m1=mean->getValue(0,i);
		//for(int j=i;j<dim;j++)
		for(int j=i;j<i+1;j++)
		{
			double m2=mean->getValue(0,j);
			double cov=0;
			if(i==j)
			{
				cov=0.00001;
			}
			double sum=0;
			for(map<string,int>::iterator dIter=speciesGenes->begin();dIter!=speciesGenes->end();dIter++)
			{
				vector<double>* geneExpr=exprMgr->getExp(dIter->first);
				int ogid=mor->getMappedOrthogroupID(dIter->first.c_str(),specName.c_str());
				Matrix* gamma_i_k_s=gammaMgr->getGamma(ogid,(string&)dIter->first,specName);
				//Need to sum over all possible ways in which we could end up in clusterID
				for(int r=0;r<gamma_i_k_s->getRowCnt();r++)
				{
					double g_i=gamma_i_k_s->getValue(r,clusterID);
					sum=sum+g_i;
					double diff=g_i*((*geneExpr)[i]-m1)*((*geneExpr)[j]-m2);
					cov=cov+diff;
				}
			}
			cov=cov/sum;

			// If we have fixed covariance, this is where we set it. (DC)
			if (fixCov)
			{
				cov=constCov;
			}
			covariance->setValue(cov,i,j);
			covariance->setValue(cov,j,i);
		}
	}
	//cout <<"Mean estimated from "<< totaldp << " in " << specName <<":" << clusterID<< endl;
//	mean->showMatrix();
	if(e->getMean()==NULL)
	{
		e->setMean(mean);
	}
	//cout <<"Covariance estimated from " << totaldp << " in " << specName << ":" << clusterID << endl;
//	covariance->showMatrix();
	if(e->getCovariance()==NULL)
	{
		e->setCovariance(covariance);
	}
	else
	{
		e->updateCovariance();
	}
	//cout <<"Estimated params for " << specName << ":"<< clusterID << " prior="<< prior << endl;
	return 0;
}

/*
 * DRMN analogue to estimateMeanCov. Estimates regulatory programs.
 * This function will update information in the Expert...
 * Actually, we need to update all the experts' datapoints before running MRTLE.
 * **SR: updated this function 
 */
/**
 * SR: This function will use MRTLE logic to find the best regulatory program for each module
	//We will follow MRTLE logic. clusterID points to the "target gene". We have a set of regulators that we will
	//iterate over to find the best regulator. We need the set of regulator IDs, which we will make the same across all
	//cell types. 
 * */
/**
* DC updating this to try to implement bookkeeping!
*/

int
SpeciesClusterManager::estimateRegProgs_PerModule(int moduleID)
{
	cout << "Updating regulatory program for " << moduleID << " ";
	time_t result = time(0);
	cout << std::asctime(localtime(&result)) << endl;

	map<string,Matrix*> cell2exp;
	map<string,Matrix*> cell2feat;
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		string s=sIter->first;
		CLUSTERSET* expertSet=sIter->second;
		Expert* expert=(*expertSet)[moduleID];
		map<string,int>& geneSet=expert->getGeneSet();

		map<string,int> temp_gene2id;
		temp_gene2id.clear();
		int tgi=0;
		for (map<string,int>::iterator gitr = geneSet.begin(); gitr!=geneSet.end(); gitr++)
		{
			string gname = gitr->first;
			temp_gene2id[gname] = tgi;
			tgi++;
		}

		int regCnt = allRegulatorIndex.size();
		int geneCnt = geneSet.size();
		cout << s << " has " << geneCnt << " genes." << endl;
		EvidenceManager* evMgr=evMgrSet[s];
		Matrix* Y = new Matrix(1,geneCnt);
		for (map<string,int>::iterator gitr = geneSet.begin(); gitr!=geneSet.end(); gitr++)
		{
			string gname = gitr->first;
			int gi=temp_gene2id[gname];
			EMAP* evidMap=evMgr->getEvidenceAt(gname);
			int varID=varNameIDMap["Expression"];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			Y->setValue(val,0,gi);
		}
		cell2exp[s] = Y;
		Matrix* X = new Matrix(regCnt,geneCnt);
		for (map<string,int>::iterator gitr = geneSet.begin(); gitr!=geneSet.end(); gitr++)
		{
			string gname = gitr->first;
			int gi = temp_gene2id[gitr->first];
			EMAP* evidMap = evMgr->getEvidenceAt(gname);
			for (map<string,int>::iterator ritr=allRegulatorIndex.begin(); ritr!=allRegulatorIndex.end();ritr++)
			{
				int ri=ritr->second;
				int varID=varNameIDMap[ritr->first];
				//Evidence* evid = (*evidMap)[varID];
				//double val = evid->getEvidVal();
				double val = 0;
				if(evidMap->find(varID)!=evidMap->end())
				{
					val=(*evidMap)[varID];
				}
				X->setValue(val,ri,gi);
			}
		}
		cell2feat[s] = X;
	}

	//cout << "Providing updated module data to expert (assign gene to expert?)" << endl;
	//As in mrtle, our score will be made up two parts: (a) how well does this regulator 
	//explain the expression of this module, and (b) how costly is it to use this regulator in terms of it being unique to 
	//a specific cell type
	//Each time don't add more than regsPerIter regulators (this is to be determined what or how)

	// regsPerIter: number of regulators that we're adding right now 
	
	double hasConv=false;
	int currIter=0;
	double currScore=0;
	SpeciesDistManager* sdMgr=gammaMgr->getSpeciesDistManager();
	//NEW: compute all covar for all mRNA and regulator features.
	// key is species, value is a base covariance matrix.
	map<string,map<int,INTDBLMAP*>> baseCovarAllSpecies; 
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		string s=sIter->first;
		Matrix* X = cell2feat[s];
		Matrix* Y = cell2exp[s];
		CLUSTERSET* expertSet=sIter->second;
		Expert* expert=(*expertSet)[moduleID];
		if (expert==NULL)
		{
			cerr << "NULL EXPERT for " << sIter->first << endl; // diagnosing badness
			//continue; // Skip this species?? no, let it seg fault!
		}

		map<int,INTDBLMAP*> myCovar;
		//computeBaseUnnormCovar(expert, s, myCovar);
		estimateCov_All(expert, X, Y, myCovar);
		baseCovarAllSpecies[s]=myCovar;
	}
	
	while(currIter<regsPerIter && !hasConv)
	{

		//Each time we want to add the best regulator. The score of the regulator includes how good the regulator 
		//explains the expression of a module and how similar are the regulators (we can tune this as an input 
		//parameter). Note, we could also use a multi-task learning framework to force selecting the same regulator
		//or use dirty group lasso. TBD	
		string bestRegulator;
		double bestScoreImprovement=0;
		map<string,int> speciesEdgeStat_bestRegulator;
		map<string,double> cellTypesToChange_bestRegulator;
		map<string,DRMNPotential*> potentialsForCellTypesToChange_bestRegulator;
		double bestPrior=0;
		string bestRegulatorModulePair;

		clock_t t=clock();
		for(map<string,int>::iterator rIter=allRegulatorSet.begin();rIter!=allRegulatorSet.end();rIter++)
		{
			//We will now go over each species' expert that corresponds to the clusterID and then we will
			//estimate the score improvement over the current regulator set; which is empty in iter=0;
			map<string,int> speciesEdgeStat;
			double currScoreImprovement=0;
			map<string,double> cellTypesToChange;
			map<string,DRMNPotential*> potentialsForCellTypesToChange;
			for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
			{
				string s=sIter->first;
				Matrix* X = cell2feat[s];
				Matrix* Y = cell2exp[s];
				CLUSTERSET* expertSet=sIter->second;
				Expert* expert=(*expertSet)[moduleID];
				if (expert==NULL)
				{
					cerr << "NULL EXPERT for " << sIter->first << endl; // diagnosing badness
					//continue; // Skip this species?? no, let it seg fault!
				}
				map<int,INTDBLMAP*> baseCovar=baseCovarAllSpecies[s];

				//We'll need to pass the predictive feature-gene associations aka motifs too here
				//SpeciesFeatureManager* sfMgr=speciesFeatSet[sIter->first];
				bool regulatorStatus=false;
				DRMNPotential* aPot=NULL;
				// DC UPDATE: try the bookkeeping version.
				//double scoreImprovement_PerCelltype=assessScoreImprovement(expert,(string&)rIter->first,(string&)sIter->first,regulatorStatus,&aPot);
				//double scoreImprovement_PerCelltype=assessScoreImprovement_Bookkeeping(expert,(string&)rIter->first,(string&)sIter->first,regulatorStatus,&aPot,baseCovar);
				double scoreImprovement_PerCelltype=assessScoreImprovement_Bookkeeping_Ali(expert,(string&)rIter->first,(string&)sIter->first,regulatorStatus,&aPot,baseCovar, X, Y);

				//Three things can happen: scoreImprovement_PerCelltype can be positive, negative or 0.
				//If it is 0, it means the regulator is already associated with this expert, but we want to
				//be doubly sure so we will pass a flag regulatorStatus that will be true if the regulator is 
				//already in this expert's regulatory program 
				//speciesEdgeStat should be set to 1
				if(regulatorStatus)
				{
					speciesEdgeStat[sIter->first]=1;
				}
				else 
				{
					if(scoreImprovement_PerCelltype>0)
					{
						currScoreImprovement=currScoreImprovement+scoreImprovement_PerCelltype;
						speciesEdgeStat[sIter->first]=1;
						cellTypesToChange[sIter->first]=scoreImprovement_PerCelltype;
						//potentialsForCellTypesToChange[sIter->first]=(*aPot);
						potentialsForCellTypesToChange[sIter->first]=aPot;
					}
					else
					{
						speciesEdgeStat[sIter->first]=0;
					}
				}
			} // end species loop

			//Now use the sharing prior to see how much this specific regulator will help
			
			double ePrior=sdMgr->getEdgeStatusProb(speciesEdgeStat,moduleID);
			char regModulePair[1024];
			sprintf(regModulePair,"%s-Module%d",rIter->first.c_str(),moduleID);
			string regModulePairKey(regModulePair); 
			double oldpriorScore=regModulePairSetPrior[regModulePairKey];
			ePrior=log(ePrior);
			currScoreImprovement=currScoreImprovement+ePrior-oldpriorScore;
			if(currScoreImprovement>0 && currScoreImprovement>bestScoreImprovement)
			{
				//cout <<"Updating best score improvement: " << currScoreImprovement <<" for edge " << regModulePairKey << endl;
				bestScoreImprovement=currScoreImprovement;
				bestRegulator.clear();
				bestRegulator.append(rIter->first);
				bestRegulatorModulePair.clear();
				bestRegulatorModulePair.append(regModulePairKey.c_str());
				bestPrior=ePrior;
				cellTypesToChange_bestRegulator.clear();
				for(map<string,double>::iterator aIter=cellTypesToChange.begin();aIter!=cellTypesToChange.end();aIter++)
				{
					cellTypesToChange_bestRegulator[aIter->first]=aIter->second;
					potentialsForCellTypesToChange_bestRegulator[aIter->first]=potentialsForCellTypesToChange[aIter->first];
				}
				for(map<string,int>::iterator sIter=speciesEdgeStat.begin();sIter!=speciesEdgeStat.end();sIter++)
				{
					speciesEdgeStat_bestRegulator[sIter->first]=sIter->second;
				}
			}
			else
			{
				cellTypesToChange.clear();
				speciesEdgeStat.clear();
				for(map<string,DRMNPotential*>::iterator aIter=potentialsForCellTypesToChange.begin();aIter!=potentialsForCellTypesToChange.end();aIter++)
				{
					delete aIter->second;
				}
				potentialsForCellTypesToChange.clear();
				speciesEdgeStat.clear();
				continue;
			}
		} // end loop over all possible regulators
	
		t=clock()-t;
		cout << "\tTime to score regulators (s) = " << ((long double) t)/CLOCKS_PER_SEC << endl;

		// Clean up the base covars (DC NEW)
		/*for (map<string,map<int,INTDBLMAP*> >::iterator bIter=baseCovarAllSpecies.begin(); bIter!=baseCovarAllSpecies.end(); bIter++)
		{
			map<int,INTDBLMAP*> gCovar=bIter->second;
			for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin() ;idIter!=gCovar.end();idIter++)
			{
				idIter->second->clear();
				delete idIter->second;
			}
		}*/

		//	if no regulator chosen, clean up
		if(bestScoreImprovement<=0)
		{
			cellTypesToChange_bestRegulator.clear();
			for(map<string,DRMNPotential*>::iterator aIter=potentialsForCellTypesToChange_bestRegulator.begin();aIter!=potentialsForCellTypesToChange_bestRegulator.end();aIter++)
			{
				delete aIter->second;
			}
			potentialsForCellTypesToChange_bestRegulator.clear();
			speciesEdgeStat_bestRegulator.clear();
		}
		else
		{	
			cout <<"Adding new regulator-module pair "<<  bestRegulatorModulePair << "\t"; //<< endl;
			// print out how long it took us to get here
			time_t result = time(0);
			cout << std::asctime(localtime(&result)) << endl;

			//At this point, we have found a module regulator configuration that we are going to go ahead with
			//So let's make the move
			for(map<string,double>::iterator aIter=cellTypesToChange_bestRegulator.begin();aIter!=cellTypesToChange_bestRegulator.end();aIter++)
			{
				CLUSTERSET* experts=speciesExpertSet[aIter->first];
				Expert* e=(*experts)[moduleID];
				e->addRegulator(bestRegulator);
				double currScore=e->getLLScore();
				e->setLLScore(currScore+aIter->second);
				DRMNPotential* aPot=e->getDRMNPotential();
				if(aPot!=NULL)
				{
					delete aPot;	
				}
				e->setDRMNPotential(potentialsForCellTypesToChange_bestRegulator[aIter->first]);
				//update the edge status now
				STRINTMAP* edgestatus=edgeCelltypeMap[bestRegulatorModulePair];
				for(map<string,int>::iterator sIter=speciesEdgeStat_bestRegulator.begin();sIter!=speciesEdgeStat_bestRegulator.end();sIter++)
				{
					(*edgestatus)[sIter->first]=sIter->second;
				}
			}
			regModulePairSetPrior[bestRegulatorModulePair]=bestPrior;
		}
		double newScore=getDRMNScoreForModule(moduleID);
		if(currIter==0)
		{
			currScore=newScore;
		}	
		else
		{
			double diff=newScore-currScore;
			currScore=newScore;	
			if(diff<0.01)//arbitrary now
			{
				hasConv=true;
			}
		}
		currIter++;
	} // end loop over reg-adding iterations for this module
	cout <<"Done learning regulators for moduleID " << moduleID <<endl;
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* experts=speciesExpertSet[sIter->first];
		Expert* e=(*experts)[moduleID];
		map<string,int>& regSet=e->getCurrentRegSet();
		if (regSet.size()==0)
		{
			cout << "No regulators for " << sIter->first << " Module" << moduleID << endl;
		}
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			cout << sIter->first<< " Module"<<moduleID<<"-" << rIter->first <<endl;
		}
	}
	for (map<string,Matrix*>::iterator sitr=cell2feat.begin(); sitr!=cell2feat.end(); sitr++)
	{
		string  s = sitr->first;
		Matrix* X = cell2feat[s];
		Matrix* Y = cell2exp[s];
		delete X;
		delete Y;
	}
	cell2feat.clear();
	cell2exp.clear();

	return 0;
}

//SR: This is probably the hardest/tricky part of all the code because we are trying to predict expression from regulator activities
//So far we had everything loaded up in EvidenceManger, but I don't think we need that. What we want is to predict mRNA of the module
//from the features.
//In order for us to do this, let's figure out what PotentialManager does for creating a potential.
//First, it associates variables with the potential
//Second, it calls potZeroInit()
//Third, it calls populatePotential().
//Fourth, it calls initMBCovMean()

//I am going to take these out of Potential.C and PotentialManager.C and try to initialize the conditional mean and variance in each expert.
double
SpeciesClusterManager::assessScoreImprovement(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential** potPtr) 
{
	//SpeciesFeatureManager* sfMgr=speciesFeatSet[cellType];
	EvidenceManager* evMgr=evMgrSet[cellType];
	map<string,int>& regSet=e->getCurrentRegSet();
	
	if(regSet.find(regName)!=regSet.end())
	{
		regStatus=true;
		return 0;
	}
	//Otherwise, suppose that this expert does not have any regulators yet. We cannot do any pre-computation upfront like 
	//the code that uses PotentialManager because the genes in the module keeps changing (unless we do soft assignment, in which case we just keep all the stuff). So how to do this.
	//Better to work with Matrix and keep populating them and destroying them.
	
	map<string,int>& geneSet=e->getGeneSet();

	// if geneset is 0, issue? this means the expert has no genes! :(
	// I think this might be happening earlier
	if (geneSet.size()==0)
	{
		cout << "geneset size 0 for expert!" << endl;
		return 0;
	}

	//First let's do a sanity check to make sure there are some genes that this regulator can regulate some genes in this set
	//int overlapCnt=sfMgr->getTargetHitCntForRegulator(regName,geneSet);
	//Must overlap at least 2% of the genes.. but might want to also know if this is a non-specific regulator.
	//double overlapFrac=((double)overlapCnt)/((double)geneSet.size());
	//double expFrac=sfMgr->getExpectedTargetCnt(regName);

	// DC -- Removing this check, but making it check for ANY overlap
	// why did I remove it? probably it was failing when I was testing.
	//if((overlapFrac/expFrac)<1) // enriched more than other modules ("spec")
	//if(overlapCnt==0) // original version
	//if(overlapFrac<0.1) // Cheap version: must overlap at least 10% of module genes. ("o10")
	//{
	//	//cout <<"Regulator " << regName <<" has too little overlap :"<< overlapFrac<<endl;
	//	return 0;
	//}
	
	/*if(strcmp(regName.c_str(),"Nfic")==0)
	{
		cout <<"Stop here" << endl;
	}*/
	//cout <<"Assessing " << regName<< " for "<< cellType << " overlapCnt: " << overlapCnt <<" among "<< geneSet.size() << endl;
	//We will add the regulator to this module's regulator set and remove it when we are done
	regSet[regName]=0;
	//We are going to use our EvidenceManager associated with this cell type to compute everything. No need to make copies of the data. Ew.
	int colID=0;
	//First we will compute the mean, and then the covariance (this is the equivalent of calling populatePotential
	//To compute the mean  we will first consider the mRNA levels 
	//Then we can calculate the coefficients using something initMBCovMean

	INTDBLMAP mean; 
	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		//first get expression
		int varID=varNameIDMap["Expression"];
		//Evidence* evid=(*evidMap)[varID];
		//double val=evid->getEvidVal();
		double val=(*evidMap)[varID];
		if(mean.find(varID)==mean.end())
		{
			mean[varID]=val;
		}
		else	
		{
			mean[varID]=mean[varID]+val;
		}
		//Now get the existing and new regulators
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			int varID=varNameIDMap[rIter->first];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			if(mean.find(varID)==mean.end())
			{
				mean[varID]=val;
			}
			else	
			{
				mean[varID]=mean[varID]+val;
			}
		}
	}
	
	//Now we just normalize to get the mean
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double val=vIter->second/geneSet.size();
		vIter->second=val;
	}
	//Mean done. We now need the covariance. 
	//We can just go over each evidence and populate the covariance entry which we will then normalize
	map<int,INTDBLMAP*> gCovar;
	int covPair=0;
	//Now the variance
	// for each gene, we will look at each of expression and regulators

	// TODO: If we are within the same DRMN iteration, we don't need to recompute the whole gCovar
	// each time we consider a regulator. 
	// We can start with expression and existing vars, and then do each variable paired with the new one.
	// So we can skip the inner loop except for one.
	// I think we also only need to normalize the new entries, since the gene set size hasn't changed.

	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++) // for each of expression and all regs
		{
			int vId=vIter->first;
			//Evidence* evid=(*evidMap)[vIter->first];
			//double vval=evid->getEvidVal();
			double vval=(*evidMap)[vIter->first];
			double vmean=mean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(INTDBLMAP_ITER uIter=vIter;uIter!=mean.end();uIter++) // for everything else
			{
				int uId=uIter->first;
				//Evidence* evid1=(*evidMap)[uIter->first];
				//double uval=evid1->getEvidVal();
				double uval=(*evidMap)[uIter->first];
				double umean=mean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					//add a small correction to avoid singularity
					(*vcov)[uId]=diffprod+0.001;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						//add a small correction to avoid singularity
						(*ucov)[vId]=diffprod+0.001;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}


	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=(0.001+vIter->second)/((double)(geneSet.size()-1));
				vIter->second=(vIter->second)/((double)(geneSet.size()-1));
				double variance=vIter->second;
			}
			else
			{
				vIter->second=vIter->second/((double)(geneSet.size()-1));
				//vIter->second=vIter->second/((double)(gCovar.size()-1));
				//vIter->second=0;
			}
		}
	}


	DRMNPotential* pot=new DRMNPotential;
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		string& vName=varIDNameMap[vIter->first];
		if(strcmp(vName.c_str(),"Expression")==0)
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::FACTOR);
		}
		else
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::MARKOV_BNKT);
		}
	}
	pot->potZeroInit();
	//populate potential using mean and covariance
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double m=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		m=mean[vIter->first];
		covar=gCovar[vIter->first];
		pot->updateMean(vIter->first,m);
		for(INTDBLMAP_ITER uIter=vIter;uIter!=mean.end();uIter++)
		{
			double cval=(*covar)[uIter->first];
			pot->updateCovariance(vIter->first,uIter->first,cval);
			pot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}
	//pot->makeValidJPD(ludecomp,perm);
	int status=pot->initMBCovMean();
	

	if(status!=0)
	{
		cout <<"Determinant too small, not going ahead with potential adding " <<regName<< " in cell "<< cellType  << endl;
		//pot->dumpPotential(cout); // DC doesn't think we need to print this right now
		map<string,int>::iterator rIter=regSet.find(regName);
		regSet.erase(rIter);
		delete pot;
		return 0;
	}

	//Now we have estimated everything! Now we can use the same logic as before to get the score
	// how slow is this part?
	double currScore=getScoreForPot(geneSet,pot,cellType);

	double oldScore=e->getLLScore();
	double impr=currScore-oldScore;
	mean.clear();
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		idIter->second->clear();
		delete idIter->second;
	}
	//potPtr=&pot;
	*potPtr=pot;
	//potPtr=pot;
	gCovar.clear();
	map<string,int>::iterator rIter=regSet.find(regName);
	regSet.erase(rIter);

	return impr;
}



/**
* DC is trying to save time by only recomputing the mean/covariance matrix for a module and its regulators
* after we change the module members. We don't need to recompute it every time we consider a new regulator.
* We will need to make the 'base' gCovar at the beginning of the iteration. This will include all existing regulators (if any).
* Then, we will start with that base gCovar, and we will only need to compute and normalize pairs of old vars with the new regulator.
*/
double
SpeciesClusterManager::assessScoreImprovement_Bookkeeping(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential** potPtr, map<int,INTDBLMAP*> baseCovar) 
{
	// fail if base covar is empty
	if (baseCovar.size()==0)
	{
		cerr << "You need to initialize the base covariance matrix first" << endl;
		return 0;
	}

	//SpeciesFeatureManager* sfMgr=speciesFeatSet[cellType];
	EvidenceManager* evMgr=evMgrSet[cellType];
	map<string,int>& regSet=e->getCurrentRegSet();
	
	if(regSet.find(regName)!=regSet.end())
	{
		regStatus=true;
		return 0;
	}
	//Otherwise, suppose that this expert does not have any regulators yet. We cannot do any pre-computation upfront like 
	//the code that uses PotentialManager because the genes in the module keeps changing (unless we do soft assignment, in which case we just keep all the stuff). So how to do this.
	//Better to work with Matrix and keep populating them and destroying them.
	
	map<string,int>& geneSet=e->getGeneSet();

	// if geneset is 0, issue? this means the expert has no genes! :(
	// I think this might be happening earlier
	if (geneSet.size()==0)
	{
		cout << "geneset size 0 for expert!" << endl;
		return 0;
	}

	//First let's do a sanity check to make sure there are some genes that this regulator can regulate some genes in this set
	//int overlapCnt=sfMgr->getTargetHitCntForRegulator(regName,geneSet);
	//Must overlap at least 2% of the genes.. but might want to also know if this is a non-specific regulator.
	//double overlapFrac=((double)overlapCnt)/((double)geneSet.size());
	//double expFrac=sfMgr->getExpectedTargetCnt(regName);

	// DC -- Removing this check, but making it check for ANY overlap
	// why did I remove it? probably it was failing when I was testing.
	//if((overlapFrac/expFrac)<1) // enriched more than other modules ("spec")
	//if(overlapCnt==0) // original version
	//if(overlapFrac<0.1) // Cheap version: must overlap at least 10% of module genes. ("o10")
	//{
	//	//cout <<"Regulator " << regName <<" has too little overlap :"<< overlapFrac<<endl;
	//	return 0;
	//}
	
	/*if(strcmp(regName.c_str(),"Nfic")==0)
	{
		cout <<"Stop here" << endl;
	}*/
	//cout <<"Assessing " << regName<< " for "<< cellType << " overlapCnt: " << overlapCnt <<" among "<< geneSet.size() << endl;
	//We will add the regulator to this module's regulator set and remove it when we are done
	regSet[regName]=0;
	//We are going to use our EvidenceManager associated with this cell type to compute everything. No need to make copies of the data. Ew.
	int colID=0;
	//First we will compute the mean, and then the covariance (this is the equivalent of calling populatePotential
	//To compute the mean  we will first consider the mRNA levels 
	//Then we can calculate the coefficients using something initMBCovMean

	INTDBLMAP mean; 
	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		// first get expression
		int varID=varNameIDMap["Expression"];
		//Evidence* evid=(*evidMap)[varID];
		//double val=evid->getEvidVal();
		double val=(*evidMap)[varID];
		if(mean.find(varID)==mean.end())
		{
			mean[varID]=val;
		}
		else	
		{
			mean[varID]=mean[varID]+val;
		}
		//Now get the existing and new regulators
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			int varID=varNameIDMap[rIter->first];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			if(mean.find(varID)==mean.end())
			{
				mean[varID]=val;
			}
			else	
			{
				mean[varID]=mean[varID]+val;
			}
		}
	}
	
	//Now we just normalize to get the mean
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double val=vIter->second/geneSet.size();
		vIter->second=val;
	}
	//Mean done. We now need the covariance. 
	//We can just go over each evidence and populate the covariance entry which we will then normalize

	// NEW STUFF: We'll copy in the base unnormalized covariance matrix.
	map<int,INTDBLMAP*> gCovar;
	copyCovar(baseCovar, gCovar); // copy in the base covariance matrix.

	// debugging to make sure it worked. I'm satisfied.
	//printCovar(baseCovar);
	//printCovar(gCovar);

	int covPair=0;
	//Now the variance
	// for each gene, we will look at each of expression and regulators

	// TODO: If we are within the same DRMN iteration, we don't need to recompute the whole gCovar
	// each time we consider a regulator. 
	// We can start with expression and existing vars, and then do each variable paired with the new one.
	// So we can skip the inner loop except for one.
	// I think we also only need to normalize the new entries, since the gene set size hasn't changed.
	int newRegVar=varNameIDMap[regName]; 


	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		int vId=newRegVar;

		//Evidence* evid=(*evidMap)[vId];
		//double vval=evid->getEvidVal();
		double vval=(*evidMap)[vId];
		double vmean=mean[vId];
		INTDBLMAP* vcov=NULL;
		if(gCovar.find(vId)==gCovar.end()) 
		{
			vcov=new INTDBLMAP;
			gCovar[vId]=vcov;
		}
		else
		{
			vcov=gCovar[vId]; 
		}

		for(INTDBLMAP_ITER uIter=mean.begin();uIter!=mean.end();uIter++) // For the other variables...
		{
			int uId=uIter->first;
			//if (uId==vId) // don't compare to itself? no, do!
			//{
			//	continue; 
			//}
			//Evidence* evid1=(*evidMap)[uIter->first];
			//double uval=evid1->getEvidVal();
			double uval=(*evidMap)[uIter->first];
			double umean=mean[uId];
			double diffprod=(vval-vmean)*(uval-umean);
			INTDBLMAP* ucov=NULL;
			if(gCovar.find(uId)==gCovar.end())
			{
				ucov=new INTDBLMAP;
				gCovar[uId]=ucov;
			}
			else
			{
				ucov=gCovar[uId];
			}
			if(vcov->find(uId)==vcov->end())
			{
				covPair++;
				//add a small correction to avoid singularity
				(*vcov)[uId]=diffprod+0.001;
			}
			else
			{
				(*vcov)[uId]=(*vcov)[uId]+diffprod;
			}
			if(uId!=vId)
			{
				if(ucov->find(vId)==ucov->end())
				{
					//add a small correction to avoid singularity
					(*ucov)[vId]=diffprod+0.001;
				}
				else
				{
					(*ucov)[vId]=(*ucov)[vId]+diffprod;
				}
			}	
		}
	}
	
	//printCovar(gCovar);

	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=(0.001+vIter->second)/((double)(geneSet.size()-1));
				vIter->second=(vIter->second)/((double)(geneSet.size()-1));
				double variance=vIter->second;
			}
			else
			{
				vIter->second=vIter->second/((double)(geneSet.size()-1));
				//vIter->second=vIter->second/((double)(gCovar.size()-1));
				//vIter->second=0;
			}
		}
	}
	//cout << "Normalized covar:" << endl;
	//printCovar(gCovar);


	DRMNPotential* pot=new DRMNPotential;
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		string& vName=varIDNameMap[vIter->first];
		if(strcmp(vName.c_str(),"Expression")==0)
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::FACTOR);
		}
		else
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::MARKOV_BNKT);
		}
	}
	pot->potZeroInit();
	//populate potential using mean and covariance
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double m=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		m=mean[vIter->first];
		covar=gCovar[vIter->first];
		pot->updateMean(vIter->first,m);
		for(INTDBLMAP_ITER uIter=vIter;uIter!=mean.end();uIter++)
		{
			double cval=(*covar)[uIter->first];
			pot->updateCovariance(vIter->first,uIter->first,cval);
			pot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}
	//pot->makeValidJPD(ludecomp,perm);
	int status=pot->initMBCovMean();
	

	if(status!=0)
	{
		cout <<"Determinant too small, not going ahead with potential adding " <<regName<< " in cell "<< cellType  << endl;
		//pot->dumpPotential(cout); // DC doesn't think we need to print this right now
		map<string,int>::iterator rIter=regSet.find(regName);
		regSet.erase(rIter);
		delete pot;
		return 0;
	}

	//Now we have estimated everything! Now we can use the same logic as before to get the score
	double currScore=getScoreForPot(geneSet,pot,cellType);
	double oldScore=e->getLLScore();
	double impr=currScore-oldScore;
	mean.clear();
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		idIter->second->clear();
		delete idIter->second;
	}
	//potPtr=&pot;
	*potPtr=pot;
	//potPtr=pot;
	gCovar.clear();
	map<string,int>::iterator rIter=regSet.find(regName);
	regSet.erase(rIter);

	return impr;
}

double
SpeciesClusterManager::assessScoreImprovement_Bookkeeping_Ali(Expert* e,string& regName,string& cellType,bool& regStatus,DRMNPotential** potPtr, map<int,INTDBLMAP*> baseCovar, Matrix* X, Matrix* Y) 
{
	// fail if base covar is empty
	if (baseCovar.size()==0)
	{
		cerr << "You need to initialize the base covariance matrix first" << endl;
		return 0;
	}

	//SpeciesFeatureManager* sfMgr=speciesFeatSet[cellType];
	EvidenceManager* evMgr=evMgrSet[cellType];
	map<string,int>& regSet=e->getCurrentRegSet();
	
	if(regSet.find(regName)!=regSet.end())
	{
		regStatus=true;
		return 0;
	}
	//Otherwise, suppose that this expert does not have any regulators yet. We cannot do any pre-computation upfront like 
	//the code that uses PotentialManager because the genes in the module keeps changing (unless we do soft assignment, in which case we just keep all the stuff). So how to do this.
	//Better to work with Matrix and keep populating them and destroying them.
	
	map<string,int>& geneSet=e->getGeneSet();

	// if geneset is 0, issue? this means the expert has no genes! :(
	// I think this might be happening earlier
	if (geneSet.size()==0)
	{
		cout << "geneset size 0 for expert!" << endl;
		return 0;
	}

	//First let's do a sanity check to make sure there are some genes that this regulator can regulate some genes in this set
	//int overlapCnt=sfMgr->getTargetHitCntForRegulator(regName,geneSet);
	//Must overlap at least 2% of the genes.. but might want to also know if this is a non-specific regulator.
	//double overlapFrac=((double)overlapCnt)/((double)geneSet.size());
	//double expFrac=sfMgr->getExpectedTargetCnt(regName);

	// DC -- Removing this check, but making it check for ANY overlap
	// why did I remove it? probably it was failing when I was testing.
	//if((overlapFrac/expFrac)<1) // enriched more than other modules ("spec")
	//if(overlapCnt==0) // original version
	//if(overlapFrac<0.1) // Cheap version: must overlap at least 10% of module genes. ("o10")
	//{
	//	//cout <<"Regulator " << regName <<" has too little overlap :"<< overlapFrac<<endl;
	//	return 0;
	//}
	
	/*if(strcmp(regName.c_str(),"Nfic")==0)
	{
		cout <<"Stop here" << endl;
	}*/
	//cout <<"Assessing " << regName<< " for "<< cellType << " overlapCnt: " << overlapCnt <<" among "<< geneSet.size() << endl;
	//We will add the regulator to this module's regulator set and remove it when we are done
	regSet[regName]=0;
	//We are going to use our EvidenceManager associated with this cell type to compute everything. No need to make copies of the data. Ew.
	int colID=0;
	//First we will compute the mean, and then the covariance (this is the equivalent of calling populatePotential
	//To compute the mean  we will first consider the mRNA levels 
	//Then we can calculate the coefficients using something initMBCovMean

	INTDBLMAP mean; 
	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		// first get expression
		int varID=varNameIDMap["Expression"];
		//Evidence* evid=(*evidMap)[varID];
		//double val=evid->getEvidVal();
		double val=(*evidMap)[varID];
		if(mean.find(varID)==mean.end())
		{
			mean[varID]=val;
		}
		else	
		{
			mean[varID]=mean[varID]+val;
		}
		//Now get the existing and new regulators
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			int varID=varNameIDMap[rIter->first];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			if(mean.find(varID)==mean.end())
			{
				mean[varID]=val;
			}
			else	
			{
				mean[varID]=mean[varID]+val;
			}
		}
	}
	
	//Now we just normalize to get the mean
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double val=vIter->second/geneSet.size();
		vIter->second=val;
	}
	//Mean done. We now need the covariance. 
	//We can just go over each evidence and populate the covariance entry which we will then normalize

	// NEW STUFF: We'll copy in the base unnormalized covariance matrix.
	map<int,INTDBLMAP*> gCovar;
	//copyCovar(baseCovar, gCovar); // copy in the base covariance matrix.
	//addCov(e, X, Y, gCovar, regName);

	// debugging to make sure it worked. I'm satisfied.
	//printCovar(baseCovar);
	//printCovar(gCovar);

	//Now the variance
	// for each gene, we will look at each of expression and regulators

	// TODO: If we are within the same DRMN iteration, we don't need to recompute the whole gCovar
	// each time we consider a regulator. 
	// We can start with expression and existing vars, and then do each variable paired with the new one.
	// So we can skip the inner loop except for one.
	// I think we also only need to normalize the new entries, since the gene set size hasn't changed.
	
	//printCovar(gCovar);

	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	//cout << "Normalized covar:" << endl;
	//printCovar(gCovar);


	DRMNPotential* pot=new DRMNPotential;
	DRMNPotential* parentPot=new DRMNPotential;
	int expVarID=-1;
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		string& vName=varIDNameMap[vIter->first];
		if(strcmp(vName.c_str(),"Expression")==0)
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::FACTOR);
			expVarID=vIter->first;
		}
		else
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::MARKOV_BNKT);
			parentPot->setAssocVariable(vIter->first,DRMNPotential::MARKOV_BNKT);
		}
	}
	pot->potZeroInit();
	parentPot->potZeroInit();
	//populate potential using mean and covariance
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double m=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		m=mean[vIter->first];
		//covar=gCovar[vIter->first];
		covar=baseCovar[vIter->first];
		pot->updateMean(vIter->first,m);
		if(vIter->first!=expVarID)
		{
			parentPot->updateMean(vIter->first,m);
		}
		for(INTDBLMAP_ITER uIter=vIter;uIter!=mean.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				cout <<"Did not find covariance entry for " << varIDNameMap[uIter->first] << " in " << varIDNameMap[vIter->first] << endl;
				exit(0);
			}
			double cval=(*covar)[uIter->first];
			pot->updateCovariance(vIter->first,uIter->first,cval);
			pot->updateCovariance(uIter->first,vIter->first,cval);
			if(uIter->first!=expVarID && vIter->first!=expVarID)
			{
				parentPot->updateCovariance(vIter->first,uIter->first,cval);
				parentPot->updateCovariance(uIter->first,vIter->first,cval);
			}
		}
	}
	//pot->makeValidJPD(ludecomp,perm);
	int status=pot->initMBCovMean();
	

	if(status!=0)
	{
		cout <<"Determinant too small, not going ahead with potential adding " <<regName<< " in cell "<< cellType  << endl;
		//pot->dumpPotential(cout); // DC doesn't think we need to print this right now
		map<string,int>::iterator rIter=regSet.find(regName);
		regSet.erase(rIter);
		delete pot;
		return 0;
	}

	//Now we have estimated everything! Now we can use the same logic as before to get the score
	//double currScore=getScoreForPot(geneSet,pot,cellType);
	double currScore=getScoreForPot_Tracetrick(geneSet,pot,parentPot,cellType);
	double oldScore=e->getLLScore();
	double impr=currScore-oldScore;
	mean.clear();
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		idIter->second->clear();
		delete idIter->second;
	}
	//potPtr=&pot;
	*potPtr=pot;
	//potPtr=pot;
	gCovar.clear();
	map<string,int>::iterator rIter=regSet.find(regName);
	regSet.erase(rIter);

	return impr;
}

/**
* Copies the contents of one covariance matrix into another.
* We assume the target is empty.
* We will use this to book keep during consideration of candidate regulators.
*/
int
SpeciesClusterManager::copyCovar(map<int,INTDBLMAP*>& source, map<int,INTDBLMAP*>& target)
{
	// loop over the first int...
	for (map<int,INTDBLMAP*>::iterator sIter=source.begin(); sIter!=source.end(); sIter++)
	{
		int i=sIter->first;
		INTDBLMAP* insides=sIter->second;

		target[i]=new INTDBLMAP; // store the first key
		for (INTDBLMAP_ITER jIter=insides->begin(); jIter!=insides->end(); jIter++)	
		{
			int j=jIter->first;
			double v=jIter->second;
			(*target[i])[j]=v;
		}
	}
	return 0;
}

/**
* Prints a covariance matrix to the console (for debugging)
*/
int
SpeciesClusterManager::printCovar(map<int,INTDBLMAP*>& gCovar)
{
	// print a header row by variable name
	cout << "Covar";
	for (map<int,INTDBLMAP*>::iterator sIter=gCovar.begin(); sIter!=gCovar.end(); sIter++)
	{
		cout << "\t" << varIDNameMap[sIter->first];
	}
	cout << endl;

	for (map<int,INTDBLMAP*>::iterator sIter=gCovar.begin(); sIter!=gCovar.end(); sIter++)
	{
		int i=sIter->first;
		cout << varIDNameMap[i];
		INTDBLMAP* insides=sIter->second;
		for (INTDBLMAP_ITER jIter=insides->begin(); jIter!=insides->end(); jIter++)	
		{
			int j=jIter->first;
			double v=jIter->second;
			cout << "\t" << v;
		}
		cout << endl;
	}
	return 0;
}

/**
* This computes the unnormalized covariance matrix for an existing module and set of regulators.
* We will use this to avoid recomputing it every time we consider a new regulator.
*/
int
SpeciesClusterManager::computeBaseUnnormCovar(Expert* e,string& cellType, map<int,INTDBLMAP*>& gCovar) 
{
	//SpeciesFeatureManager* sfMgr=speciesFeatSet[cellType];
	EvidenceManager* evMgr=evMgrSet[cellType];
	map<string,int>& regSet=e->getCurrentRegSet(); // current regulators

	// get genes
	map<string,int>& geneSet=e->getGeneSet();

	// compute mean -- I haven't implemented bookkeeping with this yet
	INTDBLMAP mean;
	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		// first get expression
		int varID=varNameIDMap["Expression"];
		//Evidence* evid=(*evidMap)[varID];
		//double val=evid->getEvidVal();
		double val=(*evidMap)[varID];
		if(mean.find(varID)==mean.end())
		{
			mean[varID]=val;
		}
		else	
		{
			mean[varID]=mean[varID]+val;
		}
		//Now get the existing regulators
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			int varID=varNameIDMap[rIter->first];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			if(mean.find(varID)==mean.end())
			{
				mean[varID]=val;
			}
			else	
			{
				mean[varID]=mean[varID]+val;
			}
		}
	}
	//Now we just normalize to get the mean
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double val=vIter->second/geneSet.size();
		vIter->second=val;
	}
	
	// if geneset is 0, issue? this means the expert has no genes! :(
	// I think this might be happening earlier
	if (geneSet.size()==0)
	{
		cout << "geneset size 0 for expert!" << endl;
		return 0;
	}


	int covPair=0;
	//Now the variance
	// for each gene, we will look at each of expression and regulators

	// TODO: If we are within the same DRMN iteration, we don't need to recompute the whole gCovar
	// each time we consider a regulator. 
	// We can start with expression and existing vars, and then do each variable paired with the new one.
	// So we can skip the inner loop except for one.
	// I think we also only need to normalize the new entries, since the gene set size hasn't changed.

	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++) // for each of expression and all regs
		{
			int vId=vIter->first;
			//Evidence* evid=(*evidMap)[vIter->first];
			//double vval=evid->getEvidVal();
			double vval=(*evidMap)[vIter->first];
			double vmean=mean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(INTDBLMAP_ITER uIter=vIter;uIter!=mean.end();uIter++) // for itself AND everything else
			{
				int uId=uIter->first;
				//Evidence* evid1=(*evidMap)[uIter->first];
				//double uval=evid1->getEvidVal();
				double uval=(*evidMap)[uIter->first];
				double umean=mean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					//add a small correction to avoid singularity
					(*vcov)[uId]=diffprod+0.001;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						//add a small correction to avoid singularity
						(*ucov)[vId]=diffprod+0.001;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}
	return 0;
}

//This could be optimized by using the covariance trick.
double
SpeciesClusterManager::getScoreForPot(map<string,int>& geneSet,DRMNPotential* apot, string& cellType)
{
	//Each row is a sample/gene
	double score=0;
	EvidenceManager* evMgr=evMgrSet[cellType];
	for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)gIter->first);
		double ll=apot->getCondPotValueFor(evidMap);
		score=score+log(ll);
	}
	double varCnt=apot->getAssocVariables().size();
	double paramCnt=2*varCnt;
	//double paramCnt=varCnt;
	paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
	double p_score=score-((paramCnt/2)*log(geneSet.size()));
	//pll=unreg_pll-(0.5*log(evMgr->getTestSet().size()));
	//return pll;
	double score_tt=apot->computeLL_Tracetrick(geneSet.size());
	return score;
	//return p_score;
}

//This could be optimized by using the covariance trick.
double
SpeciesClusterManager::getScoreForPot_Tracetrick(map<string,int>& geneSet,DRMNPotential* apot, DRMNPotential* parentPot, string& cellType)
{
	//Each row is a sample/gene
	/*double score=0;
	EvidenceManager* evMgr=evMgrSet[cellType];
	for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)gIter->first);
		double ll=apot->getCondPotValueFor(evidMap);
		score=score+log(ll);
	}
	double varCnt=apot->getAssocVariables().size();
	double paramCnt=2*varCnt;
	//double paramCnt=varCnt;
	paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
	double p_score=score-((paramCnt/2)*log(geneSet.size()));*/
	//pll=unreg_pll-(0.5*log(evMgr->getTestSet().size()));
	//return pll;
	double score_tt1=apot->computeLL_Tracetrick(geneSet.size());
	double score_tt2=parentPot->computeLL_Tracetrick(geneSet.size());
	double score_tt=score_tt1-score_tt2;
	//return score;
	return score_tt;
	//return p_score;
}

/*
 * EMINT
 * Here we are computing the emission probabilities for each gene, each expert for one species.
 * For DRMN, this will be largely the same -- maybe exactly the same? The difference is what's
 * going on inside the expert and how it computes the outputPDF.
 */
int
SpeciesClusterManager::expectationStep_Species(string& specName, CLUSTERSET* expertSet)
{
	GeneExpManager* exprMgr=speciesExprSet[specName];
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
	{
		vector<double>* exprProf=exprMgr->getExp(vIter->first);
		map<int,double> mixOutProbs;
		double sum=0;
		if(strcmp(vIter->first.c_str(),"orf19.5809")==0 || strcmp(vIter->first.c_str(),"Q0080")==0)
		{
			cout << "Found " << vIter->first << endl;
		}
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			Expert* e=eIter->second;
			double pdf=e->getOutputPDF(exprProf);
			if(std::isnan(pdf))
			{
				cout <<"PDF is nan for " << specName << " " << vIter->first << " for expert " << eIter->first << endl;
			}
			mixOutProbs[eIter->first]=pdf;
			sum=sum+pdf;
		}
		int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),specName.c_str());
		//gammaMgr->estimateLeafGamma(ogid,mixOutProbs,(string&)vIter->first,specName);
		gammaMgr->setObservedProbs(ogid,mixOutProbs,(string&)vIter->first,specName);
		
 	       mixOutProbs.clear();
	}
	return 0;
}


//DRMN expectation step will be similar to CMINT except we have new potentials. What is not clear is whether we do expectation
//at individual nodes or do this jointly at all nodes. And do we assign genes now to experts or do we maintain soft assignments
//The way the maximization step is set up is assuming hard assignment.
int
SpeciesClusterManager::expectationStep_DRMN_Species(string& specName, CLUSTERSET* expertSet)
{
	GeneExpManager* exprMgr=speciesExprSet[specName];
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	EvidenceManager* evMgr=evMgrSet[specName];
	
	for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
	{
		EMAP* evidSet=evMgr->getEvidenceAt((string&)vIter->first);
		map<int,double> mixOutProbs;
		double sum=0;
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			Expert* e=eIter->second;
			double pdf=e->getDRMNProb(evidSet);
			if(std::isnan(pdf))
			{
				cout <<"PDF is nan for " << specName << " " << vIter->first << " for expert " << eIter->first << endl;
			}
			mixOutProbs[eIter->first]=pdf;
			sum=sum+pdf;
		}
		int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),specName.c_str());
		//gammaMgr->estimateLeafGamma(ogid,mixOutProbs,(string&)vIter->first,specName);
		gammaMgr->setObservedProbs(ogid,mixOutProbs,(string&)vIter->first,specName);
		
 	       mixOutProbs.clear();
	}
	return 0;
}

int
SpeciesClusterManager::doInferenceForGene(string& geneName, string& speciesName)
{
	EvidenceManager* evMgr=evMgrSet[speciesName];
	EMAP* evidence=evMgr->getEvidenceAt(geneName);
	if(evidence==NULL)
	{
		cout <<"Did not find gene " << geneName << "for species" << speciesName << endl;
	}
	CLUSTERSET* expertSet=speciesExpertSet[speciesName];
	map<int,double> mixOutProbs;
	double sum=0;
	for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
	{
		Expert* e=eIter->second;
		double pdf=e->getDRMNProb(evidence);
		if(std::isnan(pdf))
		{
			cout <<"PDF is nan for expert " << eIter->first << endl;
		}
		mixOutProbs[eIter->first]=pdf;
		sum=sum+pdf;
	}
	int ogid=mor->getMappedOrthogroupID(geneName.c_str(),speciesName.c_str());
	//gammaMgr->estimateLeafGamma(ogid,mixOutProbs,(string&)vIter->first,specName);
	gammaMgr->setObservedProbs(ogid,mixOutProbs,geneName,speciesName);
 	mixOutProbs.clear();
	return 0;
}


int
SpeciesClusterManager::showClusters_Extant(const char* outputDir)
{
	char aFName[1024];
	sprintf(aFName,"%s/rseed.txt",outputDir);
	ofstream oFile(aFName);
	oFile << rseed << endl;
	oFile.close();
	assignGenesToExperts_FromMap();
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		char outFName[1024];
		sprintf(outFName,"%s/%s_exprtab.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusters(outFName,(string&)sIter->first);
		sprintf(outFName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusterAssignments(outFName,(string&)sIter->first);
		sprintf(outFName,"%s/%s_speciesspecnames_clusterassign.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusterAssignments_NonSrcNames(outFName,(string&)sIter->first);
	}
	return 0;
}

int
SpeciesClusterManager::showMeans(const char* outputDir)
{
	char outFName[1024];
	sprintf(outFName,"%s/clustermeans.txt",outputDir);
	ofstream oFile(outFName);
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* clustering=sIter->second;
		for(CLUSTERSET_ITER cIter=clustering->begin();cIter!=clustering->end();cIter++)
		{
			oFile <<sIter->first<<"_"<<cIter->first;
			Expert* expert=cIter->second;
			Matrix* m=expert->getMean();
			for(int c=0;c<m->getColCnt();c++)
			{
				oFile <<"\t" << m->getValue(0,c);
			}
			oFile << endl;
		}
	}
	oFile.close();
	return 0;
}


int
SpeciesClusterManager::showMeans(const char* outputDir, int iter)
{
	char outFName[1024];
	sprintf(outFName,"%s/clustermeans_%d.txt",outputDir,iter);
	ofstream oFile(outFName);
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* clustering=sIter->second;
		for(CLUSTERSET_ITER cIter=clustering->begin();cIter!=clustering->end();cIter++)
		{
			oFile <<sIter->first<<"_"<<cIter->first;
			Expert* expert=cIter->second;
			Matrix* m=expert->getMean();
			for(int c=0;c<m->getColCnt();c++)
			{
				oFile <<"\t" << m->getValue(0,c);
			}
			oFile << endl;
		}
	}
	oFile.close();
	return 0;
}

int
SpeciesClusterManager::showClusters(const char* outputDir)
{
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		char outFName[1024];
		sprintf(outFName,"%s/%s_exprtab.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusters(outFName,(string&)sIter->first);
		sprintf(outFName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusterAssignments(outFName,(string&)sIter->first);
	}
	return 0;
}




int
SpeciesClusterManager::showClusters_Ancestral(const char* outputDir)
{
	map<string,ofstream*> filePtrs;
	char geneName[32];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		map<string,int>* clusterAssign=ogIter->second;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			strcpy(geneName,cIter->first.c_str());
			char* pos=strchr(geneName,':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			*pos='\0';
			string specName(pos+1);
			if(speciesExpertSet.find(specName)!=speciesExpertSet.end())
			{	
				continue;
			}
			ofstream* oFile=NULL;
			if(filePtrs.find(specName)==filePtrs.end())
			{	
				char aFName[1024];
				sprintf(aFName,"%s/%s_clusterassign.txt",outputDir,specName.c_str());
				oFile=new ofstream(aFName);
				filePtrs[specName]=oFile;
			}
			else
			{	
				oFile=filePtrs[specName];
			}
			GeneMap* geneMap=mgrp->getSpeciesHits(srcSpecies);
			if(geneMap==NULL)
			{
				(*oFile)<< "OG"<<ogIter->first<<":"<< geneName <<"\t" << cIter->second << endl;
			}
			else
			{
				//Get the scer ortholog of this gene
				map<string,map<string,STRINTMAP*>*>& geneset=geneMap->getGeneSet();
				//Fear here is that the same scer gene is repeated within and between clusters
				for(map<string,map<string,STRINTMAP*>*>::iterator sIter=geneset.begin();sIter!=geneset.end();sIter++)
				{
					(*oFile) <<sIter->first <<"\t" << cIter->second << endl;
				}
			}
		}
	}
	return 0;
}


int
SpeciesClusterManager::displaySpeciesClusters(const char* outFName,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	GeneExpManager* expMgr=speciesExprSet[specName];
	ofstream oFile(outFName);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	vector<string>& colNames=expMgr->getColNames();
	oFile <<"Loci";
	for(int c=0;c<colNames.size();c++)
	{
		oFile<<"\t" <<colNames[c];
	}
	oFile <<endl;
	vector<double>* expr=expMgr->getExp(speciesGenes->begin()->first);
	int dim=expr->size();
	/*for(int i=0;i<dim;i++)
	{
		oFile <<"\tExp"<<i;
	}
	oFile << endl;*/
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			for(int i=0;i<dim;i++)
			{
				oFile <<"\t" <<"-100";
			}
			oFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			oFile<< gIter->first;
			vector<double>* expr=expMgr->getExp(gIter->first);
			if(expr==NULL)
			{
				cout <<"No expression for " << gIter->first << endl;
				continue;
			}
			for(int i=0;i<expr->size();i++)
			{
				oFile << "\t"<<(*expr)[i];
			}
			oFile << endl;
		}
		oFile <<"Dummy" << clusterID;
		for(int i=0;i<dim;i++)
		{
			oFile <<"\t" <<"-100";
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;
}

/*
 * Prints regulator and target files for this cluster, this species (species specific names)
 * Will be used for 
 */
/*int
SpeciesClusterManager::printRegTargetFiles(const char* outfile, int cID, string& specName)
{
	ofstream ofile(outfile);
	
	SpeciesFeatureManager* specMgr=speciesFeatSet[specName];
	

	ofile.close();
}*/


/*
* From MRTLE.
* Reads regulator OG IDs.
*/
int
SpeciesClusterManager::setRestrictedList(const char* aFName)
{
	ifstream inFile(aFName);
	string buffer;
	int r=0;
	while(inFile.good())
	{
		getline(inFile,buffer);
		if(buffer.length()<=0)
		{
			continue;
		}
		int ogid=atoi(buffer.c_str());
		inputRegulatorOGs[ogid]=0;
		r++;
	}
	inFile.close();

	if (r==0) // no regulators?
	{
		return 1;
	}

	return 0;
}

/*
* Identifies the OGID that we'll use for "expression" in MRTLE. Starting at 0, use first unused OGID.
*/
int 
SpeciesClusterManager::setUpModuleOGIDs()
{
	// We'll let an unused OGID be the expression data and give speciesname_clusterID as the "gene names".
	int og=0;

	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	while (allOgs.find(og)!=allOgs.end())
	{
		og++; 
	}
	moduleOG=og; // this is the member variable
	//cout << "Using og " << moduleOG << " for expression." << endl;
	//cout << "TO DO: Add module gene names to MOR!!!" << endl;
	return 0;
}

/*
*
* Create ogids file for the regulatory network data. Do this once per cluster just to make
* the results easier to interpret.
* It's like we're rotating the matrix: the module is the target gene and the regulatory features 
* are other genes.
* NOT USED
*/
int
SpeciesClusterManager::makeOGIDsForRegnet(const char* outDir, vector<string>& speciesList)
{
	
	// Now print out the regulator orthogroups as we read them in
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	
	

	// Make a target file (really just need one)
	char outf[1024];
	sprintf(outf, "%s/all_target_ids.txt", outDir);
	ofstream oFile(outf);
	oFile << moduleOG << endl;
	oFile.close();
	
	// now for each cluster
	for (int cID=0; cID<maxClusterCnt;cID++)
	{
	    	char outpref[1024];    // prefix for cluster-specific files
		sprintf(outpref, "%s/state%d", outDir, cID);

		char outf[1024];
		sprintf(outf, "%s_ortho_map.txt", outpref);
		cout << outf << endl;
		ofstream oFile(outf); 

		oFile << "Gene_OGID\tNAME" << endl;
		oFile << "OG" << moduleOG << "_0\t";
		for(int s=0;s<speciesList.size();s++)
		{
			oFile << speciesList[s] << "_state" << cID;
			if (s<speciesList.size()-1)
			{
				oFile << ",";
			}
		}
		oFile << endl;

		//for (map<int, MappedOrthogroup*>::iterator rIter=allOgs.begin(); rIter!=allOgs.end(); rIter++)
		for (map<int,int>::iterator rIter=inputRegulatorOGs.begin(); rIter!=inputRegulatorOGs.end();rIter++)
		{
			// get mapped orthogroup for this ID
			MappedOrthogroup* group=allOgs[rIter->first];
		
			oFile << "OG" << group->getID() << "_0\t";
			//tfile << rIter->second->getID() << endl; // tf file
			for (int s=0; s<speciesList.size(); s++)
			{
				GeneMap* geneMap=group->getSpeciesHits(speciesList[s].c_str());
				map<string,map<string,STRINTMAP*>*>& geneSet=geneMap->getGeneSet();


		                for(map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		                {
		                        if(gIter==geneSet.begin())
		                        {
		                                oFile << gIter->first;
		                        }
		                        else
		                        {
		                                oFile <<","<< gIter->first;
		                        }
		                }			


				if (s<speciesList.size()-1)
				{
					oFile <<",";
				}
			}
			oFile << endl;

		}
		oFile.close();

	}
	return 0;
}

/*
* Loop over cell types and dump cluster-specific data.
*/
int
SpeciesClusterManager::dumpFeatures(const char* outDir, vector<string>& speciesList)
{
	for(int s=0;s<speciesList.size();s++)
	{
		dumpFeaturesPerCluster(outDir, speciesList[s]);
	}
	return 0;
}

/*
* DC added
* For a species, dumps out a new feature file per cluster that incorporates both expression and feature data.
* Format:
*	gene	Expr	motif1	motif2 .... motifN
* For now, we'll print "NaN" if no value given for motif to gene in the network data.
*
* outDir : name of directory to drop the output files (one file per cluster)
* specName : species name 
*
* Also print target list (name of expression column)
*/
int
SpeciesClusterManager::dumpFeaturesPerCluster(const char* outDir,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	GeneExpManager* expMgr=speciesExprSet[specName];
	SpeciesFeatureManager* featMgr=speciesFeatSet[specName];

	// Get the possible regulators in the motif network
	map<string,int> regulators;
	featMgr->getRegulators(regulators);

	
	// one file per cluster for features, one for "targets", one for regulators
	for (int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		char outFName[1024];
		sprintf(outFName,"%s/%s_features_state%d.txt",outDir, specName.c_str(), clusterID);
		ofstream oFile(outFName); 
		cout << outFName << endl;

		// target file for cell type, cluster
		char myname[1024];
		sprintf(myname, "%s/%s_targets_state%d.txt", outDir, specName.c_str(), clusterID);
		cout << myname << endl;
		ofstream tfile(myname);
		tfile << specName << "_state" << clusterID << endl;
		tfile.close();	

		// regulator file file for cell type, cluster
		sprintf(myname, "%s/%s_regulators_state%d.txt", outDir, specName.c_str(), clusterID);
		cout << myname << endl;
		ofstream rfile(myname);
		
		// we'll get the genes from the experts
		CLUSTERSET* expertSet=speciesExpertSet[specName];
		vector<string>& colNames=expMgr->getColNames();

		oFile <<"Loci"; // Header column

		// First do "expression" column(s)
		// UPDATE: 
		if (colNames.size()>1)
		{
			cerr << "Currently only implemented for 1 expression value per cell type!!!" << endl;
		}
		for(int c=0;c<colNames.size();c++) // this version is if we have different names for expression
		{
			oFile<<"\t" <<colNames[c] << "_state" << clusterID;
		}



		// next do feature data
		// also make regulator file
		for(map<string,int>::iterator regIter=regulators.begin(); regIter!=regulators.end(); regIter++)
		{
			oFile<<"\t" << regIter->first;
			rfile << regIter->first << endl;
		}
		rfile.close(); // end regulator list for this cluster, cell type
		oFile <<endl;
		

		// and now per gene...
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			oFile<< gIter->first;
			vector<double>* expr=expMgr->getExp(gIter->first);
			if(expr==NULL)
			{
				cout <<"No expression for " << gIter->first << endl;
				continue;
			}
			for(int i=0;i<expr->size();i++)
			{
				oFile << "\t"<<(*expr)[i];
			}

			// Now the cisreg elements
			for(map<string,int>::iterator regIter=regulators.begin(); regIter!=regulators.end(); regIter++)
			{
				string reg=regIter->first;
				string target=gIter->first;
				//cout << reg << "\t" << target << endl;
				if (featMgr->hasValue(reg,target))
				{
					oFile<<"\t" << featMgr->getValue(reg,target);
				}
				else
				{
					oFile<<"\t"<< featMgr->getDefaultValue();
				}
			} 

			oFile << endl;
		}
		oFile.close();	
	}

	
	return 0;
	
	/*CLUSTERSET* expertSet=speciesExpertSet[specName];
	vector<string>& colNames=expMgr->getColNames();
	oFile <<"Loci";
	for(int c=0;c<colNames.size();c++)
	{
		oFile<<"\t" <<colNames[c];
	}
	oFile <<endl;
	vector<double>* expr=expMgr->getExp(speciesGenes->begin()->first);
	int dim=expr->size();
	/*for(int i=0;i<dim;i++)
	{
		oFile <<"\tExp"<<i;
	}
	oFile << endl;*/
	/*map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			for(int i=0;i<dim;i++)
			{
				oFile <<"\t" <<"-100";
			}
			oFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			oFile<< gIter->first;
			vector<double>* expr=expMgr->getExp(gIter->first);
			if(expr==NULL)
			{
				cout <<"No expression for " << gIter->first << endl;
				continue;
			}
			for(int i=0;i<expr->size();i++)
			{
				oFile << "\t"<<(*expr)[i];
			}
			oFile << endl;
		}
		oFile <<"Dummy" << clusterID;
		for(int i=0;i<dim;i++)
		{
			oFile <<"\t" <<"-100";
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;*/
}

/*
int
SpeciesClusterManager::dumpAllInferredClusters_Srcwise(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_clusterassign_brk.txt",outputDir);
	ofstream oFile(outputFName);
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	ofstream dFile(outputFName);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	//for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	map<string,int> shownGenes;
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			dFile <<"Dummy" << clusterID;
			for(int s=0;s<speciesList.size();s++)
			{
				oFile <<"\t" <<"-100";
				dFile <<"\t" <<"-100";
			}
			oFile <<endl;
			dFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{	
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			int ogid=mgrp->getID();
			if(ogid==3865)
			{
				cout <<"OG " << ogid << endl;
			}
			map<string,int>* clusterAssign=mappedClusterAssignment[ogid];
			map<string,int> geneAssign;
			char geneName[1024];
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				strcpy(geneName,cIter->first.c_str());
				char* pos=strchr(geneName,':');
				if(pos==NULL)
				{
					cout <<"Bad format " << endl;	
					exit(0);
				}
				*pos='\0';
				string specName(pos+1);
				string geneNameKey(geneName);
				if(strstr(specName.c_str(),"Anc")==NULL)
				{
					geneAssign[geneNameKey]=cIter->second;
				}
				else	
				{
					geneAssign[geneNameKey]=cIter->second;
				}
			}
			const char* dupAnc=gammaMgr->getDupAncestor(ogid);
			bool wgd=false;
			if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
			{
				wgd=true;
			}
			map<int,map<string,string>*>& geneSets=mgrp->getGeneSets();
			map<string,int> localKey;
			for(map<int,map<string,string>*>::iterator hIter=geneSets.begin();hIter!=geneSets.end();hIter++)
			{
				map<string,string>* gset=hIter->second;
				if(gset->find(srcSpecies)==gset->end())
				{
					oFile << "OG" << ogid <<"_" << hIter->first;
					if(wgd)
					{
						oFile <<"*";
					}
					if(geneSets.size()>1)
					{
						dFile << "OG" << ogid <<"_" << hIter->first;
						if(wgd)
						{
							dFile <<"*";
						}
					}
				}	
				else
				{
					string dispName(gnm.getCommonName((*gset)[srcSpecies].c_str()));
					if(localKey.find(dispName)!=localKey.end())
					{
						dispName.append("_2");
					}
					localKey[dispName]=0;
					oFile << dispName;
					if(wgd)
					{
						oFile <<"*";
					}
					if(geneSets.size()>1)
					{
						dFile << dispName;
						if(wgd)
						{
							dFile <<"*";
						}
					}
				}
				for(int s=0;s<speciesList.size();s++)
				{
					int geneclusterID=-2;
					if(strstr(speciesList[s].c_str(),"Anc")!=NULL)
					{
						char temp[1024];
						sprintf(temp,"%s_%d",speciesList[s].c_str(),hIter->first+1);
						string tempKey(temp);
						//if(geneAssign.find(speciesList[s])!=geneAssign.end())
						if(geneAssign.find(tempKey)!=geneAssign.end())
						{
							geneclusterID=geneAssign[tempKey];
						}
								
					}
					else
					{
						if(gset->find(speciesList[s])==gset->end())
						{
							geneclusterID=-2;
						}
						else 
						{	
							string& geneName=(*gset)[speciesList[s]];
							if(geneAssign.find(geneName)!=geneAssign.end())
							{
								geneclusterID=geneAssign[geneName];
							}
							shownGenes[geneName]=0;
						}
					}
					oFile<< "\t"<< geneclusterID;
					if(geneSets.size()>1)
					{
						dFile<< "\t"<< geneclusterID;
					}
				}
				oFile << endl;	
				if(geneSets.size()>1)
				{
					dFile<< endl;
				}
			}
			localKey.clear();
		}
		oFile <<"Dummy" << clusterID;
		dFile <<"Dummy" << clusterID;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
			dFile <<"\t" <<"-100";
		}
		oFile <<endl;
		dFile <<endl;
	}
	oFile.close();
	dFile.close();
	return 0;
}*/


int
SpeciesClusterManager::dumpAllInferredClusters_SrcwiseGrouped(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_clusterassign_brk.txt",outputDir);
	ofstream oFile(outputFName);
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	ofstream dFile(outputFName);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	map<int,map<string,int>*> allClusterAssignmentsGrouped;
	gammaMgr->getAllClusterAssignments_Grouped(allClusterAssignmentsGrouped);
	map<int,string> ogDupSpeciesMap;
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			dFile <<"Dummy" << clusterID;
			for(int s=0;s<speciesList.size();s++)
			{
				oFile <<"\t" <<"-100";
				dFile <<"\t" <<"-100";
			}
			oFile <<endl;
			dFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{	
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			int ogid=mgrp->getID();
			map<string,int>* clusterAssign=allClusterAssignmentsGrouped[ogid];
			map<string,int>* clusterAssignOrig=mappedClusterAssignment[ogid];
			if(ogid==4597)
			{
				cout <<"OG " << ogid << endl;
				for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
				{
					cout << cIter->first << " " << cIter->second << endl;
				}
				cout <<"Original assignment" << endl;
				for(map<string,int>::iterator cIter=clusterAssignOrig->begin();cIter!=clusterAssignOrig->end();cIter++)
				{
					cout << cIter->first << " " << cIter->second << endl;
				}
			}

			int groupID=0;
			map<int,map<string,int>*> geneAssignSet;
			map<int,map<string,string>*> geneNamesSet;
			char geneName[1024];
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				char tempBuff[1024];
				strcpy(tempBuff,cIter->first.c_str());
				char* tok=strtok(tempBuff,":");
				int tokCnt=0;
				string speciesName;
				string aName;
				int dup=-1;
				while(tok!=NULL)
				{
					if(tokCnt==0)
					{
						dup=atoi(tok);
					}
					else if(tokCnt==1)
					{
						speciesName.append(tok);
					}	
					else
					{
						aName.append(tok);
					}
					tok=strtok(NULL,":");
					tokCnt++;
				}
				map<string,int>* geneSet=NULL;
				map<string,string>* geneName=NULL;
				if(geneAssignSet.find(dup)==geneAssignSet.end())
				{
					geneSet=new map<string,int>;
					geneName=new map<string,string>;
					geneAssignSet[dup]=geneSet;
					geneNamesSet[dup]=geneName;	
				}
				else
				{
					geneSet=geneAssignSet[dup];
					geneName=geneNamesSet[dup];
				}
				(*geneSet)[speciesName]=cIter->second;
				(*geneName)[speciesName]=aName;
				char buffer[1024];
				string key;
				if(strstr(cIter->first.c_str(),"Anc")!=NULL)
				{
					sprintf(buffer,"%s_%d:%s",speciesName.c_str(),dup,speciesName.c_str());
				}	
				else
				{
					sprintf(buffer,"%s:%s",aName.c_str(),speciesName.c_str());
				}
				key.append(buffer);
				int otherassign=-1;
				if(clusterAssignOrig->find(key)==clusterAssignOrig->end())
				{
					cout << " No cluster for "<< ogid << " " << key.c_str() << endl;
				}
				else
				{
					otherassign=(*clusterAssignOrig)[key];
					if(otherassign!=cIter->second)
					{
						cout <<"Cluster asssignment mismatch for  " << ogid<<  " " << key.c_str() << " OLD " << otherassign << " " << cIter->second<< endl;
					}
				}
			}
			const char* dupAnc=gammaMgr->getDupAncestor(ogid);
			bool wgd=false;
			if(dupAnc!=NULL)
			{
				string dupAncStr(dupAnc);
				ogDupSpeciesMap[ogid]=dupAncStr;	
			}
			if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
			{
				wgd=true;
			}
			map<int,map<string,string>*>& geneSets=mgrp->getGeneSets();
			if(geneSets.size()!=geneAssignSet.size())
			{
				cout <<"Missed a duplicate for  " << ogid << endl;
			}
			map<string,int> localKey;
			map<string,int> localKeyDup;
			for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
			{
				map<string,int>* gset=hIter->second;
				map<string,string>* gname=geneNamesSet[hIter->first];
				int srcClusterAssign=-3;
				if(gset->find(srcSpecies)!=gset->end())
				{
					srcClusterAssign=(*gset)[srcSpecies];
				}
				
				if(gset->find(srcSpecies)==gset->end())
				{
					if(srcClusterAssign==clusterID)
					{
						oFile << "OG" << ogid <<"_" << hIter->first;
						if(wgd)
						{
							oFile <<"*";
						}
					}
					//if(geneSets.size()>1)
					if(geneAssignSet.size()>1)
					{
						dFile << "OG" << ogid <<"_" << hIter->first;
						if(wgd)
						{
							dFile <<"*";
						}
					}
				}	
				else
				{
					if(srcClusterAssign==clusterID)
					{
						string& geneName=(*gname)[srcSpecies];
						shownGenes[geneName]=0;
						string dispName(gnm.getCommonName((*gname)[srcSpecies].c_str()));
						if(localKey.find(dispName)!=localKey.end())
						{
							dispName.append("_2");
						}
						oFile << dispName;
						if(wgd)
						{
							oFile <<"*";
						}
						localKey[dispName]=0;
					}
					//if(geneSets.size()>1)
					if(geneAssignSet.size()>1)
					{
						string dispName(gnm.getCommonName((*gname)[srcSpecies].c_str()));
						if(localKeyDup.find(dispName)!=localKeyDup.end())
						{
							dispName.append("_2");
						}
						dFile << dispName;
						if(wgd)
						{
							dFile <<"*";
						}
						localKeyDup[dispName]=0;
					}
				}

				for(int s=0;s<speciesList.size();s++)
				{
					int geneclusterID=-2;
					if(gset->find(speciesList[s])==gset->end())
					{
						geneclusterID=-2;
					}
					else 
					{	
						string& geneName=(*gname)[speciesList[s]];
						geneclusterID=(*gset)[speciesList[s]];
				//		shownGenes[geneName]=0;
					}
					if(srcClusterAssign==clusterID)
					{
						oFile<< "\t"<< geneclusterID;
					}
					//if(geneSets.size()>1)
					if(geneAssignSet.size()>1)
					{
						dFile<< "\t"<< geneclusterID;
					}
				}
				if(srcClusterAssign==clusterID)
				{
					oFile << endl;	
				}
				//if(geneSets.size()>1)
				if(geneAssignSet.size()>1)
				{
					dFile<< endl;
				}
			}
		}
		oFile <<"Dummy" << clusterID;
		dFile <<"Dummy" << clusterID;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
			dFile <<"\t" <<"-100";
		}
		oFile <<endl;
		dFile <<endl;
	}
	oFile.close();
	dFile.close();
	sprintf(outputFName,"%s/dupinfo_inferred.txt",outputDir);
	ofstream dupFile(outputFName);
	for(map<int,string>::iterator ogIter=ogDupSpeciesMap.begin();ogIter!=ogDupSpeciesMap.end();ogIter++)
	{
		dupFile << ogIter->first <<"\t" << ogIter->second << endl;
	}
	dupFile.close();
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusters_LCA(const char* outputDir,vector<string>& speciesList, string& lcaName)
{	
	char outputFName[1024];
	//sprintf(outputFName,"%s/allspecies_clusterassign_lca_brk.txt",outputDir);
	sprintf(outputFName,"%s/allcelltypes_clusterassign_brk.txt",outputDir);
	ofstream oFile(outputFName);
	//sprintf(outputFName,"%s/allspecies_duplicates_clusterassign__lcabrk.txt",outputDir);
	//ofstream dFile(outputFName);
	//sprintf(outputFName,"%s/genemembers_perog.txt",outputDir);
	//ofstream geneMembersFile(outputFName);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	map<int,map<string,int>*> allClusterAssignmentsGrouped;
	gammaMgr->getAllClusterAssignments_Grouped(allClusterAssignmentsGrouped);
	map<int,string> ogDupSpeciesMap;
	int badOG_DupLCA=0;
	int badOG_NoLCA=0;
	map<int,int> shownOGs;
	oFile <<"Loci";
	for(int i=0;i<speciesList.size();i++)
	{
		oFile << "\t"<< speciesList[i];
	}
	oFile <<endl;
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		//Iterate over the allClusterAssignments but show only those groups for which the Anc14 cluster assingment matches
		//the clusterID
		for(map<int,map<string,int>*>::iterator ogIter=allClusterAssignmentsGrouped.begin();ogIter!=allClusterAssignmentsGrouped.end();ogIter++)
		{
			int ogid=ogIter->first;
			map<string,int>* clusterAssign=ogIter->second;
			map<string,int>* clusterAssignOrig=mappedClusterAssignment[ogid];
			int groupID=0;
			map<int,map<string,int>*> geneAssignSet;
			map<int,map<string,string>*> geneNamesSet;
			char geneName[1024];
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				char tempBuff[1024];
				strcpy(tempBuff,cIter->first.c_str());
				char* tok=strtok(tempBuff,":");
				int tokCnt=0;
				string speciesName;
				string aName;
				int dup=-1;
				while(tok!=NULL)
				{
					if(tokCnt==0)
					{
						dup=atoi(tok);
					}
					else if(tokCnt==1)
					{
						speciesName.append(tok);
					}	
					else
					{
						aName.append(tok);
					}
					tok=strtok(NULL,":");
					tokCnt++;
				}
				map<string,int>* geneSet=NULL;
				map<string,string>* geneName=NULL;
				if(geneAssignSet.find(dup)==geneAssignSet.end())
				{
					geneSet=new map<string,int>;
					geneName=new map<string,string>;
					geneAssignSet[dup]=geneSet;
					geneNamesSet[dup]=geneName;	
				}
				else
				{
					geneSet=geneAssignSet[dup];
					geneName=geneNamesSet[dup];
				}
				(*geneSet)[speciesName]=cIter->second;
				(*geneName)[speciesName]=aName;
				char buffer[1024];
				string key;
				if(strstr(cIter->first.c_str(),"Anc")!=NULL)
				{
					sprintf(buffer,"%s_%d:%s",speciesName.c_str(),dup,speciesName.c_str());
				}	
				else
				{
					sprintf(buffer,"%s:%s",aName.c_str(),speciesName.c_str());
				}
				key.append(buffer);
				int otherassign=-1;
				if(clusterAssignOrig->find(key)==clusterAssignOrig->end())
				{
					cout << " No cluster for "<< ogid << " " << key.c_str() << endl;
				}
				else
				{
					otherassign=(*clusterAssignOrig)[key];
					if(otherassign!=cIter->second)
					{
						cout <<"Cluster asssignment mismatch for  " << ogid<<  " " << key.c_str() << " OLD " << otherassign << " " << cIter->second<< endl;
					}
				}
			}
			const char* dupAnc=gammaMgr->getDupAncestor(ogid);
			bool wgd=false;
			if(dupAnc!=NULL)
			{
				string dupAncStr(dupAnc);
				ogDupSpeciesMap[ogid]=dupAncStr;	
			}
			if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
			{
				wgd=true;
			}
			map<string,int> localKey;
			map<string,int> localKeyDup;
			//Also check if Anc14 has more than 1 copy
			//Show the orthogroup if and only if one of the copies has the cluster ID as current cluster ID.
			bool toShowOG=false;
			for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
			{
				map<string,int>* gset=hIter->second;
				int ancClusterAssign=-3;
				if(gset->find(lcaName)!=gset->end())
				{
					ancClusterAssign=(*gset)[lcaName];
					if(ancClusterAssign==clusterID)
					{
						toShowOG=true;
					}
				}
				
			}
			if(!toShowOG)
			{
				continue;
			}
			shownOGs[ogid]=0;
			int lcacnt=0;
			for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
			{
				map<string,string>* gname=geneNamesSet[hIter->first];
				map<string,int>* gset=hIter->second;
				int ancClusterAssign=-3;
				if(gset->find(lcaName)!=gset->end())
				{
					ancClusterAssign=(*gset)[lcaName];
					lcacnt++;
				}
				//Now use the scername to show this og
				string dispName;
				if(gset->find(srcSpecies)==gset->end())
				{
					char tempName[1024];	
					sprintf(tempName,"OG%d_%d",ogid,hIter->first);
					dispName.append(tempName);
					if(wgd)
					{
						dispName.append("*");
					}
				}	
				else
				{
					string geneName((*gname)[srcSpecies]);
					shownGenes[geneName]=0;
					dispName.append(gnm.getCommonName((*gname)[srcSpecies].c_str()));
					if(localKey.find(dispName)!=localKey.end())
					{
						dispName.append("_2");
					}
					localKey[dispName]=0;
					if(wgd)
					{
						dispName.append("*");
					}
				}
				//if(ancClusterAssign==clusterID)
				//{
					oFile <<dispName;
					/*if(geneAssignSet.size()>1)
					{
						dFile << dispName;
					}*/
					//Show the gene members
				//}
				string members;

				for(int s=0;s<speciesList.size();s++)
				{
					int geneclusterID=-2;
					if(gset->find(speciesList[s])==gset->end())
					{
						geneclusterID=-2;
					}
					else 
					{	
						geneclusterID=(*gset)[speciesList[s]];
					}
					if(members.length()>0)
					{
						members.append(";");
					}	
					if(gname->find(speciesList[s])==gname->end())
					{
						members.append("NULL");
					}
					else
					{
						members.append((*gname)[speciesList[s]]);
					}		
					//if(ancClusterAssign==clusterID)
				//	{
						oFile<< "\t"<< geneclusterID;
				//	}
					/*if(geneAssignSet.size()>1)
					{
						dFile<< "\t"<< geneclusterID;
					}*/
				}
				//geneMembersFile <<"OG"<<ogid <<"_" << hIter->first <<"\t" << members << endl;
				//if(ancClusterAssign==clusterID)
				//{
					oFile << endl;	
				//}
				/*if(geneAssignSet.size()>1)
				{
					dFile<< endl;
				}*/
			}
			if(lcacnt>1)
			{
				badOG_DupLCA++;
				cout <<"LCA has multiple copies in "<< ogid << endl;
			}
		}
		oFile <<"Dummy" << clusterID;
	//	dFile <<"Dummy" << clusterID;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
	//		dFile <<"\t" <<"-100";
		}
		oFile <<endl;
	//	dFile <<endl;
	}
	oFile.close();
	//dFile.close();
	//geneMembersFile.close();
	/*sprintf(outputFName,"%s/dupinfo_inferred.txt",outputDir);
	ofstream dupFile(outputFName);
	for(map<int,string>::iterator ogIter=ogDupSpeciesMap.begin();ogIter!=ogDupSpeciesMap.end();ogIter++)
	{
		dupFile << ogIter->first <<"\t" << ogIter->second << endl;
	}
	dupFile.close();
	cout <<" BadOG: Dup_LCA_OGs " << badOG_DupLCA  << endl;*/
	cout <<"Shown " << shownOGs.size () << " out of " << allClusterAssignmentsGrouped.size() <<  " total OGs"  << endl;
	return 0;
}




/*
int
SpeciesClusterManager::dumpAllInferredClusters_Srcwise_Dup(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	ofstream oFile(outputFName);
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			if(mgrp->getCnt()<2)
			{
				continue;
			}
			map<string,int>* clusterAssign=mappedClusterAssignment[ogid];
			GeneMap* dupspecies_GeneMap=mgrp->getSpeciesHits(srcName);
			if(dupspecies_GeneMap==NULL)
			{
				cout <<"Warning! No " << srcSpecies << " gene in OGid "<<  << endl;
				continue;
			}
			bool dupInScer=false;
			//get the number of species that have 2 genes
			map<string,int> speciesWithDupgenes;
			for(map<string,int>::iterator aIter=speciesList.begin();aIter!=speciesList.end();aIter++)
			{
				GeneMap* specRep=mgrp->getSpeciesHits(aIter->first.c_str());
				if(specRep==NULL)
				{
					continue;
				}	
				if(specRep->getGeneSet().size()>=2)
				{
					speciesWithDupGenes[aIter->first]=0;
				}
			}
			if(speciesWithDupGenes.find(srcSpecies)!=speciesWithDupGenes.end())
			{
				dupInScer=true;
			}
			string dupSpecies(srcName);
			if(!dupInScer)
			{
				dupSpecies.clear();
				dupSpecies.append(speciesWithDupgenes.begin()->first.c_str());
				dupspecies_GeneMap=mgrp->getSpeciesHits(dupSpecies.c_str());
			}
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				strcpy(geneName,cIter->first.c_str());
				char* pos=strchr(geneName,':');
				if(pos==NULL)
				{
					cout <<"Bad format " << endl;	
					exit(0);
				}
				*pos='\0';
				string specName(pos+1);
				geneAssign[geneName]=cIter->second;
			}
			map<string,map<string,STRINTMAP*>*>& dupGenes=dupSpecies_GeneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator hIter=dupGenes.begin();hIter!=dupGenes.end();hIter++)
			{
				oFile << gnm.getCommonName(hIter->first.c_str());
				shownGenes[hIter->first]=0;
				map<string,int> specAssign;
				int clusterID=-1;
				if(geneAssign.find(hIter->first)!=geneAssign.end())
				{
					clusterID=geneAssign[hIter->first];
				}
				specAssign[srcSpecies]=clusterID;
				map<string,STRINTMAP*>* hitsInOthers=hIter->second;
				for(int s=0;s<speciesList.size();s++)
				{
					if(strcmp(speciesLis[s].c_str(),"Scer")==0)
					{
						continue;
					}
					//This means the gene is altogether missing in the species
					int clusterID=-2;
					if(hitsInOthers->find(speciesList[s])==hitsInOthers->end())
					{
						clusterID=-2;
					}
					else
					{
						map<string,int>* genesInOthers=(*hitsInOthers)[speciesList[s]];
						if(genesInOthers->size()>1)
						{
							cout <<"Warning! Multiple hits to " << hIter->first << " from " << speciesList[s] << endl;
						}
						string& orthogene=genesInOthers->begin();
						if(geneAssign.find(orthogene)!=geneAssignAssign.end())
						{	
							clusterID=geneAssign[orthogene];
						}	
						else
						{
							clusterID=-1;
						}
						shownGenes[orthogene]=1;
					}
					specAssign[speciesList[s]]=clusterID;
				}
				for(map<string,int>::iterator sIter=specAssign.begin();sIter!=specAssign.end();sIter++)
				{
					oFile <<"\t" << sIter->second;
				}
				oFile << endl;	
				specAssign.clear();
			}
			
			geneAssign.clear();
		}
		oFile <<"Dummy" << cIter->first;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;
}*/


int
SpeciesClusterManager::displaySpeciesClusterAssignments(const char* outFName,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	ofstream oFile(outFName);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),specName.c_str());
			const string& genename=gIter->first;
			if(strcmp(specName.c_str(),srcSpecies)==0)
			{
				oFile<< genename <<"\t" << cIter->first << endl;
			}
			else
			{
				//Get the scer ortholog of this gene
				map<string,int>* scerOrtholog=mor->getOrtholog(specName.c_str(),genename.c_str(),srcSpecies);
				if(scerOrtholog==NULL)
				{
					oFile <<"OG"<< mgrp->getID() << ":"<< genename<<"\t"<< cIter->first << endl;
				}
				else
				{
					for(map<string,int>::iterator sIter=scerOrtholog->begin();sIter!=scerOrtholog->end();sIter++)
					{
						oFile <<sIter->first<<"\t"<< cIter->first << endl;
					}
				}
			}
		}
	}
	oFile.close();
	return 0;
}

//The main differenc here and the displaySpeciesClusterAssignments
//is that here we display the original species-specific name of the species
int
SpeciesClusterManager::displaySpeciesClusterAssignments_NonSrcNames(const char* outFName,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	ofstream oFile(outFName);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			oFile<< gIter->first<<"\t" << cIter->first << endl;
		}
	}
	oFile.close();
	return 0;
}


int
SpeciesClusterManager::assignGenesToExperts()
{
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* expertSet=sIter->second;
		for(CLUSTERSET_ITER eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetAssignedGenes();
		}
	}
	for(map<string,GeneExpManager*>::iterator sIter=speciesExprSet.begin();sIter!=speciesExprSet.end();sIter++)
	{
		map<string,int>* speciesGenes=speciesClusterSet_Genewise[sIter->first];
		CLUSTERSET* expertSet=speciesExpertSet[sIter->first];	
		for(map<string,int>::iterator gIter=speciesGenes->begin();gIter!=speciesGenes->end();gIter++)
		{
			int ogid=mor->getMappedOrthogroupID(gIter->first.c_str(),sIter->first.c_str());
			//This gamma matrix is a kXk matrix which stores the joint probability of being in cluster i
			//given its ancestor was in some other species. 
			Matrix* normterm=gammaMgr->getNormTerm(ogid,(string&)gIter->first,(string&)sIter->first);
			int maxClusterID=-1;	
			double maxProb=0;
			for(int c=0;c<normterm->getColCnt();c++)
			{
				double prob=normterm->getValue(0,c);
				if(prob> maxProb)
				{
					maxProb=prob;
					maxClusterID=c;
				}
			}
			if(maxClusterID==-1)
			{
				cout <<"Gene " << gIter->first << " cannot be assigned to any cluster of " << sIter->first.c_str() << endl;
				continue;
			}
			Expert* e=(*expertSet)[maxClusterID];
			e->assignGeneToExpert(gIter->first.c_str());
		}	
	}
	return 0;
}


int
SpeciesClusterManager::assignGenesToExperts_FromMap()
{
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* expertSet=sIter->second;
		for(CLUSTERSET_ITER eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetAssignedGenes();
		}
	}
	char geneName[1024];
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		map<string,int>* geneAssignments=ogIter->second;
		int ogid=ogIter->first;
		for(map<string,int>::iterator gIter=geneAssignments->begin();gIter!=geneAssignments->end();gIter++)
		{
			strcpy(geneName,gIter->first.c_str());
			char* pos=strchr(geneName,':');
			if(pos==NULL)
			{
				cout <<"Bad format" << endl;
				exit(0);
			}
			if(ogid==4597)
			{
				cout << gIter->first << "\t" << gIter->second << endl;
			}
			*pos='\0';
			if(strcmp(geneName,"orf19.5809")==0)
			{	
	//			cout << "Stop here " << geneName << " " << gIter->first << " " << gIter->second<< endl;
			}
			string specName(pos+1);
			if(speciesExpertSet.find(specName)==speciesExpertSet.end())
			{
				continue;
			}
			CLUSTERSET* expertSet=speciesExpertSet[specName];
			if(expertSet->find(gIter->second)==expertSet->end())
			{
				cout << "No cluster " << gIter->second<< " for " << gIter->first  << endl;
				continue;
			}
			Expert* e=(*expertSet)[gIter->second];
			e->assignGeneToExpert(geneName);
		}	
	}
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir)
{
	char outputFName[1024];
	//sprintf(outputFName,"%s/clusterassign_multspecies.txt",outputDir);
	//ofstream oFile(outputFName);
	//gammaMgr->getAllClusterAssignments_Conditional(mappedClusterAssignment);
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,true);
	gammaMgr->showClusterFlipCnts();
	gammaMgr->reestimateTransitionProbability();
	/*char outputDirCmd[1024];
	sprintf(outputDirCmd,"mkdir %s/prior_pp/",outputDir);
	system(outputDirCmd);
	char outputSubDir[1024];
	sprintf(outputSubDir,"%s/prior_pp",outputDir);
	showClusters_Extant(outputSubDir);
	int iter=0;
	double currScore=0;
	bool convergence=false;
	//while(iter<50 && !convergence)
	while(iter<2 && !convergence)
	{
		showMeans(outputDir,iter);
		if(iter>0)
		{
			expectationStep();
		}
		maximizationStep();
		double newScore=getScore();
		double diff=fabs(newScore-currScore);
		if((iter>0) && (diff<0.5))
		{
			convergence=true;
		}
		cout <<"PPIter: " << iter << " score: " << newScore << " diff " << diff << endl;
		currScore=newScore;
		iter++;
	}*/
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,false);
	//gammaMgr->showClusterFlipCnts();
	gammaMgr->reestimateTransitionProbability();
	
	//gammaMgr->getAllClusterAssignments_Conditional(mappedClusterAssignment);
	/*map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		//oFile<< ogIter->first;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		GeneMap* geneMap=mgrp->getSpeciesHits(srcSpecies);
		if(geneMap==NULL)
		{
			oFile <<" --";
		}
		else
		{
			map<string,map<string,STRINTMAP*>*>& geneSet=geneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				if(gIter==geneSet.begin())
				{
					oFile << gIter->first;
				}
				else
				{
					oFile <<","<< gIter->first;
				}
			}
		}
		map<string,int>* clusterAss=ogIter->second;
		for(map<string,int>::iterator cIter=clusterAss->begin();cIter!=clusterAss->end();cIter++)
		{
			const char* pos=strchr(cIter->first.c_str(),':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			string specName(pos+1);
			oFile << "\t"<<specName<<"=" << cIter->second;
		}
		oFile << endl;
	}
	oFile.close();*/
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir,int clusterID)
{
	char outputFName[1024];
	sprintf(outputFName,"%s/clusterassign_multspecies_%d.txt",outputDir,clusterID);
	ofstream oFile(outputFName);
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,true);
	//gammaMgr->getAllClusterAssignments_Conditional(mappedClusterAssignment);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		//oFile<< ogIter->first;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		GeneMap* geneMap=mgrp->getSpeciesHits(srcSpecies);
		if(geneMap==NULL)
		{
			oFile <<" --";
		}
		else
		{
			map<string,map<string,STRINTMAP*>*>& geneSet=geneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				if(gIter==geneSet.begin())
				{
					oFile << gIter->first;
				}
				else
				{
					oFile <<","<< gIter->first;
				}
			}
		}
		map<string,int>* clusterAss=ogIter->second;
		for(map<string,int>::iterator cIter=clusterAss->begin();cIter!=clusterAss->end();cIter++)
		{
			const char* pos=strchr(cIter->first.c_str(),':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			string specName(pos+1);
			oFile << "\t"<<specName<<"=" << cIter->second;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

//Show all the cluster assignments per ancestral species too
int
SpeciesClusterManager::dumpAllInferredClusterGammas(const char* outputDir,vector<string>& speciesList)
{
	char outputFName[1024];
	gammaMgr->getAllClusterGammas(mappedClusterGamma);
	sprintf(outputFName,"%s/clustergamma_multspecies.txt",outputDir);
	ofstream oFile(outputFName);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	oFile <<"Loci";
	for(int s=0;s<speciesList.size();s++)
	{
		CLUSTERSET* speciesExperts=speciesExpertSet[speciesList[s]];
		for(CLUSTERSET_ITER eIter=speciesExperts->begin();eIter!=speciesExperts->end();eIter++)
		{
			oFile<< "\t" <<speciesList[s]<<"_"<< eIter->first;
		}	
	}
	oFile <<endl;
	for(map<int,map<string,map<int,double>*>*>::iterator ogIter=mappedClusterGamma.begin();ogIter!=mappedClusterGamma.end();ogIter++)
	{
		//oFile<< ogIter->first;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		int ogid=mgrp->getID();
		map<string,map<int,double>*> geneAssign;
		map<string,map<int,double>*>* clusterAssign=mappedClusterGamma[ogid];
		char geneName[1024];
		for(map<string,map<int,double>*>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			strcpy(geneName,cIter->first.c_str());
			char* pos=strchr(geneName,':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			*pos='\0';
			string specName(pos+1);
			string geneNameKey(geneName);
			if(strstr(specName.c_str(),"Anc")==NULL)
			{
				geneAssign[geneNameKey]=cIter->second;
			}
			else	
			{
				geneAssign[specName]=cIter->second;
			}
		}
		map<int,map<string,string>*>& geneSets=mgrp->getGeneSets();
		for(map<int,map<string,string>*>::iterator hIter=geneSets.begin();hIter!=geneSets.end();hIter++)
		{
			map<string,string>* gset=hIter->second;
			if(gset->find(srcSpecies)==gset->end())
			{
				continue;
			}
			oFile << gnm.getCommonName((*gset)[srcSpecies].c_str());
			for(int s=0;s<speciesList.size();s++)
			{
				map<int,double>* geneclusterID=NULL;
				if(strstr(speciesList[s].c_str(),"Anc")!=NULL)
				{
					if(geneAssign.find(speciesList[s])!=geneAssign.end())
					{
						geneclusterID=geneAssign[speciesList[s]];
					}
								
				}
				else
				{
					if(gset->find(speciesList[s])==gset->end())
					{
						geneclusterID=NULL;
					}
					else 
					{	
						string& geneName=(*gset)[speciesList[s]];
						if(geneAssign.find(geneName)!=geneAssign.end())
						{
							geneclusterID=geneAssign[geneName];
						}
					}
				}
				if(geneclusterID!=NULL)
				{
					for(map<int,double>::iterator aIter=geneclusterID->begin();aIter!=geneclusterID->end();aIter++)
					{
						oFile <<"\t"<< aIter->second;
					}
				}
				else
				{
					for(int k=0;k<maxClusterCnt;k++)
					{
						oFile <<"\t-1";
					}
				}
			}
			oFile << endl;
		}
	}
	return 0;
}


map<string,int>* 
SpeciesClusterManager::getGenesForSpecies(string& specName)
{
	if(speciesClusterSet_Genewise.find(specName)==speciesClusterSet_Genewise.end())
	{	
		return NULL;
	}
	return speciesClusterSet_Genewise[specName];
}

int
SpeciesClusterManager::generateData(const char* outputDir, string& lcaName,vector<string>& speciesList)
{
	map<int,map<string,int>*> clusterAssignments;
	gammaMgr->sampleAllClusterAssignments(clusterAssignments);
	char fName[1024];
	char dirName[1024];
	map<string,ofstream*> filePtrSet;
	map<string,ofstream*> filePtrClusteredSet;
	map<string,ofstream*> clusterFilePtrSet;
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		sprintf(fName,"%s/%s_samples.txt",outputDir,sIter->first.c_str());
		ofstream* oFile=new ofstream(fName);
		filePtrSet[sIter->first]=oFile;
		GeneExpManager* expMgr=speciesExprSet[sIter->first];
		vector<string>& colNames=expMgr->getColNames();
		(*oFile) <<"Loci";
		for(int c=0;c<colNames.size();c++)
		{
			(*oFile) <<"\t" <<colNames[c];
		}
		(*oFile) <<endl;
		sprintf(fName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
		ofstream* cFile=new ofstream(fName);
		clusterFilePtrSet[sIter->first]=cFile;
	}
	vector<double> sampleValues;
	char clusterAssignmentFName[1024];
	char allassignfname[1024];
	sprintf(allassignfname,"%s/allcelltypes_clusterassign.txt",outputDir);
	ofstream allAssignFile(allassignfname);
	allAssignFile <<"Loci";
	for(int s=0;s<speciesList.size();s++)
	{
		allAssignFile <<"\t" << speciesList[s];
	}
	allAssignFile <<endl;
	map<string,int> shownGenes;
	for(map<int,map<string,int>*>::iterator gIter=clusterAssignments.begin();gIter!=clusterAssignments.end();gIter++)
	{
		
		int ogid=gIter->first;
		map<string,int>* assign=gIter->second;
		map<int,map<string,int>*> geneAssignSet;
		map<int,map<string,string>*> geneNamesSet;
		for(map<string,int>::iterator cIter=assign->begin();cIter!=assign->end();cIter++)
		{
			char tempBuff[1024];
			strcpy(tempBuff,cIter->first.c_str());
			char* tok=strtok(tempBuff,":");
			int tokCnt=0;
			string speciesName;
			string aName;
			int dup=-1;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					dup=atoi(tok);
				}
				else if(tokCnt==1)
				{
					speciesName.append(tok);
				}	
				else
				{
					aName.append(tok);
				}
				tok=strtok(NULL,":");
				tokCnt++;
			}
			map<string,int>* geneSet=NULL;
			map<string,string>* geneName=NULL;
			if(geneAssignSet.find(dup)==geneAssignSet.end())
			{
				geneSet=new map<string,int>;
				geneName=new map<string,string>;
				geneAssignSet[dup]=geneSet;
				geneNamesSet[dup]=geneName;	
			}
			else
			{
				geneSet=geneAssignSet[dup];
				geneName=geneNamesSet[dup];
			}
			(*geneSet)[speciesName]=cIter->second;
			(*geneName)[speciesName]=aName;
		}
		for(map<int,map<string,int>*>::iterator gIter=geneAssignSet.begin();gIter!=geneAssignSet.end();gIter++)
		{
			map<string,int>* assignment=gIter->second;
			map<string,string>* geneName=geneNamesSet[gIter->first];
			for(map<string,int>::iterator sIter=assignment->begin();sIter!=assignment->end();sIter++)
			{
				if(clusterFilePtrSet.find(sIter->first)==clusterFilePtrSet.end())
				{
					sprintf(fName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
					ofstream* cFile=new ofstream(fName);
					clusterFilePtrSet[sIter->first]=cFile;
				}
				ofstream* cFile=clusterFilePtrSet[sIter->first];
				string& aName=(*geneName)[sIter->first];
				if(strstr(sIter->first.c_str(),"Anc")!=NULL)
				{
					aName=(*geneName)[srcSpecies];
				}
				(*cFile) <<aName <<"\t" << sIter->second<<endl;
				if(speciesExpertSet.find(sIter->first)==speciesExpertSet.end())
				{
					continue;	
				}
				CLUSTERSET* clusterSet=speciesExpertSet[sIter->first];
				Expert* e=(*clusterSet)[sIter->second];
				e->generateSample(r,sampleValues);
				ofstream* oFile=filePtrSet[sIter->first];
				(*oFile) << aName;
				for(int j=0;j<sampleValues.size();j++)
				{
					(*oFile) <<"\t" << sampleValues[j];
				}
				(*oFile) << endl;
				sampleValues.clear();
			}
		}
		int lcacnt=0;
		map<string,int> localKey;
		const char* dupAnc=gammaMgr->getDupAncestor(ogid);
		bool wgd=false;
		if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
		{
			wgd=true;
		}
		for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
		{
			map<string,string>* gname=geneNamesSet[hIter->first];
			map<string,int>* gset=hIter->second;
			int ancClusterAssign=-3;
			if(gset->find(lcaName)!=gset->end())
			{
				ancClusterAssign=(*gset)[lcaName];
				lcacnt++;
			}
			//Now use the scername to show this og
			string dispName;
			if(gset->find(srcSpecies)==gset->end())
			{
				char tempName[1024];	
				sprintf(tempName,"OG%d_%d",ogid,hIter->first);
				dispName.append(tempName);
				if(wgd)
				{
					dispName.append("*");
				}
			}	
			else
			{
				string geneName((*gname)[srcSpecies]);
				shownGenes[geneName]=0;
				dispName.append(gnm.getCommonName((*gname)[srcSpecies].c_str()));
				if(localKey.find(dispName)!=localKey.end())
				{
					dispName.append("_2");
				}
				localKey[dispName]=0;
				if(wgd)
				{
					dispName.append("*");
				}
			}
			allAssignFile <<dispName;
			string members;
			for(int s=0;s<speciesList.size();s++)
			{
				int geneclusterID=-2;
				if(gset->find(speciesList[s])==gset->end())
				{
					geneclusterID=-2;
				}
				else 
				{	
					geneclusterID=(*gset)[speciesList[s]];
				}
				if(members.length()>0)
				{
					members.append(";");
				}	
				if(gname->find(speciesList[s])==gname->end())
				{
					members.append("NULL");
				}
				else
				{
					members.append((*gname)[speciesList[s]]);
				}		
				allAssignFile<< "\t"<< geneclusterID;
			}
			allAssignFile << endl;	
		}
	}
	for(map<string,ofstream*>::iterator fIter=clusterFilePtrSet.begin();fIter!=clusterFilePtrSet.end();fIter++)
	{
		fIter->second->close();
		if(filePtrSet.find(fIter->first)==filePtrSet.end())
		{
			continue;
		}
		ofstream* file=filePtrSet[fIter->first];
		file->close();
	}
	
	return 0;
}



/*
* Prints out all the files to a results directory.
* Useful for printing out the EMINT-initialization results before running DRMN.
*/
int
SpeciesClusterManager::printCurrentData(const char* allOut)
{
	cerr << "scmgr::printCurrentData THIS FUNCTION ISN'T IN USE " << endl;
	return 0;

	// Dump EMINT stuff
	// emint output
	/*char outputDir[1024];
	sprintf(outputDir,"%s/emint",allOut);

	dumpAllInferredClusterAssignments(outputDir);
	showClusters_Extant(outputDir);
	showClusters_Ancestral(outputDir);
	showMeans(outputDir);

	SpeciesDistManager* sdMgr=gammaMgr->getSpeciesDistManager();

	vector<string> speciesList;
	sdMgr->getSpeciesListPrefix(speciesList); // this one is in LCA order
	//mor->getSpeciesOrder(speciesList);

	dumpAllInferredClusters_LCA(outputDir,speciesList,sdMgr->getRoot()->name);
	dumpAllInferredClusterGammas(outputDir,speciesList);
	sdMgr->showInferredConditionals(outputDir);

	// Dump feature data for DRMN and regulator/target lists 
	char featDirName[1024];	
	sprintf(featDirName,"mkdir -p %s/features",allOut);
	system(featDirName);
	sprintf(featDirName,"%s/features",allOut);
	
	dumpFeatures(featDirName, speciesList);

	// make OGIds files
	cout << "Printing OGIDs file, regulators, target files " << endl;
	makeOGIDsForRegnet(featDirName, speciesList);  

	// Print species trees per cluster
	// for each cluster, print the maintenance probs
	// and print OGIDS file for MRTLE
	for (int k=0; k<maxClusterCnt; k++)
	{
		cout << "Gain/loss probabilities for cluster " << k << endl;
		sdMgr->showTreeForCluster(k);

		cout << "Printing tree" << endl;
		char myname[1024];
		sprintf(myname, "%s/state%d_tree.txt", featDirName, k);
		sdMgr->printTreeForCluster(k,myname);
	}
	return 0;*/
}

int
SpeciesClusterManager::estimateCov(Expert* e, Matrix* X, Matrix* Y, map<int,INTDBLMAP*>& gCovar)
{
	map<string,int>& regSet=e->getCurrentRegSet(); // current regulators
	for(map<string,int>::iterator vIter=regSet.begin();vIter!=regSet.end();vIter++)
	{
		int vXID = allRegulatorIndex[vIter->first];
		int vID  = varNameIDMap[vIter->first];
		gsl_vector_view Dv = X->getRowView(vXID);

		INTDBLMAP* vcov=NULL;
		if(gCovar.find(vID)==gCovar.end())
		{
			vcov=new INTDBLMAP;
			gCovar[vID]=vcov;
		}
		else
		{
			vcov=gCovar[vID];
		}

		for(map<string,int>::iterator uIter=vIter;uIter!=regSet.end();uIter++)
		{
			int uXID = allRegulatorIndex[uIter->first];
			int uID  = varNameIDMap[uIter->first];
			gsl_vector_view Du = X->getRowView(uXID);
			double uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
			if (std::isnan(uvcov))
			{
				cout << "vID:" << vID << ",uID:" << uID << " is nan." << endl;
				uvcov = 0;
			}

			INTDBLMAP* ucov=NULL;
			if(gCovar.find(uID)==gCovar.end())
			{
				ucov=new INTDBLMAP;
				gCovar[uID]=ucov;
			}
			else
			{
				ucov=gCovar[uID];
			}
			if (uID == vID)
			{
				uvcov += 0.001;
			}
			(*vcov)[uID]=uvcov;
			(*ucov)[vID]=uvcov;
		}
		int uID=varNameIDMap["Expression"];
		gsl_vector_view Du = Y->getRowView(0);
		double uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
		if (std::isnan(uvcov))
		{
			cout << "vID:" << vID << ",uID:" << uID << " is nan." << endl;
			uvcov = 0;
		}

		INTDBLMAP* ucov=NULL;
		if(gCovar.find(uID)==gCovar.end())
		{
			ucov=new INTDBLMAP;
			gCovar[uID]=ucov;
		}
		else
		{
			ucov=gCovar[uID];
		}
		(*vcov)[uID]=uvcov;
		(*ucov)[vID]=uvcov;
	}
	int uID=varNameIDMap["Expression"];
	gsl_vector_view Du = Y->getRowView(0);
	double uvar  = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Du.vector.data, Du.vector.stride, Du.vector.size); 
	INTDBLMAP* ucov=NULL;
	if(gCovar.find(uID)==gCovar.end())
	{
		ucov=new INTDBLMAP;
		gCovar[uID]=ucov;
	}
	else
	{
		ucov=gCovar[uID];
	}
	(*ucov)[uID]=uvar+0.001;
	return 0;
}


int
SpeciesClusterManager::estimateCov_All(Expert* e, Matrix* X, Matrix* Y, map<int,INTDBLMAP*>& gCovar)
{
	for(map<string,int>::iterator vIter=allRegulatorIndex.begin();vIter!=allRegulatorIndex.end();vIter++)
	{
		int vXID = allRegulatorIndex[vIter->first];
		int vID  = varNameIDMap[vIter->first];
		gsl_vector_view Dv = X->getRowView(vXID);

		INTDBLMAP* vcov=NULL;
		if(gCovar.find(vID)==gCovar.end())
		{
			vcov=new INTDBLMAP;
			gCovar[vID]=vcov;
		}
		else
		{
			vcov=gCovar[vID];
		}

		//for(map<string,int>::iterator uIter=vIter;uIter!=regSet.end();uIter++)
		for(map<string,int>::iterator uIter=vIter;uIter!=allRegulatorIndex.end();uIter++)
		{
			int uXID = allRegulatorIndex[uIter->first];
			int uID  = varNameIDMap[uIter->first];
			gsl_vector_view Du = X->getRowView(uXID);
			double uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
			if (std::isnan(uvcov))
			{
				cout << "vID:" << vID << ",uID:" << uID << " is nan." << endl;
				uvcov = 0;
			}

			INTDBLMAP* ucov=NULL;
			if(gCovar.find(uID)==gCovar.end())
			{
				ucov=new INTDBLMAP;
				gCovar[uID]=ucov;
			}
			else
			{
				ucov=gCovar[uID];
			}
			if (uID == vID)
			{
				uvcov += 0.001;
			}
			(*vcov)[uID]=uvcov;
			(*ucov)[vID]=uvcov;
		}
		int uID=varNameIDMap["Expression"];
		gsl_vector_view Du = Y->getRowView(0);
		double uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
		if (std::isnan(uvcov))
		{
			cout << "vID:" << vID << ",uID:" << uID << " is nan." << endl;
			uvcov = 0;
		}

		INTDBLMAP* ucov=NULL;
		if(gCovar.find(uID)==gCovar.end())
		{
			ucov=new INTDBLMAP;
			gCovar[uID]=ucov;
		}
		else
		{
			ucov=gCovar[uID];
		}
		(*vcov)[uID]=uvcov;
		(*ucov)[vID]=uvcov;
	}
	int uID=varNameIDMap["Expression"];
	gsl_vector_view Du = Y->getRowView(0);
	double uvar  = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Du.vector.data, Du.vector.stride, Du.vector.size); 
	INTDBLMAP* ucov=NULL;
	if(gCovar.find(uID)==gCovar.end())
	{
		ucov=new INTDBLMAP;
		gCovar[uID]=ucov;
	}
	else
	{
		ucov=gCovar[uID];
	}
	(*ucov)[uID]=uvar+0.001;
	return 0;
}


int
SpeciesClusterManager::addCov(Expert* e, Matrix* X, Matrix* Y, map<int,INTDBLMAP*>& gCovar, string regName)
{
	map<string,int>& regSet=e->getCurrentRegSet(); // current regulators
	int vXID = allRegulatorIndex[regName];
	int vID  = varNameIDMap[regName];
	gsl_vector_view Dv = X->getRowView(vXID);

	INTDBLMAP* vcov=NULL;
	if(gCovar.find(vID)==gCovar.end())
	{
		vcov=new INTDBLMAP;
		gCovar[vID]=vcov;
	}
	else
	{
		vcov=gCovar[vID];
	}

	for(map<string,int>::iterator uIter=regSet.begin();uIter!=regSet.end();uIter++)
	{
		int uXID = allRegulatorIndex[uIter->first];
		int uID  = varNameIDMap[uIter->first];
		gsl_vector_view Du = X->getRowView(uXID);
double uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
		if (std::isnan(uvcov))
		{
			cout << "vID:" << vID << ",uID:" << uID << " is nan." << endl;
			uvcov = 0;
		}

		INTDBLMAP* ucov=NULL;
		if(gCovar.find(uID)==gCovar.end())
		{
			ucov=new INTDBLMAP;
			gCovar[uID]=ucov;
		}
		else
		{
			ucov=gCovar[uID];
		}
		if (uID == vID)
		{
			uvcov += 0.001;
		}
		if(vcov->find(uID)!=vcov->end())
		{
			cout <<"vcov already has an entry " << (*vcov)[uID] << " and the new val is " << uvcov << endl;
		}
		else
		{
			(*vcov)[uID]=uvcov;
			(*ucov)[vID]=uvcov;
		}
	}
	int uID=varNameIDMap["Expression"];
	gsl_vector_view Du = Y->getRowView(0);
double uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
	if (std::isnan(uvcov))
	{
		cout << "vID:" << vID << ",uID:" << uID << " is nan." << endl;
		uvcov = 0;
	}

	INTDBLMAP* ucov=NULL;
	if(gCovar.find(uID)==gCovar.end())
	{
		ucov=new INTDBLMAP;
		gCovar[uID]=ucov;
	}
	else
	{
		ucov=gCovar[uID];
	}
	if(vcov->find(uID)!=vcov->end())
	{
		cout <<"vcov already has expression entry " << (*vcov)[uID] << " and the new val is " << uvcov << endl;
	}
	else
	{
		(*vcov)[uID]=uvcov;
		(*ucov)[vID]=uvcov;
	}
	return 0;
}

int
SpeciesClusterManager::setMode(LearnMode lm, double r1, double r2, double r3)
{
	learnMode = lm;
	rho1 = r1;
	rho2 = r2;
	rho3 = r3;
	return 0;
}

int 
SpeciesClusterManager::learnLEAST(vector<Task_T*>* allt,vector<Matrix*>& allW)
{
	map<int,vector<int>> tree;
	for (int i=0;i<allt->size()-1;i++)
	{
		vector<int> e;
		e.push_back(i+1);
		tree[i] = e;
	}
	GenericLearner* ll;
	if (learnMode == LEASTDIRTY)
	{
		cout << "LEASTDIRTY!" << endl;
	//ll = new LeastDirty(newlambda,200);
		vector<Matrix*> allP;
		vector<Matrix*> allQ;
		//ll = new LeastDirty(newlambda,newlambda);
		ll = new LeastDirty(rho1,rho2);
		ll->doAGM(allt,allW,allP,allQ);
		delete ll;
		clearAllW(allP);
		clearAllW(allQ);
	}
	else if (learnMode == LEASTL21)
	{
		cout << "LEASTL21!" << endl;
		ll = new LeastL21(rho1,rho2);
		ll->doAGM(allt,allW);
		delete ll;
	}
	else if (learnMode == LEASTFUSED)
	{
		cout << "LEASTFUSED!" << endl;
		ll = new LeastCFGLasso(rho1,rho2,rho3);
		ll->doAGM(allt,allW,tree);
		delete ll;
	}
	else
	{
		cerr << "UNKNOWN LEARNING METHOD, shouldnt happen" << endl;
		exit(0);
	}
	return 0;
}

double
SpeciesClusterManager::makeOnePotential(Expert* e,vector<string>& regNames,string& cellType,bool& regStatus,DRMNPotential** potPtr, Matrix* X, Matrix* Y, int moduleID, int iter) 
{
	//SpeciesFeatureManager* sfMgr=speciesFeatSet[cellType];
	EvidenceManager* evMgr=evMgrSet[cellType];
	map<string,int>& regSet=e->getCurrentRegSet();
	
	map<string,int>& geneSet=e->getGeneSet();

	// if geneset is 0, issue? this means the expert has no genes! :(
	// I think this might be happening earlier
	if (geneSet.size()==0)
	{
		cout << "geneset size 0 for expert!" << endl;
		return 0;
	}

	//We will add the regulator to this module's regulator set and remove it when we are done
	regSet.clear();
	for(int i=0;i<regNames.size();i++)
	{
		string regName = regNames[i];
		regSet[regName]=0;
	}
	//We are going to use our EvidenceManager associated with this cell type to compute everything. No need to make copies of the data. Ew.
	int colID=0;
	//First we will compute the mean, and then the covariance (this is the equivalent of calling populatePotential
	//To compute the mean  we will first consider the mRNA levels 
	//Then we can calculate the coefficients using something initMBCovMean

	cout << "estimateMean:" << endl;
	INTDBLMAP mean; 
	for(map<string,int>::iterator geneIter=geneSet.begin();geneIter!=geneSet.end();geneIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt((string&)geneIter->first);
		// first get expression
		int varID=varNameIDMap["Expression"];
		//Evidence* evid=(*evidMap)[varID];
		//double val=evid->getEvidVal();
		double val=(*evidMap)[varID];
		if(mean.find(varID)==mean.end())
		{
			mean[varID]=val;
		}
		else	
		{
			mean[varID]=mean[varID]+val;
		}
		//Now get the existing and new regulators
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			int varID=varNameIDMap[rIter->first];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			if(mean.find(varID)==mean.end())
			{
				mean[varID]=val;
			}
			else	
			{
				mean[varID]=mean[varID]+val;
			}
		}
	}
	
	//Now we just normalize to get the mean
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double val=vIter->second/geneSet.size();
		vIter->second=val;
	}
	//Mean done. We now need the covariance. 
	//We can just go over each evidence and populate the covariance entry which we will then normalize

	// NEW STUFF: We'll copy in the base unnormalized covariance matrix.
	map<int,INTDBLMAP*> gCovar;
	//computeBaseUnnormCovar(e, cellType, gCovar);
	cout << "estimateCov:" << endl;
	estimateCov(e, X, Y, gCovar);
	cout << "gCovar:" << gCovar.size() << endl;

	//{
	//	char tempoutname[1024];
	//	sprintf(tempoutname,"%s_%d_%d.txt",cellType.c_str(),iter,moduleID);
	//	ofstream tempOF(tempoutname);
	//	for (map<int,INTDBLMAP*>::iterator iitr=gCovar.begin();iitr!=gCovar.end();iitr++)
	//	{
	//		int curID = iitr->first;
	//		string curName = varIDNameMap[curID];
	//		tempOF << curName;
	//		INTDBLMAP* curcov = iitr->second;
	//		for (map<int,double>::iterator jitr=curcov->begin();jitr!=curcov->end();jitr++)
	//		{
	//			tempOF << "\t" << jitr->second;
	//		}
	//		tempOF << endl;
	//	}
	//	tempOF.close();
	//}
	

	DRMNPotential* pot=new DRMNPotential;
	DRMNPotential* parentPot=new DRMNPotential;
	int expVarID=-1;
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		string& vName=varIDNameMap[vIter->first];
		if(strcmp(vName.c_str(),"Expression")==0)
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::FACTOR);
			expVarID=vIter->first;
		}
		else
		{
			pot->setAssocVariable(vIter->first,DRMNPotential::MARKOV_BNKT);
			parentPot->setAssocVariable(vIter->first,DRMNPotential::MARKOV_BNKT);
		}
	}
	pot->potZeroInit();
	parentPot->potZeroInit();
	//populate potential using mean and covariance
	for(INTDBLMAP_ITER vIter=mean.begin();vIter!=mean.end();vIter++)
	{
		double m=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		m=mean[vIter->first];
		
		covar=gCovar[vIter->first];
		pot->updateMean(vIter->first,m);
		if(vIter->first!=expVarID)
		{
			parentPot->updateMean(vIter->first,m);
		}
		for(INTDBLMAP_ITER uIter=vIter;uIter!=mean.end();uIter++)
		{
			double cval=(*covar)[uIter->first];
			pot->updateCovariance(vIter->first,uIter->first,cval);
			pot->updateCovariance(uIter->first,vIter->first,cval);
			if(uIter->first!=expVarID && vIter->first!=expVarID)
			{
				parentPot->updateCovariance(vIter->first,uIter->first,cval);
				parentPot->updateCovariance(uIter->first,vIter->first,cval);
			}
		}
	}
	//pot->makeValidJPD(ludecomp,perm);
	int status=pot->initMBCovMean();
	

	if(status!=0)
	{
		cout <<"Determinant too small, not going ahead with making potential in cell " << cellType  << endl;
		delete pot;
		return 0;
	}

	//Now we have estimated everything! Now we can use the same logic as before to get the score
	//double currScore=getScoreForPot(geneSet,pot,cellType);
	double currScore=getScoreForPot_Tracetrick(geneSet,pot,parentPot,cellType);
	double oldScore=e->getLLScore();
	double impr=currScore-oldScore;
	mean.clear();
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		idIter->second->clear();
		delete idIter->second;
	}
	gCovar.clear();
	//potPtr=&pot;
	*potPtr=pot;
	//potPtr=pot;

	return impr;
}
int
SpeciesClusterManager::countAllW(vector<Matrix*> allW)
{
	int cnt = 0;
	for(int i=0;i<allW.size();i++)
	{
		Matrix* W = allW[i];
		for (int j=0;j<W->getRowCnt();j++)
		{
			double v = W->getValue(j,0);
			//cout << "# " << v << endl;
			//printf("# %.20f\n", v);
			//if ( fabs(v)>0.001 )
			if ( v != 0 )
			{
				cnt ++;
			}
		}
	}
	cnt = cnt/allW.size();
	return cnt;
}

int
SpeciesClusterManager::clearAllW(vector<Matrix*>& allW)
{
	for(int i=0;i<allW.size();i++)
	{
		Matrix* W = allW[i];
		delete W;
	}
	allW.clear();
	return 0;
}

int
SpeciesClusterManager::estimateRegProgs_PerModule_LASSO(int moduleID,int iter)
{
	cout << "Updating regulatory program for " << moduleID << " ";
	time_t result = time(0);
	cout << std::asctime(localtime(&result)) << endl;
	if (learnMode == LEASTDIRTY)
	{
		cout << "Doing LeastDIRTY!" << endl;
	}
	else if (learnMode == LEASTL21)
	{
		cout << "Doing LeastL21!" << endl;
	}
	else if (learnMode == LEASTFUSED)
	{
		cout << "Doing LeastFUSED!" << endl;
	}
	
	vector<Task_T*>* allt = new vector<Task_T*>;
	map<string,Matrix*> cell2exp;
	map<string,Matrix*> cell2feat;

	vector<string> speciesList;
	mor->getSpeciesOrder(speciesList);
	for(int sindex=0;sindex<speciesList.size();sindex++)
	{
		string s = speciesList[sindex];
		map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.find(s);
		//string s=sIter->first;
		cerr << "in estimateRegProgs_PerModule_LASSO, index: " << sindex << ", species: " << s << endl;
		//sindex++;

		CLUSTERSET* expertSet=sIter->second;
		Expert* expert=(*expertSet)[moduleID];
		map<string,int>& geneSet=expert->getGeneSet();
		//SpeciesFeatureManager* sfMgr=speciesFeatSet[s];

		map<string,int> temp_gene2id;
		temp_gene2id.clear();
		int tgi=0;
		for (map<string,int>::iterator gitr = geneSet.begin(); gitr!=geneSet.end(); gitr++)
		{
			string gname = gitr->first;
			temp_gene2id[gname] = tgi;
			tgi++;
		}

		int regCnt = allRegulatorIndex.size();
		int geneCnt = geneSet.size();
		cout << s << " has " << geneCnt << " genes." << endl;
		EvidenceManager* evMgr=evMgrSet[s];
		Matrix* Y = new Matrix(geneCnt,1);
		map<string,int>* YNames = new map<string,int>;
		for (map<string,int>::iterator gitr = geneSet.begin(); gitr!=geneSet.end(); gitr++)
		{
			string gname = gitr->first;
			int gi=temp_gene2id[gname];
			EMAP* evidMap=evMgr->getEvidenceAt(gname);
			int varID=varNameIDMap["Expression"];
			//Evidence* evid=(*evidMap)[varID];
			//double val=evid->getEvidVal();
			double val=(*evidMap)[varID];
			Y->setValue(val,gi,0);
			(*YNames)[gname]=gi;
		}
		Matrix* YT = Y->transMatrix();
		//YT->rowStandardize();
		cell2exp[s] = YT;
		Matrix* X = new Matrix(regCnt,geneCnt);
		for (map<string,int>::iterator gitr = geneSet.begin(); gitr!=geneSet.end(); gitr++)
		{
			string gname = gitr->first;
			int gi = temp_gene2id[gitr->first];
			EMAP* evidMap = evMgr->getEvidenceAt(gname);
			for (map<string,int>::iterator ritr=allRegulatorIndex.begin(); ritr!=allRegulatorIndex.end();ritr++)
			{
				int ri=ritr->second;
				int varID=varNameIDMap[ritr->first];
				//Evidence* evid = (*evidMap)[varID];
				//double val = evid->getEvidVal();
				double val = (*evidMap)[varID];
				X->setValue(val,ri,gi);
			}
		}
		Matrix* XX = X->copyMe();
		//XX->rowStandardize();
		cell2feat[s] = XX;
		X->rowStandardize();
		Y->colStandardize();
		Task_T* t = new Task_T;
		t->YNames = YNames;
		t->X = X;
		t->Y = Y;
		t->name = s;
		allt->push_back(t);
	}
	vector<Matrix*> allW;
	learnLEAST(allt,allW);
	int cnt = countAllW(allW);
	cout << "We got (average): " << cnt << " regulators." << endl;
	/*
	double oldlambda=0;
	double newlambda=100;
	while (true)
	{
		cout << "Using lambda: " << newlambda << endl;
		learnLEAST(allt,allW,newlambda);
		int cnt = countAllW(allW);
		cout << "We got (average): " << cnt << " regulators." << endl;
		//if (cnt==0)
		if (cnt<minRegSize)
		{
			newlambda = (newlambda+oldlambda)/2;
			clearAllW(allW);
		}
		else if (cnt > maxRegSize)
		{
			if (newlambda>100000)
			{
				cout << "We failed to find the right number of regulators, stopping at lambda=" << newlambda << endl;
				break;
			}
			oldlambda = newlambda;
			newlambda = 2*oldlambda;
			clearAllW(allW);
		}
		else
		{
			break;
		}
	}
	*/


	result = time(0);
	cout << "Finished LeastL21" << " ";
	cout << std::asctime(localtime(&result)) << endl;
	
	result = time(0);
	cout << "Start adding" << " ";
	cout << std::asctime(localtime(&result)) << endl;

	for(int i=0;i<allW.size();i++)
	{
		Task_T* t = allt->at(i);
		string  s = t->name;
		Matrix* W = allW[i];
		Matrix* X = cell2feat[s];
		int gcnt = X->getColCnt();
		Matrix* Y = cell2exp[s];

		CLUSTERSET* expertSet=speciesExpertSet[s];
		Expert* expert=(*expertSet)[moduleID];
		//SpeciesFeatureManager* sfMgr=speciesFeatSet[s];
		vector<string> regNames;
		vector<int>    regIds;
		int v0cnt = 0;
		int skipCnt = 0;
		for (map<string,int>::iterator ritr=allRegulatorIndex.begin(); ritr!=allRegulatorIndex.end();ritr++)
		{
			int ri=ritr->second;
			string rname = ritr->first;
			double v = W->getValue(ri,0);
			//Add this one
			//if ( fabs(v)>0.001 )
			if ( v != 0 )
			{
				v0cnt++;
				//check repeat
				int skip = 0;
				int skipr = 0;
				for (int j=0;j<regIds.size();j++)
				{
					int rj = regIds[j];
					int eq = 1;
					for (int k=0;k<gcnt;k++)
					{
						double vi = X->getValue(ri,k);
						double vj = X->getValue(rj,k);
						if (vi != vj)
						{
							eq = 0;
							break;
						}
					}
					if (eq == 1)
					{
						skip = 1;
						skipr = j;
						break;
					}
				}
				if (skip != 1)
				{
					regNames.push_back(rname);
					regIds.push_back(ri);
					//cout << ".";
				}
				else
				{
					skipCnt++;
					cout << "SKIP! in module " << moduleID << " Reg " << rname << " is the same as " << regNames[skipr] << endl;
				}
			}
		}
		cout << "I had " << v0cnt << ", skipped " << skipCnt << endl;
		cout << "I have " << regNames.size() << " regulators. Let's make a POTENTIAL!" << endl;
		if (regNames.size() == 0)
		{
			continue;
		}
		bool regulatorStatus=false;
		DRMNPotential* aPot=NULL;
		double scoreImprovement_PerCelltype=makeOnePotential(expert,regNames,s,regulatorStatus,&aPot,X,Y,moduleID,iter);

		if (aPot == NULL)
		{
			cout << "Something went wrong , skipping this one" << endl;
			//exit(0);
		}
		else
		{
			double currScore=expert->getLLScore();
			expert->setLLScore(currScore+scoreImprovement_PerCelltype);
			DRMNPotential* tPot=expert->getDRMNPotential();
			if(tPot!=NULL)
			{
				delete tPot;
			}
			expert->setDRMNPotential(aPot);
		}
		cout << endl;
	}
	
	cout << "Done adding" << " ";
	cout << std::asctime(localtime(&result)) << endl;

	cout <<"Done learning regulators for moduleID " << moduleID <<endl;
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* experts=speciesExpertSet[sIter->first];
		Expert* e=(*experts)[moduleID];
		map<string,int>& regSet=e->getCurrentRegSet();
		if (regSet.size()==0)
		{
			cout << "No regulators for Module" << moduleID << endl;
		}
		for(map<string,int>::iterator rIter=regSet.begin();rIter!=regSet.end();rIter++)
		{
			cout << sIter->first<< " Module"<<moduleID<<"-" << rIter->first <<endl;
		}
	}

	for (int i=0;i<allt->size();i++)
	{
		Task_T* t = allt->at(i);
		string s = t->name;
		Matrix* XX = cell2feat[s];
		delete XX;
		Matrix* Y  = cell2exp[s];
		delete Y;
		delete t;
		Matrix* W = allW[i];
		delete W;
	}
	allt->clear();
	delete allt;
	allW.clear();
	cell2feat.clear();
	cell2exp.clear();

	return 0;
}

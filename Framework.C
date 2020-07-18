/*
CMINT: An algorithm to cluster functional genomics data for multiple cell types
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
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "Matrix.H"
#include "SpeciesDistManager.H"
#include "GeneTree.H"
#include "GeneTreeManager.H"
#include "Expert.H"
#include "Gamma.H"
#include "GammaManager.H"
#include "GeneNameMapper.H"
#include "SpeciesClusterManager.H"
#include "Error.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "Framework.H"


Framework::Framework()
{
}

Framework::~Framework()
{
}


int 
Framework::readSpeciesData(const char* aFName, const char* rand)
{
	if(strcmp(rand,"none")==0)
	{
		scMgr.setRandom(false);
	}
	else if(strcmp(rand,"yes")==0)
	{
		scMgr.setRandom(true);
	}
	else if(isdigit(rand[0]))
	{
		scMgr.setRandom(true);
		scMgr.setRandSeed(atoi(rand));
	}
	scMgr.setMaxClusterCnt(maxClusterCnt);

	// get species list from sdMgr -- this way we don't read data for other species
	// the RNG for scMgr gets set up in readSpeciesData
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);
	int dataOK=scMgr.readSpeciesData(aFName, speciesList);
	if (dataOK != 0)
	{
		return 1;
	}
	
	// this is an RNG for Framework that is different from scMgr's
	randnum=gsl_rng_alloc(gsl_rng_default);
	return 0;
}

int 
Framework::readSpeciesTree(int clusterCnt, const char* aFName)
{
	maxClusterCnt=clusterCnt;
	sdMgr.setMaxClusters(clusterCnt);
	sdMgr.readSpeciesTree(aFName);
	sdMgr.assignLevel();
	gammaMgr.setSpeciesDistManager(&sdMgr);
	return 0;
}

int 
Framework::readOrthology(const char* specOrder, const char* orthomapfile)
{
	int orderOK = mor.readSpeciesMapping(specOrder);
	if (orderOK != 0)
	{
		cerr << "Problem reading cell order file: " << specOrder << endl;
		return 1;
	}
	mor.readFile(orthomapfile);
	scMgr.setOrthogroupReader(&mor);
	gammaMgr.setOrthogroupReader(&mor);
	gammaMgr.setMaxClusterCnt(maxClusterCnt);
	scMgr.setGammaManager(&gammaMgr);
	return 0;
}

/*
* Reads regulator OG list.
*/
int Framework::readRegulatorOGIds(const char* regfile)
{
	int success=scMgr.setRestrictedList(regfile);
	if (success != 0)
	{
		cerr << "Problem reading regulator list file: " << regfile << endl;
	}
	return success;
}

int
Framework::setSrcSpecies(const char* specName)
{
	strcpy(srcSpecies,specName);
	scMgr.setSrcSpecies(specName);
	return 0;
}

/**
* This is the main function.
*
*/
int
Framework::startClustering(const char* aDir)
{
	strcpy(outputDir,aDir);
	scMgr.initExperts();
	cout <<"Total updated parent nodes "<< gammaMgr.getTotalUpdatedParentCnt() << endl; 
	gammaMgr.showTotalUpdatedParents();
	initClusterTransitionProb(); // initialize from input

	
	// emint output
	char dirName[1024];
	sprintf(dirName,"mkdir -p %s/emint",outputDir);
	int errcode = system(dirName);
	if (errcode != 0)
	{
		cerr << "Could not make directory " << dirName << endl;
		return errcode;
	}

	sprintf(dirName,"%s/emint",outputDir);


	cout << "Running EMINT; putting results in " << dirName << endl;

	//SR: We will rename this to initializeExperts
	scMgr.estimateExpertParameters(dirName); // Run EMINT until convergence.
	// print the emint stuff
	printResults(dirName);
	// print current data and DRMN data to subdirs emint/ and features/
	//scMgr.printCurrentData(outputDir); // not used

	// for each cluster, print the maintenance probs
	/*for (int k=0; k<maxClusterCnt; k++)
	{
		cout << "Maintenance probs for cluster " << k << endl;
		sdMgr.showTreeForCluster(k);
		
	}*/
	

	// species list
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);


	//return 0;

	// Get maximal cluster assignments and update transition probs
	// die if cannot create dir
	sprintf(dirName,"mkdir -p %s/drmn",outputDir);
	errcode = system(dirName);
	if (errcode != 0)
	{
		cerr << "Could not create output directory " << outputDir << endl;
		return errcode; 
	}

	sprintf(dirName,"%s/drmn",outputDir);
	cout << "Will put drmn results in " << dirName << endl;
	char drmnOutputDir[1024];
	sprintf(drmnOutputDir,"%s/drmn",outputDir);
	scMgr.setMaxAssignments(); // Redundant with the printing function, but OK -- need to make sure we do this.
	
	// check for errors
	int ran = scMgr.estimateDRMN(outputDir);
	if (ran != 0)
	{
		cerr << "Could not estimate DRMN." << endl;
		return 1;
	}

	double newScore=scMgr.getScore();
	scMgr.dumpAllInferredClusterAssignments(outputDir);
	
	double newScore_PP=scMgr.getScore();
	//cout <<"Score before PP " << newScore << "\t" << " Score after PP " << newScore_PP << endl;
	
	scMgr.showClusters_Extant(outputDir);
	scMgr.showClusters_Ancestral(outputDir);
	scMgr.showMeans(outputDir);
	//This is only for visualization purposes
	/*vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);
	for(int i=0;i<speciesList.size();i++)
	{	
		cout << speciesList[i] << endl;
	}*/

	//scMgr.dumpAllInferredClusters_SrcwiseGrouped(outputDir,speciesList);
	scMgr.dumpAllInferredClusters_LCA(outputDir,speciesList,sdMgr.getRoot()->name);
	scMgr.dumpAllInferredClusterGammas(outputDir,speciesList);
	
	int printed = scMgr.showDRMNResults(drmnOutputDir);
	if (printed !=0)
	{
		cerr << "Could not print final results." << endl;
		return printed;
	}
	

	//sdMgr.showInferredConditionals_ML(outputDir);
	return 0;
}

//This function will implement the logic of startClustering for different folds of the data
// note that folds are determined by OGIDs, which may be a larger set than the initial cluster assignments
// int myFold: which fold (out of total folds) to run this time. If < 0, run all folds.
// int foldSeed: seed to RNG for ordering examples into folds. if < 0, use PID. 
// If foldSeed < 0 and myFold >=0, then we are doing sampling with replacement.

int
Framework::startClustering_CV(const char* aFName, const char* rand, const char* aDir,int folds, int myFold, int foldSeed)
{

	if(strcmp(rand,"none")==0)
	{
		scMgr.setRandom(false); // no randomization of clusterIDs
	}
	else if(strcmp(rand,"yes")==0)
	{ 
		scMgr.setRandom(true); // set cluster IDs randomly
	}
	else if(isdigit(rand[0]))
	{
		scMgr.setRandom(true);  
		scMgr.setRandSeed(atoi(rand)); // set cluster ID randomizer seed
	}
	strcpy(outputDir,aDir);

	// This is the RNG for CV. 
	randnum=gsl_rng_alloc(gsl_rng_default);  
	if (foldSeed >= 0) // user provided seed so that we can make folds without replacement
	{
		gsl_rng_set(randnum, foldSeed);
		cout << "Set RNG seed to: " << foldSeed << endl;
	}

	scMgr.setMaxClusterCnt(maxClusterCnt);

	EvidenceManager dummyEvManager;
	map<int,MappedOrthogroup*>& ogSet=mor.getMappedOrthogroups();
	vector<int> ogSet_IDs;
	for(map<int,MappedOrthogroup*>::iterator oIter=ogSet.begin();oIter!=ogSet.end();oIter++)
	{
		ogSet_IDs.push_back(oIter->first);	
	}
	int* randInds=new int[ogSet.size()];
	// this is where we populate randInds with the order of examples
	dummyEvManager.populateRandIntegers(randnum,randInds,ogSet.size());

	int testBegin=0;
	int foldSize=(int)(ogSet.size()/folds); 
	for(int f=0;f<folds;f++)
	{
		int testEnd=testBegin+foldSize;
		map<int,int> trainIDs;
		map<int,int> testIDs;
		for(int i=0;i<ogSet.size();i++)
		{
			int randid=randInds[i];
			int ogid=ogSet_IDs[randid];
			if(i>=testBegin && i<testEnd)
			{
				testIDs[ogid]=0;
			}
			else
			{
				trainIDs[ogid]=0;
			}
		}
		testBegin=testEnd;

		// if we have specified one fold to run and it's not this, skip.
		if (myFold >= 0 && myFold != f)
		{
			cout << "SKIPPING fold " << f << endl;
			continue;	
		}
		cout << "NOW RUNNING fold " << f << endl;

		scMgr.setTrainOGIDs(trainIDs);
		
		// species list
		vector<string> mySpecies;
		sdMgr.getSpeciesListPrefix(mySpecies);
		int dataOK=scMgr.readSpeciesData(aFName, mySpecies);

		if (dataOK != 0)
		{
			cerr << "Problem reading species data in fold " << f << endl;
			return 1;
		}
		scMgr.initExperts();
		cout <<"Total updated parent nodes "<< gammaMgr.getTotalUpdatedParentCnt() << endl; 
		gammaMgr.showTotalUpdatedParents();
		initClusterTransitionProb(); // initialize from input
	
		// emint output
		char dirName[1024];
		sprintf(dirName,"mkdir -p %s/fold%d/emint",outputDir,f);		
		int errcode = system(dirName);
		if (errcode != 0)
		{
			cerr << "Could not create output directory " << outputDir << endl;
			return 1;
		}

		sprintf(dirName,"%s/fold%d/emint",outputDir,f);
		cout << "Running EMINT; putting results in " << dirName << endl;
		
		//SR: We will rename this to initializeExperts
		scMgr.estimateExpertParameters(dirName); // Run EMINT until convergence.
		
		// print the emint stuff
		printResults(dirName);

		// species list
		vector<string> speciesList;
		sdMgr.getSpeciesListPrefix(speciesList);

		// Get maximal cluster assignments and update transition probs
		sprintf(dirName,"mkdir -p %s/fold%d/drmn",outputDir,f);
		errcode = system(dirName);
		if (errcode != 0)
		{
			cerr << "Could not create output directory " << outputDir << endl;
			return 1;
		}

		sprintf(dirName,"%s/fold%d/drmn",outputDir,f);
		cout << "Will put drmn results in " << dirName << endl;
		char drmnOutputDir[1024];
		sprintf(drmnOutputDir,"%s/fold%d/drmn",outputDir,f);
		scMgr.setMaxAssignments(); // Redundant with the printing function, but OK -- need to make sure we do this.
		
		int success=scMgr.estimateDRMN(drmnOutputDir);
		if (success != 0)
		{
			cerr << "Could not estimate DRMN" << endl;
			return success;
		}
		
		// collect training data likelihood
		double trainUnpen=0;
		double trainPen=0;
		scMgr.getDRMNScore_test(trainIDs, trainUnpen, trainPen);

		scMgr.dumpAllInferredClusterAssignments(drmnOutputDir);
		scMgr.showClusters_Extant(drmnOutputDir);
		scMgr.showClusters_Ancestral(drmnOutputDir);
		scMgr.showMeans(drmnOutputDir);
		scMgr.dumpAllInferredClusters_LCA(drmnOutputDir,speciesList,sdMgr.getRoot()->name);
		scMgr.dumpAllInferredClusterGammas(drmnOutputDir,speciesList);

		//Then we will generate test predictions	
		scMgr.getTestPerformance(drmnOutputDir,testIDs);

		// collect test data likelihood
		double unpenalized;
		double penalized;
		scMgr.getDRMNScore_test(testIDs, unpenalized, penalized);

		// print to console for now
		cout << "TrainDRMN UnpenalizedLL " << trainUnpen << " PenalizedLL " << trainPen << endl;
		cout << "TestDRMN UnpenalizedLL " << unpenalized << " PenalizedLL " << penalized << endl;

		// print likelihoods to file
		char llName[1024];
		sprintf(llName,"%s/fold%d/drmn/likelihood.txt",outputDir, f);
		ofstream oFile(llName);
		oFile << "Data\tUnpenalizedLL\tPenalizedLL" << endl;
		oFile << "Train_" << f << "\t" << trainUnpen << "\t" << trainPen << endl;
		oFile << "Test_" << f << "\t" << unpenalized << "\t" << penalized << endl;
		oFile.close();
		cout << "Wrote likelihood scores to " << llName << endl;

		scMgr.clear();
		trainIDs.clear();
		testIDs.clear();
	}
	return 0;
}
/*
* Prints out all the files to a results directory.
* Useful for printing out the EMINT-initialization results before running DRMN.
*/
int
Framework::printResults(const char* outputDir)
{

	scMgr.dumpAllInferredClusterAssignments(outputDir);
	scMgr.showClusters_Extant(outputDir);
	scMgr.showClusters_Ancestral(outputDir);
	scMgr.showMeans(outputDir);

	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);

	scMgr.dumpAllInferredClusters_LCA(outputDir,speciesList,sdMgr.getRoot()->name);
	scMgr.dumpAllInferredClusterGammas(outputDir,speciesList);
	sdMgr.showInferredConditionals(outputDir);
	return 0;
}


int
Framework::generateData(const char* outputDir)
{
	scMgr.initExperts();
	//Lets do one round of learning shall we.
	initClusterTransitionProb();
	scMgr.estimateExpertParameters(outputDir);
	sdMgr.showInferredConditionals(outputDir);
	scMgr.showMeans(outputDir);
	char dirName[1024];
	sprintf(dirName,"mkdir -p %s/samples",outputDir);
	system(dirName);
	sprintf(dirName,"%s/samples",outputDir);
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);

	scMgr.generateData(dirName,sdMgr.getRoot()->name,speciesList);
	return 0; // ENDS HERE


	/*char fName[1024];
	//Assume we want to generate the same number of genes as in the original data
	SpeciesDistManager::Species* root=sdMgr.getRoot();
	map<string,CLUSTERSET*>& extantSpeciesSet=scMgr.getExtantSpeciesClusters();
	map<string,ofstream*> filePtrSet;
	map<string,ofstream*> filePtrClusteredSet;
	map<string,ofstream*> clusterFilePtrSet;
	map<string,map<int,int>*> speciesClusterDist;
	map<string,map<int,map<string,int>*>*> speciesClusterMembers;
	for(map<string,CLUSTERSET*>::iterator sIter=extantSpeciesSet.begin();sIter!=extantSpeciesSet.end();sIter++)
	{
		sprintf(fName,"%s/%s_samples.txt",outputDir,sIter->first.c_str());
		ofstream* oFile=new ofstream(fName);
		filePtrSet[sIter->first]=oFile;
		sprintf(dirName,"mkdir -p %s/%s",outputDir,sIter->first.c_str());
		system(dirName);
		sprintf(fName,"%s/%s/clusterassign.txt",outputDir,sIter->first.c_str());
		ofstream* cFile=new ofstream(fName);
		clusterFilePtrSet[sIter->first]=cFile;
	}
	map<int,MappedOrthogroup*>& ogSet=mor.getMappedOrthogroups();
	vector<double> sampleValues;
	char clusterAssignmentFName[1024];
	sprintf(clusterAssignmentFName,"%s/clusterassign_multspecies.txt",outputDir);
	ofstream caFile(clusterAssignmentFName);
	string scer(srcSpecies);
	map<string,int>* scerGenes=scMgr.getGenesForSpecies(scer);
	int shown=0;
	for(map<string,int>::iterator gIter=scerGenes->begin();gIter!=scerGenes->end();gIter++)
	{
		int cId=sampleAncestralCluster(randnum,root);
		map<string,int>* clusterAssign=new map<string,int>;
		clusterAssignments[gIter->first]=clusterAssign;
		(*clusterAssign)[root->name]=cId;
        
		//sampleChildCluster(randnum,root->leftchild,cId,*clusterAssign);
		//sampleChildCluster(randnum,root->rightchild,cId,*clusterAssign);
        
        for (int i=0; i<root->children.size(); i++) {
            sampleChildCluster(randnum,root->children[i],cId,*clusterAssign);
        }
        
		caFile << gIter->first;
		for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			cout<<" " <<cIter->first <<"=" << cIter->second;
			caFile<<"\t" <<cIter->first <<"=" << cIter->second;
			map<int,int>* clusterCnt=NULL;
			if(speciesClusterDist.find(cIter->first)==speciesClusterDist.end())
			{
				clusterCnt=new map<int,int>;
				speciesClusterDist[cIter->first]=clusterCnt;
			}
			else
			{
				clusterCnt=speciesClusterDist[cIter->first];
			}
			if(clusterCnt->find(cIter->second)==clusterCnt->end())
			{
				(*clusterCnt)[cIter->second]=1;
			}
			else
			{
				(*clusterCnt)[cIter->second]=(*clusterCnt)[cIter->second]+1;
			}
		}
		caFile<< endl;
		cout << endl;
		int ogid=mor.getMappedOrthogroupID(gIter->first.c_str(),srcSpecies);
		MappedOrthogroup* mgrp=ogSet[ogid];
		//Then for all extant species draw the expression vector
		for(map<string,CLUSTERSET*>::iterator sIter=extantSpeciesSet.begin();sIter!=extantSpeciesSet.end();sIter++)
		{
			CLUSTERSET* clusterSet=sIter->second;
			int specClustId=(*clusterAssign)[sIter->first];
			Expert* e=(*clusterSet)[specClustId];
			e->generateSample(randnum,sampleValues);
			GeneMap* gMap=mgrp->getSpeciesHits(sIter->first.c_str());
			ofstream* oFile=filePtrSet[sIter->first];
			ofstream* cFile=clusterFilePtrSet[sIter->first];
			const string& geneName=gMap->getGeneSet().begin()->first;
			(*cFile) <<geneName <<"\t" << specClustId <<endl;
			(*oFile) << geneName;
			for(int j=0;j<sampleValues.size();j++)
			{
				(*oFile) <<"\t" << sampleValues[j];
			}
			(*oFile) << endl;
			sampleValues.clear();
		}
		shown++;
		//clusterAssign.clear();
	}
	for(map<string,ofstream*>::iterator fIter=filePtrSet.begin();fIter!=filePtrSet.end();fIter++)
	{
		fIter->second->close();
		ofstream* cFile=clusterFilePtrSet[fIter->first];
		cFile->close();
	}
	caFile.close();
	scMgr.showClusters_Extant(outputDir);
	
	//Cluster size dist
	for(map<string,map<int,int>*>::iterator sIter=speciesClusterDist.begin();sIter!=speciesClusterDist.end();sIter++)
	{
		cout <<sIter->first;
		map<int,int>* sizeDist=sIter->second;
		for(map<int,int>::iterator aIter=sizeDist->begin();aIter!=sizeDist->end();aIter++)
		{
			cout <<" " << aIter->first<<":"<< aIter->second;
		}
		cout << endl;
	}

	return 0;*/
} // end "generate"


int
Framework::redisplay(const char* outputDir)
{
	scMgr.initExperts();
	scMgr.showClusters(outputDir);
	//Only redisplay the data
	return 0;
}


int 
Framework::setPdiagonalLeaf(double aval)
{
	p_diagonal_leaf=aval;
	return 0;
}

int 
Framework::setPdiagonalNonLeaf(double aval)
{
	p_diagonal_nonleaf=aval;
	return 0;
}

//Need to initialize a conditional distribution for every branch which species
//the probability of transitioning from cluster k to cluster j
//TODO: initialize using the ribosomal clusters or some known cluster membership across species
//TODO: Use the rate and the branch length information
int 
Framework::initClusterTransitionProb()
{
	SpeciesDistManager::Species* root=sdMgr.getRoot();
	Matrix* conditional=root->getParams();
	int colcnt=conditional->getColCnt();
	for(int i=0;i<colcnt;i++)
	{
		double aval=1/((double)colcnt);
		conditional->setValue(aval,0,i);
	}
    for (int i=0; i<root->children.size(); i++) {
        initClusterTransitionProb(root->children[i]);
    }
	//initClusterTransitionProb(root->leftchild);
	//initClusterTransitionProb(root->rightchild);
	return 0;
}

int
Framework::initClusterTransitionProb(SpeciesDistManager::Species* anode)
{
	cout <<"Transitions for " << anode->name << endl;
	//if(anode->leftchild==NULL)
    if (anode->children.empty())
	{
		if(clusterTransitionProb.find(anode->name)==clusterTransitionProb.end())
		{
			cout <<"No cluster transition prob for  " << anode->name << endl;
			exit(0);
		}
		double pval=clusterTransitionProb[anode->name];
		initTransitionProb(anode->getParams(),pval);
	}
	else
	{
		if(clusterTransitionProb.find(anode->name)==clusterTransitionProb.end())
		{
			cout <<"No cluster transition prob for  " << anode->name << endl;
			exit(0);
		}
		double pval=clusterTransitionProb[anode->name];
		initTransitionProb(anode->getParams(),pval);
	}
    /*
	if(anode->leftchild!=NULL)
	{
		cout <<"Transitions for " << anode->leftchild->name << endl;
		initClusterTransitionProb(anode->leftchild);
		cout <<"Transitions for " << anode->rightchild->name << endl;
		initClusterTransitionProb(anode->rightchild);
	}*/
    for (int i=0; i<anode->children.size(); i++) {
        cout <<"Transitions for " << anode->children[i]->name << endl;
		initClusterTransitionProb(anode->children[i]);
    }
	return 0;
}

//The matrix is supposed to be a transition matrix of going from cluster k to cluster l
int
Framework::initTransitionProb(Matrix* m,double initval)
{
	int rowcnt=m->getRowCnt();
	int colcnt=m->getColCnt();
	for (int i=0;i<rowcnt;i++)
	{
		double s=0;
		for(int j=0;j<colcnt;j++)
		{
			double aval=0;
			if(i==j)
			{
				aval=initval;
			}
			else
			{
				aval=(1-initval)/((double) (rowcnt-1));
			}
			// returns random value between 0 and 0.01
			// randnum is the RNG
			double err=gsl_ran_flat(randnum,0,0.01);
			m->setValue(aval+err,i,j);
			s=s+err+aval;
		}
		for(int j=0;j<colcnt;j++)
		{
			double aval=m->getValue(i,j);	
			aval=aval/s;
			m->setValue(aval,i,j);
		}
	}
	m->showMatrix(1e-5);
	return 0;
}


int 
Framework::inferAncestralClusters(map<int,map<string,int>*>& clusterAssignments)
{
	cout <<"Inferring ancestral clusters" << endl;
	map<string,map<int,int>*> speciesClusterSizeDist;
	int disp=0;
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* extantClustering=cIter->second;
		map<string,int> ancestralClustering;
		sdMgr.getAncestralClustering(*extantClustering,ancestralClustering);
		for(map<string,int>::iterator aIter=ancestralClustering.begin();aIter!=ancestralClustering.end();aIter++)
		{
			(*extantClustering)[aIter->first]=aIter->second;
		}
		for(map<string,int>::iterator aIter=extantClustering->begin();aIter!=extantClustering->end();aIter++)
		{
			map<int,int>* clusterCnts=NULL;
			if(speciesClusterSizeDist.find(aIter->first)==speciesClusterSizeDist.end())
			{
				clusterCnts=new map<int,int>;
				speciesClusterSizeDist[aIter->first]=clusterCnts;
			}
			else
			{
				clusterCnts=speciesClusterSizeDist[aIter->first];
			}
			if(clusterCnts->find(aIter->second)==clusterCnts->end())
			{
				(*clusterCnts)[aIter->second]=1;
			}
			else
			{
				(*clusterCnts)[aIter->second]=(*clusterCnts)[aIter->second]+1;
			}
		}
		
		/*if(cIter==clusterAssignments.begin())
		{
			cout<<"Gene";
			for(map<string,int>::iterator aIter=extantClustering->begin();aIter!=extantClustering->end();aIter++)
			{
				cout <<" " << aIter->first;
			}
			cout << endl;
		}*/
		if(disp<10)
		{
			cout<<cIter->first;
			for(map<string,int>::iterator aIter=extantClustering->begin();aIter!=extantClustering->end();aIter++)
			{	
				cout <<" " << aIter->first<<"="<< aIter->second;
			} 
			cout << endl;
		}
		disp++;
	
	}
	for(map<string,map<int,int>*>::iterator sIter=speciesClusterSizeDist.begin();sIter!=speciesClusterSizeDist.end();sIter++)
	{
		cout <<sIter->first;
		map<int,int>* csize=sIter->second;
		for(map<int,int>::iterator cIter=csize->begin();cIter!=csize->end();cIter++)
		{
			cout <<" " << cIter->first <<":"<< cIter->second;
		}
		cout << endl;
	}
	return 0;
}

int
Framework::inferExtantClusters(map<int,map<string,int>*>& clusterAssignments)
{
	cout <<"Inferring extant clusters" << endl;
	map<string,map<int,int>*> speciesClusterSizeDist;
	int disp=0;
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* ancestralClustering=cIter->second;
		map<string,int> extantClustering;
		sdMgr.getExtantClustering(*ancestralClustering,extantClustering);
		for(map<string,int>::iterator aIter=extantClustering.begin();aIter!=extantClustering.end();aIter++)
		{
			(*ancestralClustering)[aIter->first]=aIter->second;
		}
		for(map<string,int>::iterator aIter=ancestralClustering->begin();aIter!=ancestralClustering->end();aIter++)
		{
			map<int,int>* clusterCnts=NULL;
			if(speciesClusterSizeDist.find(aIter->first)==speciesClusterSizeDist.end())
			{
				clusterCnts=new map<int,int>;
				speciesClusterSizeDist[aIter->first]=clusterCnts;
			}
			else
			{
				clusterCnts=speciesClusterSizeDist[aIter->first];
			}
			if(clusterCnts->find(aIter->second)==clusterCnts->end())
			{
				(*clusterCnts)[aIter->second]=1;
			}
			else
			{
				(*clusterCnts)[aIter->second]=(*clusterCnts)[aIter->second]+1;
			}
		}
		/*if(cIter==clusterAssignments.begin())
		{
			cout<<"Gene";
			for(map<string,int>::iterator aIter=ancestralClustering->begin();aIter!=ancestralClustering->end();aIter++)
			{
				cout <<" " << aIter->first;
			}
			cout << endl;
		}*/
		if(disp<10)
		{
			cout<<cIter->first;
			for(map<string,int>::iterator aIter=ancestralClustering->begin();aIter!=ancestralClustering->end();aIter++)
			{	
				cout <<" " << aIter->first<<"="<< aIter->second;
			} 
			cout << endl;
		}
		disp++;
	}
	for(map<string,map<int,int>*>::iterator sIter=speciesClusterSizeDist.begin();sIter!=speciesClusterSizeDist.end();sIter++)
	{
		cout <<sIter->first;
		map<int,int>* csize=sIter->second;
		for(map<int,int>::iterator cIter=csize->begin();cIter!=csize->end();cIter++)
		{
			cout <<" " << cIter->first <<":"<< cIter->second;
		}
		cout << endl;
	}
	return 0;
}

int 
Framework::estimateClusterTransProb(map<int,map<string,int>*>& clusterAssignments)
{
	SpeciesDistManager::Species* root=sdMgr.getRoot();
	//Estimate prior of the root
	Matrix* rootparam=root->conditional;
	rootparam->setAllValues(1e-5);
	for(map<int,map<string,int>*>::iterator aIter=clusterAssignments.begin();aIter!=clusterAssignments.end();aIter++)
	{
		int cid=(*aIter->second)[root->name];
		double aval=rootparam->getValue(cid,0);
		aval=aval+1;
		rootparam->setValue(aval,cid,0);
	}
	for(int i=0;i<rootparam->getRowCnt();i++)
	{
		double aval=rootparam->getValue(i,0)/clusterAssignments.size();
		rootparam->setValue(aval,i,0);
	}
	rootparam->showMatrix(1e-5);
	//estimateClusterTransProb(root,root->leftchild,clusterAssignments);
	//estimateClusterTransProb(root,root->rightchild,clusterAssignments);
    for (int i=0; i<root->children.size(); i++) {
        estimateClusterTransProb(root,root->children[i],clusterAssignments);
    }
	return 0;
}

int
Framework::estimateClusterTransProb(SpeciesDistManager::Species* parent, SpeciesDistManager::Species* child, map<int,map<string,int>*>& clusterAssignments)
{
	estimateTransitionMatrix(parent->name,child->name,child,clusterAssignments);
	/*if(child->leftchild!=NULL)
	{
		estimateClusterTransProb(child,child->leftchild,clusterAssignments);
		estimateClusterTransProb(child,child->rightchild,clusterAssignments);
	}*/
    for (int i=0; i<child->children.size(); i++) {
        estimateClusterTransProb(child,child->children[i],clusterAssignments);
    }

	return 0;
}

int
Framework::estimateTransitionMatrix(string& parentname,string& childname, SpeciesDistManager::Species* child, map<int,map<string,int>*>& clusterAssignments)
{
	Matrix* param=child->getParams();
	param->setAllValues(1);
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* assignments=cIter->second;
		int ancid=(*assignments)[parentname];
		int childid=(*assignments)[childname];
		double currval=param->getValue(ancid,childid);
		currval=currval+1;
		param->setValue(currval,ancid,childid);
	}
	cout <<"New Params for " << childname << " before norm " << endl;
	param->showMatrix(1e-5);
	//Row is for ancestral cluster id. Cols in a row must add to 1
	for(int i=0;i<param->getRowCnt();i++)
	{
		double sum=0;
		for(int j=0;j<param->getColCnt();j++)
		{
			sum=sum+param->getValue(i,j);
		}
		for(int j=0;j<param->getColCnt();j++)
		{
			double prob=param->getValue(i,j);
			prob=prob/sum;
			param->setValue(prob,i,j);
		}
	}
	cout <<"After normalzation" << childname << endl;
	param->showMatrix(1e-5);
	return 0;
}


bool
Framework::checkConvergence(map<int,map<string,int>*>& clusterAssignments,map<int,map<string,int>*>& oldclusterAssignments)
{
	bool convergence=false;
	int changes=0;
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* newassign=cIter->second;
		map<string,int>* oldassign=oldclusterAssignments[cIter->first];
		
		for(map<string,int>::iterator nIter=newassign->begin();nIter!=newassign->end();nIter++)
		{
			if(oldassign->find(nIter->first)==oldassign->end())
			{
				cout <<"Fatal error " << endl;
				exit(0);
			}
			int currassignval=(*oldassign)[nIter->first];
			if(currassignval!=nIter->second)
			{
				changes++;
			}
		}
	}
	if(changes==0)
	{
		convergence=true;
	}
	return convergence;

}

int
Framework::saveClusterAssignments(map<int,map<string,int>*>& clusterAssignments,map<int,map<string,int>*> &oldclusterAssignments)
{
	for(map<int,map<string,int>*>::iterator aIter=oldclusterAssignments.begin();aIter!=oldclusterAssignments.end();aIter++)
	{
		//delete aIter->second;
		//aIter->second->clear();
	}
	//oldclusterAssignments.clear();
	for(map<int,map<string,int>*>::iterator bIter=clusterAssignments.begin();bIter!=clusterAssignments.end();bIter++)
	{
		map<string,int>* newassign=NULL;
		if(oldclusterAssignments.find(bIter->first)==oldclusterAssignments.end())
		{
			newassign=new map<string,int>;
			oldclusterAssignments[bIter->first]=bIter->second;
		}
		else
		{
			newassign=oldclusterAssignments[bIter->first];
		}
		map<string,int>* oldassign=bIter->second;
		for(map<string,int>::iterator dIter=oldassign->begin();dIter!=oldassign->end();dIter++)
		{
			//cout <<"updating " << bIter->first <<" "<< dIter->first << " " << dIter->second << endl;
			(*newassign)[dIter->first]=dIter->second;
		}
		oldclusterAssignments[bIter->first]=newassign;
	}
	//clusterAssignments.clear();
	return 0;

}


int 
Framework::sampleAncestralCluster(gsl_rng* r,SpeciesDistManager::Species* root)
{
	int clustId=-1;
	double pval=gsl_ran_flat(r,0,1);
	Matrix* params=root->getParams();
	vector<int>* sortedClustIDs=root->getSortedClusterIDs(0);
	if(sortedClustIDs==NULL)
	{
		sortedClustIDs=new vector<int>;
		for(int i=0;i<params->getColCnt();i++)
		{
			sortedClustIDs->push_back(i);
		}
		for(int i=0;i<params->getColCnt();i++)
		{
			for(int j=i+1;j<params->getColCnt();j++)
			{
				double v1=params->getValue(0,(*sortedClustIDs)[i]);
				double v2=params->getValue(0,(*sortedClustIDs)[j]);
				if(v1<v2)
				{
					int oldval=(*sortedClustIDs)[i];
					(*sortedClustIDs)[i]=(*sortedClustIDs)[j];
					(*sortedClustIDs)[j]=oldval;
				}
			}
		}
		root->setSortedClusterIDs(0,sortedClustIDs);
	}
	double cdf=params->getValue(0,(*sortedClustIDs)[0]);
	clustId=0;
	while(pval>cdf)
	{
		clustId++;
		cdf=cdf+params->getValue(0,(*sortedClustIDs)[clustId]);
	}
	int actualClustID=(*sortedClustIDs)[clustId];
	return actualClustID;
}

int 
Framework::sampleChildCluster(gsl_rng* r,SpeciesDistManager::Species* child,int,map<string,int>& clusterAssign)
{
	int parentID=-1;
	if(clusterAssign.find(child->parent->name)==clusterAssign.end())
	{
		cout <<"No parent cluster id for " << child->parent->name << endl;
		exit(0);
	}
	if((strcmp(child->name.c_str(),"Scer")==0) || (strcmp(child->name.c_str(),"Sbay")==0))
	{
	//	cout <<"At child " << child->name << endl;
	}
	parentID=clusterAssign[child->parent->name];
	vector<int>* sortedClustIDs=child->getSortedClusterIDs(parentID);
	if(sortedClustIDs==NULL)
	{
		sortedClustIDs=new vector<int>;
		sortIndices(child->getParams(),parentID,sortedClustIDs);
		child->setSortedClusterIDs(parentID,sortedClustIDs);
	}
	int childID=sampleChildCluster(r,parentID,child->getParams(),sortedClustIDs);
	clusterAssign[child->name]=childID;
    /*
	if(child->leftchild!=NULL)
	{
		sampleChildCluster(r,child->leftchild,childID,clusterAssign);
		sampleChildCluster(r,child->rightchild,childID,clusterAssign);
	}*/
    for (int i=0; i<child->children.size(); i++) {
        sampleChildCluster(r,child->children[i],childID,clusterAssign);
    }
	return 0;
}

int 
Framework::sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>* sortedClustIDs)
{
	int clustId=-1;
	double pval=gsl_ran_flat(r,0,1);
	double cdf=params->getValue(parentClusterId,(*sortedClustIDs)[0]);
	clustId=0;
	while(pval>cdf)
	{
		clustId++;
		cdf=cdf+params->getValue(parentClusterId,(*sortedClustIDs)[clustId]);
	}
	int actualClustID=(*sortedClustIDs)[clustId];
	return actualClustID;
}
	
int
Framework::sortIndices(Matrix* params,int parentID,vector<int>* sortedClustIDs)
{
	int colCnt=params->getColCnt();
	for(int i=0;i<colCnt;i++)
	{
		sortedClustIDs->push_back(i);
	}
	for(int i=0;i<params->getColCnt();i++)
	{
		for(int j=i+1;j<params->getColCnt();j++)
		{
			double v1=params->getValue(parentID,(*sortedClustIDs)[i]);
			double v2=params->getValue(parentID,(*sortedClustIDs)[j]);
			if(v1<v2)
			{
				int oldval=(*sortedClustIDs)[i];
				(*sortedClustIDs)[i]=(*sortedClustIDs)[j];
				(*sortedClustIDs)[j]=oldval;
			}
		}
	}
	return 0;
}


int 
Framework::dumpInferredAssignments(const char* outputDir,const char* suffix )
{
	char aFName[1024];
	sprintf(aFName,"%s/clusterassign_multspecies_%s_%d.txt",outputDir,suffix,learnIter);
	ofstream oFile(aFName);
	for(map<string,map<string,int>*>::iterator oIter=clusterAssignments.begin();oIter!=clusterAssignments.end();oIter++)
	{
		oFile <<oIter->first;
		map<string,int>* assignments=oIter->second;
		for(map<string,int>::iterator aIter=assignments->begin();aIter!=assignments->end();aIter++)
		{
			oFile <<"\t" << aIter->first<<"=" << aIter->second;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int 
Framework::setClusterTransProb(double ctransProb)
{
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);
	for(int s=0;s<speciesList.size();s++)
	{
		clusterTransitionProb[speciesList[s]]=ctransProb;	
	}
	return 0;
}

int 
Framework::setClusterTransProb(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
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
		double pval;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				speciesName.append(tok);
			}
			else if(tokCnt==1)
			{
				pval=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		clusterTransitionProb[speciesName]=pval;
	}
	inFile.close();
	return 0;
}

/*
 * Sets constant covariance.
 *
 */
int
Framework::setConstCov(const char* val)
{
        constCov=atof(val);
        scMgr.setConstCov(constCov);
        return 0;
}

/*
* Chooses whether to init experts to own time points (true)
* or to source species (False, default option)
*/
int
Framework::setInitExpertsPerSpecies(bool setPerSpecies)
{
	scMgr.setInitExpertsPerSpecies(setPerSpecies);
	return 0;
}

int 
Framework::setMode(LearnMode lm, double r1, double r2, double r3)
{
	scMgr.setMode(lm,r1,r2,r3);
	return 0;
}

/**
* Args:
*	(1) celltype_order.txt file 
*	(2) OGIDs file (can be subset of expr data)
*	(3) list of OGIDs for regulators -- NOT USED
*	(4) number of modules (k)
*	(5) cell lineage tree file
*   (6) config filename
*   (7) random mode for clusters part: none|yes|<int>
*	(8) output directory
*	(9) running mode: learn, learnCV, generate, visualize
*	(10) name of source cell type/time point
*	(11) init type for transition probs: uniform or branchlength (I think only uniform applies here)
*	(12) p_diagonal_nonleaf : init prob of staying in same module
*	(13) optional: constant covariance value (useful depending on distro of data)
*	
*/

int
main(int argc, const char** argv)
{
	LearnMode learnMode;
	double rho1=0;
	double rho2=0;
	double rho3=0;
	// DC adds option to fix cluster variance
	// also to init clusters from self or from source species
	// also orthology map for cisregulatory elements

	//if(argc<13 || argc>15)
	if(argc<13 || argc>18)
	{
		cout <<"Usage: ./learnDRMN celltype_order genegroup null maxk celllineage config rand[none|yes|<int>] ";
		cout << "outputDir mode[learn|learnCV|learnCV:<int>:<int>:<int>|generate|visualize] srcnode inittype[uniform|branchlength] p_diagonal_nonleaf [const_cov(double), selfInit]" << endl;
		cout <<"Inittype is for specifying how the cluster transition probabilities will be initialized. " << endl;
		cout << "If inittype is uniform then set to everything to the same and if branchlength then set everything from the file" << endl;
		cout <<"[const_cov] is an optional argument to fix cluster variance." << endl;
		cout <<"[selfInit (string)] is an optional argument to initialize experts per species separately, rather than to source species params." << endl;
		cout << "Note: this version initializes cluster means to the source cell type." << endl;
		cout << "Cluster randomization options:" << endl;
		cout << "\tnone (keep init clusters)" << endl;
		cout << "\tyes (randomize init clusters with PID as seed)" << endl;
		cout << "\t<int> (user-provided seed used to RNG)" << endl;
		cout << "If plain learnCV, will run all folds." << endl;
		cout << "If learnCV:<foldInt>:<seedInt>:<nFolds>, will initialize RNG with seedInt and set up for fold foldInt out of nFolds." << endl;
		cout << "If learnCV:<foldInt>:<seedInt>, will initialize RNG with seedInt and set up for fold foldInt out of 3 folds." << endl;
		return 0;
	}
    
    time_t result = time(0);
	cout << "Starting program: " << std::asctime(localtime(&result)) << endl;
    

    clock_t t=clock();  // time the whole program

	Framework fw;
    //k and tree
	fw.readSpeciesTree(atoi(argv[4]),argv[5]);
    //Order and Mapping
	int success = fw.readOrthology(argv[1],argv[2]);
	if (success !=0)
	{
		cerr << "Problem with reading orthology maps." << endl;
		return 1;
	}

	// we never use this
	/*success=fw.readRegulatorOGIds(argv[3]);
	if (success !=0 )
	{
		cerr << "Problem reading cisregulatory element OGID list." << endl;
		return 1;
	}*/ 

	fw.setSrcSpecies(argv[10]);
	// constant covariance specified? 
	// init experts to species-specific data?
	// iterate over optional arg choices
	for (int i=13; i<argc; i++) 
	{
		if (strcmp(argv[i],"selfInit")==0)
		{
			fw.setInitExpertsPerSpecies(true);
			cout << "We will initialize experts separately per species." << endl;
		}
		//else //if (argc>13)
		//{
		//	cerr << "SHOULDNT GET HERE!" << endl;
		//	int success=fw.setConstCov(argv[i]);
		//}
	}
	if (argc>=15)
	{
		cout << "here we go:" << argv[14] << endl;
		if (strcmp(argv[14],"LEASTL21")==0)
		{
			cout << "we got LEASTL21" << endl;
			learnMode = LEASTL21;
		}
		else if (strcmp(argv[14],"LEASTDIRTY")==0)
		{
			cout << "we got LEASTDIRTY" << endl;
			learnMode = LEASTDIRTY;
		}
		else if (strcmp(argv[14],"LEASTFUSED")==0)
		{
			cout << "we got LEASTFUSED" << endl;
			learnMode = LEASTFUSED;
		}
		else
		{
			cout << "we got GREEDY" << endl;
			learnMode = GREEDY;
		}
		//if (learnMode == LEASTDIRTY || learnMode == LEASTL21)
		if (learnMode != GREEDY)
		{
			if (argc>=16)
			{
				rho1 = atof(argv[15]);
			}
			if (argc>=17)
			{
				rho2 = atof(argv[16]);
			}
			if (argc>=18)
			{
				rho3 = atof(argv[17]);
			}
			fw.setMode(learnMode,rho1,rho2,rho3);
		}
	}

	if(strcmp(argv[11],"uniform")==0)
	{
		fw.setClusterTransProb(atof(argv[12]));
	}
	else 
	{
		fw.setClusterTransProb(argv[12]);
	}

	// argv9 can be: learn, learnCV, generate, or display.
	if(strcmp(argv[9],"learn")==0)
	{
    	//  RNG gets set in here
		success=fw.readSpeciesData(argv[6],argv[7]);
		if (success !=0)
		{
			cerr << "Problem reading cell-type specific data from config file." << endl;
			return 1;
		}

        success=fw.startClustering(argv[8]);
        if (success !=0)
		{
			cerr << "Problem while running fw.startClustering." << endl;
			return 1;
		}
	}
	//else if(strcmp(argv[9],"learnCV")==0)
	// starts with learnCV
	else if (strncmp(argv[9],"learnCV", strlen("learnCV"))==0)
	{

		int myFold=-1; // do all folds by default
		int foldSeed=-1; // use standard RNG behavior by default
		int foldCV=3; // by default do 3 folds

		// did we encode the fold ID and/or fold seed here?
		char* mytok=0;
		mytok=strtok(strdup(argv[9]),":");
		int tokID=0;
		while (mytok != NULL)
		{
			if (tokID==1) myFold=atoi(mytok);
			else if (tokID==2) foldSeed=atoi(mytok);
			else if (tokID==3) foldCV=atoi(mytok);
			tokID++;
			mytok=strtok(NULL,":");
		}

		cout << "We chose to do fold " << myFold << " with seed " << foldSeed << endl;
		
		//At this stage we will first split the OGIDs into train and test sets. Set the training set OGIDs before we call readSpeciesData. Then do the usual stuff
		//for now hardcode the number of folds of CV
		if (myFold >= foldCV)
		{
			cout << "Request fold ID " << myFold << " exceeds number of folds: " << foldCV << endl;
			return 2;

		}
		success=fw.startClustering_CV(argv[6],argv[7],argv[8],foldCV, myFold, foldSeed);
		if (success !=0)
		{
			cerr << "Problem while running fw.startClustering_CV." << endl;
			return 1;
		} 
	}
	else if(strcmp(argv[9],"generate")==0)
	{
		fw.generateData(argv[8]);
	}
	else 
	{
		cout << "Unrecognized running option: " << argv[9] << endl;
		return 2;
		//fw.redisplay(argv[8]);
	}

	t=clock()-t;
	cout << "Total run time (m) = " << ((long double) t/CLOCKS_PER_SEC)/60.0 << endl;
	result=time(0);
	cout << "Ending time: " << std::asctime(localtime(&result)) << endl;

	return 0;
}

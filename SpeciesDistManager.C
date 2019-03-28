
/*
Arboretum: An algorithm to cluster functional genomesomics data from multiple species
    Copyright (C) 2013 Sushmita Roy sushroy@gmail.com

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
#include <iostream>
#include <fstream>
#include <iostream>
#include <string.h>
#include "Matrix.H"
#include "SpeciesDistManager.H"

SpeciesDistManager::SpeciesDistManager()
{
	root=NULL;
	maxClusterCnt=0;
}

SpeciesDistManager::~SpeciesDistManager()
{
}

int
SpeciesDistManager::setMaxClusters(int k)
{
	maxClusterCnt=k;
	return 0;
}

int 
SpeciesDistManager::readSpeciesTree(const char* aFName)
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
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string childSpeciesName;
		string parentSpeciesName;
		//string childtype;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				childSpeciesName.append(tok);	
			}
			//else if(tokCnt==1)
			//{
			//	childtype.append(tok);
			//}
			else if(tokCnt==1)
			{
				parentSpeciesName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		SpeciesDistManager::Species* childSpecies=NULL;
		SpeciesDistManager::Species* parentSpecies=NULL;
		if(strcmp(childSpeciesName.c_str(),"NULL")!=0) 
		{
			if(speciesSet.find(childSpeciesName)==speciesSet.end())
			{
				childSpecies=new SpeciesDistManager::Species;
				childSpecies->name.append(childSpeciesName.c_str());
				childSpecies->parent=NULL;
				//childSpecies->children=NULL;
				speciesSet[childSpeciesName]=childSpecies;
				childSpecies->conditional=new Matrix(maxClusterCnt,maxClusterCnt);
			}
			else
			{
				childSpecies=speciesSet[childSpeciesName];
				if(childSpecies->conditional==NULL)
				{
					childSpecies->conditional=new Matrix(maxClusterCnt,maxClusterCnt);
				}
			}
		}
		if(speciesSet.find(parentSpeciesName)==speciesSet.end())
		{
			parentSpecies=new SpeciesDistManager::Species;
			parentSpecies->name.append(parentSpeciesName.c_str());
			parentSpecies->parent=NULL;
			//parentSpecies->children=NULL;
			speciesSet[parentSpeciesName]=parentSpecies;
		}
		else
		{
			parentSpecies=speciesSet[parentSpeciesName];
		}
        	if(childSpecies!=NULL)
		{
			childSpecies->parent=parentSpecies;
			parentSpecies->children.push_back(childSpecies);
		}
			
		if(parentSpecies->parent==NULL)
		{
			root=parentSpecies;
		}
	}
	root->conditional=new Matrix(1,maxClusterCnt);

	inFile.close();
	return 0;
}


int
SpeciesDistManager::assignLevel()
{
	levelFromRoot[root->name]=0;
	//assignLevel(root->leftchild,1);
	//assignLevel(root->rightchild,1);
	for (int i=0; i<root->children.size(); i++) {
        assignLevel(root->children[i],1);
    }
    return 0;
}

int
SpeciesDistManager::assignLevel(SpeciesDistManager::Species* node,int level)
{
	levelFromRoot[node->name]=level;
	
	for (int i=0; i<node->children.size(); i++)
    {
        if(node->children[i]!=NULL)
            assignLevel(node->children[i],level+1);
    }
	return 0;
}


int
SpeciesDistManager::getLevelFromRoot(const char* nodeName)
{
	string nKey(nodeName);
	if(levelFromRoot.find(nKey)==levelFromRoot.end())
	{	
		cout <<"No node with name " << nKey << endl;
		return -1;
	}
	return levelFromRoot[nKey];
}


SpeciesDistManager::Species*
SpeciesDistManager::getRoot()
{
	return root;
}

int 
SpeciesDistManager::getSpeciesListPrefix(vector<string>& specList)
{
	getSpeciesListPrefix(root,specList);
	return 0;
}


map<string,SpeciesDistManager::Species*>& 
SpeciesDistManager::getAllSpecies()
{	
	return speciesSet;
}

SpeciesDistManager::Species*
SpeciesDistManager::getSpecies(string& specKey)
{
	if(speciesSet.find(specKey)==speciesSet.end())
	{
		return NULL;
	}
	return speciesSet[specKey];
}

int
SpeciesDistManager::getSpeciesListPrefix(SpeciesDistManager::Species* node, vector<string>& specList)
{
	for (int i=0; i<node->children.size(); i++)
    {
        if(node->children[i]!=NULL)
		getSpeciesListPrefix(node->children[i],specList);
	}
    specList.push_back(node->name);
	return 0;
}

int 
SpeciesDistManager::resetTransitionProbability()
{
	root->conditional->setAllValues(0.001);
	
    for (int i=0; i<root->children.size(); i++)
    {
        resetTransitionProbability(root->children[i]);
	}
	return 0;
}

int
SpeciesDistManager::resetTransitionProbability(SpeciesDistManager::Species* species)
{
	species->conditional->setAllValues(0.001);
	
    for (int i=0; i<species->children.size(); i++)
    {
        if(species->children[i]!=NULL)
        resetTransitionProbability(species->children[i]);
	}
	return 0;
}

int 
SpeciesDistManager::normalizeTransitionMatrix()
{
	normalizeTransitionMatrix(root);
	return 0;
}

int
SpeciesDistManager::normalizeTransitionMatrix(SpeciesDistManager::Species* node)
{
	Matrix* conditional=node->conditional;
	double globalsum=0;
	//cout <<"Transition for " << node->name << " before norm" << endl;
	//node->conditional->showMatrix(); // DEBUG
	for(int r=0;r<conditional->getRowCnt();r++)
	{
		double sum=0;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double aval=conditional->getValue(r,c);
			sum=sum+aval;
		}
		globalsum=globalsum+sum;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double prob=conditional->getValue(r,c);
			prob=prob/sum;
			conditional->setValue(prob,r,c);
		}
	}
	//cout <<"Transition for " << node->name << endl;
	//node->conditional->showMatrix();
    
	
    for (int i=0; i<node->children.size(); i++)
    {
        if(node->children[i]!=NULL)
	{
           // resetTransitionProbability(node->children[i]);
           // SR:fix
	   normalizeTransitionMatrix(node->children[i]);
	}
    }
	return 0;
}

int
SpeciesDistManager::initTransitionMatrix_ML()
{
	//Start at the root and move downwards
	root->conditional_ml=new Matrix(1,maxClusterCnt);
	root->conditional_ml->setAllValues(0.001);
	
    for (int i=0; i<root->children.size(); i++)
    {
        if(root->children[i]!=NULL)
            initTransitionMatrix_ML(root->children[i]);
	}
	return 0;
}

int
SpeciesDistManager::initTransitionMatrix_ML(SpeciesDistManager::Species* node)
{
	node->conditional_ml=new Matrix(maxClusterCnt,maxClusterCnt);
	node->conditional_ml->setAllValues(0.001);
	
    for (int i=0; i<node->children.size(); i++)
    {
        if(node->children[i]!=NULL)
            initTransitionMatrix_ML(node->children[i]);
	}
	return 0;
}

int 
SpeciesDistManager::normalizeTransitionMatrix_ML()
{
	normalizeTransitionMatrix_ML(root);
	return 0;
}

int
SpeciesDistManager::normalizeTransitionMatrix_ML(SpeciesDistManager::Species* node)
{
	Matrix* conditional=node->conditional_ml;
	double globalsum=0;
	//cout <<"Transition for " << node->name << " before norm" << endl;
	//conditional->showMatrix();
	for(int r=0;r<conditional->getRowCnt();r++)
	{
		double sum=0;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double aval=conditional->getValue(r,c);
			sum=sum+aval;
		}
		globalsum=globalsum+sum;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double prob=conditional->getValue(r,c);
			prob=prob/sum;
			conditional->setValue(prob,r,c);
		}
	}
	//cout <<"Transition for " << node->name << endl;
	//conditional->showMatrix();
	
    for (int i=0; i<node->children.size(); i++)
    {
        if(node->children[i]!=NULL)
            normalizeTransitionMatrix_ML(node->children[i]);
	}
	return 0;
}


Matrix*
SpeciesDistManager::getTransitionMatrix_ML(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	return species->conditional_ml;
}




Matrix*
SpeciesDistManager::getConditional(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	return species->conditional;
}

//The probability that a gen maintains or switches its module assignment from its parent
double 
SpeciesDistManager::getConditionalProb(string& spName, int parentCluster,int childCluster)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	double prob=species->conditional->getValue(parentCluster,childCluster);
	return prob;
}

//The probability that an edge is maintained in a child given that the edge is present in the ancestor to be computed in a module-specific manner.
double 
SpeciesDistManager::getConditionalProbForEdge(string& spName, int parentEdgeStatus,int childEdgeStatus,int moduleID)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	//We always only care about the diagonal element
	double prob=species->conditional->getValue(moduleID,moduleID);
	if(parentEdgeStatus!=childEdgeStatus)
	{
		prob=1-prob;			
	}
	return prob;
}

/*
* For one cluster, print tree with prob that edge is maintained in child given edge in parent.
* DC add for DRMN
* prints out:
*  child \t parent \t maintenance_prob
*/
int
SpeciesDistManager::showTreeForCluster(int clustID)
{
	// start at root
    	for (int i=0; i<root->children.size(); i++)
	{
        	showTreeForCluster(root->children[i], clustID, i);
	}
	return 0;
}

/*
* For one cluster, print tree with prob that edge is maintained in child given edge in parent.
* DC add for DRMN
* prints out:
*  child \t parent \t maintenance_prob
*/
int
SpeciesDistManager::printTreeForCluster(int clustID, const char* outFName)
{
	//open file
	ofstream ofile(outFName);
	// start at root
    	for (int i=0; i<root->children.size(); i++)
	{
        	printTreeForCluster(root->children[i], clustID, i, ofile);
	}

	ofile.close();

	return 0;
}

/**
* Helper function for showTreeForCluster(int clustID)
* TO DO: Update MRTLE code to handle arbitrary number of children
*/
int 
SpeciesDistManager::showTreeForCluster(Species* node, int clustID, int childID)
{	
	// get the probability of maintaining an edge in same cluster between parent and child
	double prob = getConditionalProb(node->name, clustID,clustID);

	// For MRTLE, need to specify left or right child (need to update this!!)
	string childtype;
	if (childID==0)
	{
		childtype.append("left");
	} else if(childID==1)
	{
		childtype.append("right");	
	} else 
	{	
		cerr << "SORRY, nonbinary tree not implemented yet" << endl;
	}

	// child \t parent \t gain \t loss
	// for now, let's say gain==loss==(1-maintenance)
	cout << node->name << "\t" << childtype << "\t" << node->parent->name << "\t" << (1.0-prob) << "\t" << (1.0-prob) << endl;

	for (int i=0; i<node->children.size(); i++)
	{
        	showTreeForCluster(node->children[i], clustID, i);
	}
	return 0;
}

/**
* Helper function for showTreeForCluster(int clustID)
* TO DO: Update MRTLE code to handle arbitrary number of children
*/
int 
SpeciesDistManager::printTreeForCluster(Species* node, int clustID, int childID, ofstream& ofile)
{	
	// get the probability of maintaining an edge in same cluster between parent and child
	double prob = getConditionalProb(node->name, clustID,clustID);

	// For MRTLE, need to specify left or right child (need to update this!!)
	string childtype;
	if (childID==0)
	{
		childtype.append("left");
	} else if(childID==1)
	{
		childtype.append("right");	
	} else 
	{	
		cerr << "SORRY, nonbinary tree not implemented yet" << endl;
	}

	// child \t parent \t gain \t loss
	// for now, let's say gain==loss==(1-maintenance)
	ofile << node->name << "\t" << childtype << "\t" << node->parent->name << "\t" << (1.0-prob) << "\t" << (1.0-prob) << endl;

	for (int i=0; i<node->children.size(); i++)
	{
        	printTreeForCluster(node->children[i], clustID, i, ofile);
	}
	return 0;
}


double
SpeciesDistManager::getEdgeStatusProb(map<string,int>& edgeStatus, int clusterID)
{
	/*for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
	{
		cout <<eIter->first <<"\t" << eIter->second << endl;
	}*/
	double edgeprior=0;
	//Start with the root
	int edgeInRoot=edgeStatus[root->name];
	//Don't really know what the probability of edge presence and absence should be at the root node, 
	//but let's assume that the probability of edge to be associated with a module is dependent on the 
	//prior probability of the module
	Matrix* probInRoot=root->conditional;
	double score=probInRoot->getValue(0,clusterID);
        for(int j=0; j<root->children.size(); j++)
        {
		score*=getSubTree(edgeInRoot,root->children[j],edgeStatus,clusterID);
       	}
	edgeprior=score;
	return edgeprior;
}


int 
SpeciesDistManager::showInferredConditionals(const char* outputDir)
{
	showConditionals(outputDir,root);
	return 0;
}

int
SpeciesDistManager::showConditionals(const char* outputDir,SpeciesDistManager::Species* species)
{
	char output[1024];
	sprintf(output,"%s/%s_transprob.txt",outputDir,species->name.c_str());
	ofstream oFile(output);	
	species->conditional->showMatrix(oFile);
	//if(species->leftchild!=NULL)
	//{
	for (int j=0; j<species->children.size(); j++)
        {
            showConditionals(outputDir,species->children[j]);
        }
        //showConditionals(outputDir,species->leftchild);
		//showConditionals(outputDir,species->rightchild);
	//}
	return 0;
}

int 
SpeciesDistManager::showInferredConditionals_ML(const char* outputDir)
{
	showConditionals_ML(outputDir,root);
	return 0;
}

int
SpeciesDistManager::showConditionals_ML(const char* outputDir,SpeciesDistManager::Species* species)
{
	char output[1024];
	sprintf(output,"%s/ml_%s",outputDir,species->name.c_str());
	ofstream oFile(output);	
	species->conditional_ml->showMatrix(oFile);
	//if(species->leftchild!=NULL)
	//{
        for (int j=0; j<species->children.size(); j++)
        {
            showConditionals_ML(outputDir,species->children[j]);
        }
        //showConditionals_ML(outputDir,species->leftchild);
		//showConditionals_ML(outputDir,species->rightchild);
	//}
	return 0;
}


double
SpeciesDistManager::getSubTree(int parentEdgeStatus, Species* child, map<string,int>& edgeStatus, int moduleID)
{
	double score=0;
	int myEdgeStatus=edgeStatus[child->name];
	//if(child->leftchild==NULL && child->rightchild==NULL)
	if(child->children.empty())
	{
		//This is a leaf node
		score=getConditionalProbForEdge(child->name,parentEdgeStatus,myEdgeStatus,moduleID);
	}
	else
	{
		//Since all cell types are observed, we don't need to worry about marginalization 
		double children_score=1;
		for (int j=0; j<child->children.size(); j++)
		{	
			children_score*=getSubTree(myEdgeStatus,child->children[j],edgeStatus,moduleID);
		}
		double prob=getConditionalProbForEdge(child->name,parentEdgeStatus,myEdgeStatus,moduleID);
		score=(prob*children_score);
	}
	return score;
}

double
SpeciesDistManager::getAncestralClustering(map<string,int>& extantClusterAssign,map<string,int>& ancestralClusterAssign)
{
	//Start with the root
	int maxParentCluster=-1;
	double maxScore=0;
	map<int,double> assignProb;
	map<string,int> tempAssign;
	for(int i=0;i<maxClusterCnt;i++)
	{
        double children_score=1;
        for (int j=0; j<root->children.size(); j++)
        {
            children_score*=maxSubTree(i,root->children[j],extantClusterAssign,tempAssign);
        }
        //double leftScore=maxSubTree(i,root->leftchild,extantClusterAssign,tempAssign);
		//double rightScore=maxSubTree(i,root->rightchild,extantClusterAssign,tempAssign);
		double prior=root->conditional->getValue(i,0);
		//double score=prior*leftScore*rightScore;
        double score=prior*children_score;
		assignProb[i]=score;
		if(score>maxScore)	
		{
			maxScore=score;
			maxParentCluster=i;
			for(map<string,int>::iterator aIter=tempAssign.begin();aIter!=tempAssign.end();aIter++)
			{
				ancestralClusterAssign[aIter->first]=aIter->second;
			}
		}
		/*cout << score <<" " << maxScore;
		for(map<string,int>::iterator aIter=ancestralClusterAssign.begin();aIter!=ancestralClusterAssign.end();aIter++)
		{
			cout <<" " << aIter->first<<"=" << aIter->second;
		}
		cout << endl;*/
	}
	tempAssign.clear();
	for(map<int,double>::iterator aIter=assignProb.begin();aIter!=assignProb.end();aIter++)
	{
		if(aIter->second>=maxScore && aIter->first!=maxParentCluster)
		{
			cout<< root->name<<" " <<aIter->first <<"=" << aIter->second << endl;
		}
	}
	ancestralClusterAssign[root->name]=maxParentCluster;
	return maxScore;
}

double 
SpeciesDistManager::getExtantClustering(map<string,int>& ancestralClusterAssign, map<string,int>& extantClusterAssign)
{
	int maxParentCluster=-1;
	assignExtantClustering(ancestralClusterAssign,root,extantClusterAssign);
	return 0;
}

int
SpeciesDistManager::assignExtantClustering(map<string,int>& ancestralClusterAssign,Species* node, map<string,int>& extantClusterAssign)
{
	if(node->parent!=NULL && node->children.empty())
	{
		Species* parent=node->parent;
		int clusterid=ancestralClusterAssign[parent->name];
		int maxchildclusterid=getMaxClusterAssignForChild(clusterid,node);
		extantClusterAssign[node->name]=maxchildclusterid;
	}
	else
	{
		for (int j=0; j<node->children.size(); j++)
        	{	
  	          assignExtantClustering(ancestralClusterAssign,node->children[j],extantClusterAssign);
        	}
        
      		  //assignExtantClustering(ancestralClusterAssign,node->leftchild,extantClusterAssign);
		//assignExtantClustering(ancestralClusterAssign,node->rightchild,extantClusterAssign);
	}
	return 0;
}

int
SpeciesDistManager::getMaxClusterAssignForChild(int parentClustID,Species* child)
{
	double maxchildprob=0;
	int maxchildcluster=-1;
	for(int i=0;i<maxClusterCnt;i++)
	{
		double prob=getConditionalProb(child->name,parentClustID,i);
		if(prob>maxchildprob)
		{
			maxchildprob=prob;
			maxchildcluster=i;
		}
	}
	return maxchildcluster;
}



double
SpeciesDistManager::scoreAssignment(map<string,int>& jointAssign)
{
	double assignpdf=1;
	for(map<string,int>::iterator aIter=jointAssign.begin();aIter!=jointAssign.end();aIter++)
	{
		Species* node=speciesSet[aIter->first];
		int clusterid=aIter->second;
		if(node->parent==NULL)
		{
			assignpdf=assignpdf*root->conditional->getValue(clusterid,0);;
		}
		else
		{
			Species* parent=node->parent;
			if(jointAssign.find(parent->name)==jointAssign.end())
			{
				cout <<"No assignment for " << parent->name << endl;
			}
			int parentAssign=jointAssign[parent->name];
			double score=getConditionalProb((string&)aIter->first,parentAssign,aIter->second);
			assignpdf=assignpdf*score;
		}
	}
	return assignpdf;
}


double
SpeciesDistManager::maxSubTree(int parentCluster, Species* child, map<string,int>& extantCluster,map<string,int>& ancestralCluster)
{
	double score=0;
	//if(child->leftchild==NULL && child->rightchild==NULL)
    if(child->children.empty())
	{
		//This is a leaf node
		int clusterassign=extantCluster[child->name];
		score=getConditionalProb(child->name,parentCluster,clusterassign);
	}
	else
	{
		//This is an interior node
		int maxParentCluster=-1;
		double maxScore=0;
		map<int,double> assignProb;
		for(int i=0;i<maxClusterCnt;i++)
		{
			double children_score=1;
            for (int j=0; j<child->children.size(); j++)
            {
                children_score*=maxSubTree(i,child->children[j],extantCluster,ancestralCluster);
            }
            //double leftScore=maxSubTree(i,child->leftchild,extantCluster,ancestralCluster);
			//double rightScore=maxSubTree(i,child->rightchild,extantCluster,ancestralCluster);
			
            double prob=getConditionalProb(child->name,parentCluster,i);
			//double currScore=prob*leftScore*rightScore;
            double currScore=prob*children_score;
			if(currScore>maxScore)
			{
				maxScore=currScore;
				maxParentCluster=i;
			}
			assignProb[i]=currScore;
		}
		for(map<int,double>::iterator aIter=assignProb.begin();aIter!=assignProb.end();aIter++)
		{
			
		}
		score=maxScore;
		if(strcmp(child->name.c_str(),"Anc3")==0)
		{
	//		cout <<"Update: " << child->parent->name<<"="<< parentCluster<<" "<< child->name<<"=" << maxParentCluster <<" Score=" <<maxScore << endl;
		}
		ancestralCluster[child->name]=maxParentCluster;
		assignProb.clear();
	}
	return score;
}



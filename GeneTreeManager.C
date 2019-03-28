
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

#include <iostream>
#include <fstream>
#include <map>
#include <string.h>
#include <stdlib.h>

#include "SpeciesDistManager.H"
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "GeneTree.H"
//#include "NewickReader.H"
#include "GeneTreeManager.H"

GeneTreeManager::GeneTreeManager()
{
	totalAdditionalParents=0;
}

GeneTreeManager::~GeneTreeManager()
{

}

int 
GeneTreeManager::setSpeciesTree(SpeciesDistManager::Species* node)
{
	speciestreeRoot=node;
	return 0;
}

int
GeneTreeManager::setSpeciesDistManager(SpeciesDistManager* aPtr)
{
	sdMgr=aPtr;
	return 0;
}

int
GeneTreeManager::setGeneTreeDir(const char* treeDirName)
{
	strcpy(dirName,treeDirName);
	return 0;
}



//We will assume that there is no more than 1 duplication event
GeneTree*
GeneTreeManager::generateTree(MappedOrthogroup* mor)
{
	//We will descend down the species tree until we reach the leaf nodes.
	//At the leaf node we will ask how many copies/genes does the species node have
	SpeciesDistManager::Species* dupSpecies=NULL;
	identifyAncestorWithDuplicates(speciestreeRoot,mor,&dupSpecies);
    dupSpecies=NULL;
	GeneTree* tree=createTree(dupSpecies,speciestreeRoot,mor,1);
	int ogid=mor->getID();
	genetreeSet[ogid]=tree;
	if(ogid==1339|| ogid==13692)
	{
		//showTree(tree);
	}
	return tree;
}
/**********************************
GeneTree*
GeneTreeManager::generateTreeFromFile(MappedOrthogroup* mor)
{
	int ogid=mor->getID();
	char treePath[1024];
	sprintf(treePath,"%s/%d.tre",GENETREEPATH,ogid);
	if(ogid==1339|| ogid==13692 || ogid==4597)
	{
		cout << "Stop here" << endl;
	}
	GeneTree* tree=readTreeFromFile(treePath);
	//cout << "Pruning " << ogid  << endl;
	//showTree(tree);
	pruneTree(tree);
	//cout << "After pruning " << ogid << endl;
	//showTree(tree);
	//getSpeciesName(tree);
	map<string,int> ancestors;
	getSpeciesName(tree,ancestors);
	//cout <<"After renaming "<< ogid << endl;
	//showTree(tree);
	//Move the root down if we have not been able to find a match
	while(strstr(tree->name.c_str(),"Internal")!=NULL)
	{
		//if(tree->leftchild!=NULL && tree->rightchild!=NULL)
        
        if(tree->children.empty())
		{
			cout << "Could not rename root " << tree->name << " of " << ogid << " even though children are present" << endl;
			exit(0);
		}
        /*
		if(tree->leftchild!=NULL)
		{
			tree=tree->leftchild;
			tree->parent=NULL;
		}
		else if(tree->rightchild!=NULL)
		{
			tree=tree->rightchild;
			tree->parent=NULL;
		}
        tree=tree->children[0];
        tree->parent=NULL;
	}
//	cout <<"Now inserting " << ogid << endl;
	insertTree(tree);
//	showTree(tree);
	//Need to make sure that we get to the parent
	while(tree->parent!=NULL)
	{
		tree=tree->parent;
	}
	if(ogid==1339|| ogid==13692 || ogid==4597)
	{
		showTree(tree);
	}
	SpeciesDistManager::Species* specroot=sdMgr->getRoot();
	//Add additional parents
	bool addAdditional=false;
	string currentRoot(tree->species);
	while(strcmp(tree->species.c_str(),specroot->name.c_str())!=0)
	{
		SpeciesDistManager::Species* currSpecies=sdMgr->getSpecies(tree->species);
		SpeciesDistManager::Species* parentSpecies=currSpecies->parent;
		GeneTree* anode=new GeneTree;
		anode->exclude=true;
		anode->species.append(parentSpecies->name.c_str());
		anode->name.append(parentSpecies->name.c_str());
        /*
		if(parentSpecies->leftchild==currSpecies)
		{
			anode->leftchild=tree;
			anode->rightchild=NULL;
		}
		else
		{
			anode->rightchild=tree;
			anode->leftchild=NULL;
		}
 
        vector<SpeciesDistManager::Species>::it=find(parentSpecies.begin(),parentSpecies.end(),currSpecies);
        anode->children.push_back(tree);
        
		tree->parent=anode;
		tree=anode;
		addAdditional=true;
	}
	if(addAdditional)
	{
		if(addedParents.find(currentRoot)==addedParents.end())
		{
			addedParents[currentRoot]=1;
		}	
		else
		{
			addedParents[currentRoot]=addedParents[currentRoot]+1;
		}
		totalAdditionalParents++;
	}
	setNodeType(tree);
	//showTree(tree);
	//Check if there are any Kpol genes here
	GeneMap* kpol=mor->getSpeciesHits("Kpol");
	if(kpol!=NULL)
	{
		map<string,map<string,STRINTMAP*>*>& kpolGenes=kpol->getGeneSet();
		if(kpolGenes.size()>1)
		{
			cout <<"Too many kpol genes in OG" << ogid << "!" << endl;
			return tree;
		}
		addKpolGene(tree,(string&)kpolGenes.begin()->first);
	}
	return tree;
}
***********************************/
int
GeneTreeManager::showTotalAdditionalParentDist()
{
	for(map<string,int>::iterator pIter=addedParents.begin();pIter!=addedParents.end();pIter++)
	{
		cout << pIter->first << "\t" << pIter->second << endl;
	}
	return 0;
}

int 
GeneTreeManager::identifyAncestorWithDuplicates(SpeciesDistManager::Species* speciesnode,MappedOrthogroup* mor,SpeciesDistManager::Species** duplicateSpecies)
{
	//if(speciesnode->leftchild==NULL)
    if(speciesnode->children.empty())
	{
		GeneMap* genemap=mor->getSpeciesHits(speciesnode->name.c_str());
		if(genemap==NULL)
		{
			return 0;
		}
		int copyno=genemap->getGeneSet().size();
		if(copyno==2)
		{
			*duplicateSpecies=speciesnode;
		}
		return copyno;
	}
	//int leftcopyno=identifyAncestorWithDuplicates(speciesnode->leftchild,mor,duplicateSpecies);
	//int rightcopyno=identifyAncestorWithDuplicates(speciesnode->rightchild,mor,duplicateSpecies);
	//if(leftcopyno>=2 || rightcopyno>=2)
    int flag=0;
    for (int i=0; i<speciesnode->children.size(); i++) {
        int rightcopyno=identifyAncestorWithDuplicates(speciesnode->children[i],mor,duplicateSpecies);
        if (rightcopyno>=2) {
            flag=1;
        }
        
    }
    if(flag==1)
	{
		*duplicateSpecies=speciesnode;
		return 2;
	}
	
	return 1;
}

int
GeneTreeManager::showTree(GeneTree* anode)
{
	/*GeneTree* lchild=anode->leftchild;
        if(lchild!=NULL)
        {
		if(anode->nodeType==1)
		{
                	cout <<anode->species <<"|" <<anode->name<<"\tleft\t" << lchild->species<<"|" << lchild->name << endl;
		}
		else if(anode->nodeType==2)
		{
                	cout << "(Duplicate) " <<anode->species <<"|" <<anode->name<<"\tleft\t" << lchild->species<<"|" << lchild->name << endl;
		}
		else
		{
			cout <<"Unrecognized node type " << anode->nodeType << endl;
		}
        }
        GeneTree* rchild=anode->rightchild;
        if(rchild!=NULL)
        {
		if(anode->nodeType==1)
		{
                	cout <<anode->species <<"|" <<anode->name<<"\tright\t" << rchild->species<<"|" << rchild->name << endl;
		}
		else if(anode->nodeType==2)
		{
			cout <<"(Duplicate) "  <<anode->species <<"|" <<anode->name<<"\tleft\t" << lchild->species<<"|" << lchild->name << endl; 
		}
		else
		{
			cout <<"Unrecognized node type " << anode->nodeType << " for " << anode->species << " | " << anode->name << endl;
		}
        }
        if(lchild!=NULL)
        {
                showTree(lchild);
        }
        if(rchild!=NULL)
        {
                showTree(rchild);
        }
     */
    for (int i=0; i<anode->children.size(); i++) {
        if(anode->children[i]!=NULL)
        {
            if(anode->nodeType==1)
            {
                cout <<anode->species <<"|" <<anode->name<<"\t"<<i<<"\t" << anode->children[i]->species<<"|" << anode->children[i]->name << endl;
            }
            else if(anode->nodeType==2)
            {
                cout << "(Duplicate) " <<anode->species <<"|" <<anode->name<<"\tleft\t" << anode->children[i]->species<<"|" << anode->children[i]->name << endl;
            }
            else
            {
                cout <<"Unrecognized node type " << anode->nodeType << endl;
            }
            showTree(anode->children[i]);
        }
    }
	return 0;
}

GeneTree*
GeneTreeManager::createTree(SpeciesDistManager::Species* duplicateSpecies,SpeciesDistManager::Species* node, MappedOrthogroup* mor, int copyno)
{
	//At each node make a GeneTree node
	//if(node->leftchild==NULL)
	if(node->children.empty())
	{
		GeneTree* aNode=new GeneTree;
		GeneMap* genemap=mor->getSpeciesHits(node->name.c_str());
		if(genemap==NULL)
		{
			return NULL;
		}
		map<string,map<string,STRINTMAP*>*>& geneSet=genemap->getGeneSet();
		map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();
        /*
		if(node==duplicateSpecies)
		{
			aNode->nodeType=2;
			aNode->name.append(node->name);
			aNode->name.append("_2");
			aNode->species.append(node->name);
			GeneTree* leftchild=new GeneTree;
			GeneTree* rightchild=new GeneTree;
			leftchild->nodeType=1;
			leftchild->name.append(gIter->first);
			gIter++;
			rightchild->name.append(gIter->first);
			rightchild->nodeType=1;
			aNode->leftchild=leftchild;
			aNode->rightchild=rightchild;
			leftchild->parent=aNode;
			rightchild->parent=aNode;
		}
		else
         */
		{
			aNode->nodeType=1;
			aNode->species.append(node->name);
			if(copyno==1)
			{
				aNode->name.append(gIter->first);
			}
			else
			{
				gIter++;
				aNode->name.append(gIter->first);
			}
		}
		
		return aNode;
	}
	GeneTree* geneTreeNode=new GeneTree;
	geneTreeNode->species.append(node->name);
    /*
	if(duplicateSpecies==node)
	{
		geneTreeNode->nodeType=2;
		string newName1(node->name);
		newName1.append("_1");
		GeneTree* copy1=new GeneTree;
		copy1->name.append(newName1);
		copy1->species.append(node->name);
		geneTreeNode->leftchild=copy1;
		GeneTree* newleft1=createTree(duplicateSpecies,node->leftchild,mor,1);
		GeneTree* newright1=createTree(duplicateSpecies,node->rightchild,mor,1);
		copy1->leftchild=newleft1;
		copy1->rightchild=newright1;
		if(newleft1!=NULL)
		{
			newleft1->parent=copy1;
		}
		if(newright1!=NULL)
		{
			newright1->parent=copy1;
		}

		string newName2(node->name);
		newName2.append("_2");
		GeneTree* copy2=new GeneTree;
		copy2->name.append(newName2);
		copy2->species.append(node->name);
		geneTreeNode->rightchild=copy2;
		GeneTree* newleft2=createTree(duplicateSpecies,node->leftchild,mor,2);
		GeneTree* newright2=createTree(duplicateSpecies,node->rightchild,mor,2);
		copy2->leftchild=newleft2;
		copy2->rightchild=newright2;
		if(newleft2!=NULL)
		{
			newleft2->parent=copy2;
		}	
		if(newright2!=NULL)
		{
			newright2->parent=copy2;
		}
		
		return geneTreeNode;
	    
    }*/
	//else
	{
		geneTreeNode->nodeType=1;
		char newName[256];
            
		GeneMap* genemap=mor->getSpeciesHits(node->name.c_str());
		if(genemap==NULL)
		{
            sprintf(newName,"%s",node->name.c_str()); //,copyno);
		    geneTreeNode->name.append(newName);
		
        }else{
		    map<string,map<string,STRINTMAP*>*>& geneSet=genemap->getGeneSet();
		    map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();
		    geneTreeNode->name.append(gIter->first);
        }
        /*
        GeneTree* newleft1=createTree(duplicateSpecies,node->leftchild,mor,1);
        GeneTree* newright1=createTree(duplicateSpecies,node->rightchild,mor,1);
        geneTreeNode->leftchild=newleft1;
        geneTreeNode->rightchild=newright1;
    
        
        if(newleft1!=NULL)
        {
            newleft1->parent=geneTreeNode;
        }
        if(newright1!=NULL)
        {
            newright1->parent=geneTreeNode;
        }*/
        
        for (int i=0; i<node->children.size(); i++) {
            GeneTree* newleft1=createTree(duplicateSpecies,node->children[i],mor,1);
            if(newleft1!=NULL)
            {
                geneTreeNode->children.push_back(newleft1);
                newleft1->parent=geneTreeNode;
            }
        }

	}
	return geneTreeNode;
}
/*
GeneTree*
GeneTreeManager::readTreeFromFile(const char* treeFName)
{
	NewickReader nreader;
	GeneTree* gtree=nreader.readTree(treeFName);
	
	return gtree;
}*/

/*
int
GeneTreeManager::pruneTree(GeneTree* node)
{
	//cout <<"Pruning " << node->species <<"|" <<node->name << endl;
	//If this is a leaf node and is a species which is not of interest delete it
	
    if((strstr(node->name.c_str(),"Internal")==NULL) && (strstr(node->name.c_str(),"Anc")==NULL))
	{
		map<string,SpeciesDistManager::Species*>& specSet=sdMgr->getAllSpecies();
		if(specSet.find(node->species)!=specSet.end())
		{
			return 0;
		}
		GeneTree* parent=node->parent;
		if(parent->leftchild==node)
		{
			parent->leftchild=NULL;
		}
		else if(parent->rightchild==node)
		{
			parent->rightchild=NULL;
		}
		delete node;
		return 0;
	}
	if(strcmp(node->name.c_str(),"Internal8")==0)
	{
	//	cout << "Stop here" << endl;
	}
	//If this is an intermediate node then call 
	int leftPrune=0;
	int rightPrune=0;
	if(node->leftchild!=NULL)
	{
		pruneTree(node->leftchild);
		if(node->leftchild==NULL)
		{
			leftPrune=1;
		}
	}	
	if(node->rightchild!=NULL)
	{
		pruneTree(node->rightchild);
		if(node->rightchild==NULL)
		{
			rightPrune=1;
		}
	}
	//Prune only if the species is not of interest
	if(strstr(node->species.c_str(),"Internal")!=NULL)
	{
		if(leftPrune==1 && rightPrune==1)
		{
	//		cout << node->species <<"|" << node->name << " was left and right child pruned " << endl;
			GeneTree* parent=node->parent;
			if(parent->leftchild==node)
			{
				parent->leftchild=NULL;
			}
			else if(parent->rightchild==node)
			{
				parent->rightchild=NULL;
			}
		}
		else if(leftPrune==1)
		{
	//		cout << node->species <<"|" << node->name << " was left child pruned " << endl;
			GeneTree* parent=node->parent;
			if(parent!=NULL)
			{
				if(parent->leftchild==node)
				{
					parent->leftchild=node->rightchild;
					node->rightchild->parent=parent;
				}
				else if(parent->rightchild==node)
				{
					parent->rightchild=node->rightchild;
					node->rightchild->parent=parent;
				}
			}
		}
		//If only the left branch was pruned then the rightchild of this node becomes the left child of the parent
		else if(rightPrune==1)
		{	
			//cout << node->species <<"|" << node->name << " was right child pruned " << endl;
			GeneTree* parent=node->parent;
			if(parent!=NULL)
			{
				if(parent->leftchild==node)
				{
					parent->leftchild=node->leftchild;
					node->leftchild->parent=parent;
				}
				else if(parent->rightchild==node)
				{
					parent->rightchild=node->leftchild;
					node->leftchild->parent=parent;
				}
			}
		}
		//If only the right branch was pruned then leftchild becomes the leftchild
		if(rightPrune==1 || leftPrune==1)
		{
			if(node->parent!=NULL)
			{
			//	cout <<"Deleting " << node->species <<"|" << node->name << endl;
				delete node;
			}
		}
	}
	return 0;
}

/***************************
int
GeneTreeManager::insertTree(GeneTree* node)
{
	if(node->leftchild!=NULL && node->rightchild!=NULL)
	{
		int leftlevel=sdMgr->getLevelFromRoot(node->leftchild->species.c_str());
		int rightlevel=sdMgr->getLevelFromRoot(node->rightchild->species.c_str());
		if(leftlevel>rightlevel)
		{
			insertTree(node->leftchild);
			if(node->rightchild!=NULL)
			{
				insertTree(node->rightchild);
			}
		}
		else
		{
			insertTree(node->rightchild);
			if(node->leftchild!=NULL)
			{
				insertTree(node->leftchild);
			}
		}
		if(strcmp(node->leftchild->species.c_str(),node->rightchild->parent->species.c_str())==0)
		{
			insertTree(node->rightchild);
			insertTree(node->leftchild);
		}
		else if(strcmp(node->rightchild->species.c_str(),node->leftchild->parent->species.c_str())==0)
		{
			insertTree(node->leftchild);
			insertTree(node->rightchild);
		}
		else
		{
			insertTree(node->leftchild);
			insertTree(node->rightchild);
		}
	}
	else if(node->leftchild!=NULL)
	{
		insertTree(node->leftchild);
	}
	else if(node->rightchild!=NULL)
	{
		insertTree(node->rightchild);
	}
	//Need to make sure that two siblings of a node are sibling species in the species tree
	//if(strstr(node->name.c_str(),"Internal")!=NULL)
	//{
		int left=0;
		int right=0;
		//Check if the sibling of this node matches the species tree
		GeneTree* parent=node->parent;
		if(parent==NULL)
		{
			return 0;
		}
		GeneTree* sibling=NULL;
		if(parent->leftchild==node)
		{
			sibling=parent->rightchild;
			left=1;
		}
		else if(parent->rightchild==node)
		{
			sibling=parent->leftchild;
			right=1;
		}
		SpeciesDistManager::Species* speciesnode=sdMgr->getSpecies(node->species);
		if(speciesnode==NULL)
		{
			return 0;
		}
		if(sibling==NULL)
		{
			//Just make sure that this node is the correct (left or right) child of its parent as specified by the species tree
			SpeciesDistManager::Species* speciesparent=speciesnode->parent;
			if(speciesparent==NULL)
			{
				return 0;
			}
			if(speciesparent->leftchild==speciesnode && parent->rightchild==node)
			{
				parent->leftchild=node;
				parent->rightchild=NULL;
			}
			else if(speciesparent->rightchild==speciesnode && parent->leftchild==node)
			{
				parent->rightchild=node;
				parent->leftchild=NULL;
			}
			return 0;
		}
		//If this happens then sibling and this node are duplicates
		if(strcmp(sibling->species.c_str(),node->species.c_str())==0)
		{
			return 0;	
		}
		SpeciesDistManager::Species* speciesparent=speciesnode->parent;
		if(speciesparent==NULL)
		{
			return 0;
		}
		SpeciesDistManager::Species* speciessibling=NULL;
		int leftspecies=0;
		int rightspecies=0;
		if(speciesparent->leftchild==speciesnode)
		{
			speciessibling=speciesparent->rightchild;
			leftspecies=1;
		}
		else if(speciesparent->rightchild==speciesnode)
		{
			speciessibling=speciesparent->leftchild;
			rightspecies=1;
		}
		if(strcmp(sibling->species.c_str(),speciessibling->name.c_str())==0)
		{
			//make sure the parent name in the gene tree is the same as the species tree
			if(strcmp(parent->species.c_str(),speciesparent->name.c_str())!=0)
			{
				parent->species.clear();
				parent->species.append(speciesparent->name);
				parent->name.clear();
				parent->name.append(speciesparent->name);
			}
			return 0;
		}
		else if(strcmp(sibling->species.c_str(),parent->name.c_str())!=0)
		{
			//Not sure if it is this node that needs to be modified or the sibling 
			if(strcmp(speciesparent->name.c_str(),parent->species.c_str())==0)
			{
				return 0;
			}
		}
		//cout <<"Inserting parent " << speciesparent->name << " to " << node->name << endl;
		GeneTree* newparent=new GeneTree;
		newparent->species.append(speciesparent->name.c_str());
		newparent->name.append(speciesparent->name.c_str());
		newparent->leftchild=NULL;
		newparent->rightchild=NULL;
		if(leftspecies==1)
                {
                        newparent->leftchild=node;
                }
                else if(rightspecies==1)
                {
                        newparent->rightchild=node;
                }
		int parentlevel=sdMgr->getLevelFromRoot(parent->species.c_str());
		int newparentlevel=sdMgr->getLevelFromRoot(newparent->species.c_str());
		if(parentlevel <= newparentlevel)
		{
			if(left==1)
			{
				//Here we need to see where the newparent is added.
				parent->leftchild=newparent;
			}
			else if(right==1)
			{
				parent->rightchild=newparent;
			}
			newparent->parent=parent;
		}
		else 
		{
			if(parent->parent!=NULL)
			{
				if(parent->parent->leftchild==parent)
				{
					parent->parent->leftchild=newparent;
				}
				else if(parent->parent->rightchild==parent)
				{
					parent->parent->rightchild=newparent;
				}
			}
			if(parentlevel > newparentlevel)
			{
				if(newparent->leftchild==NULL)
				{
					newparent->leftchild=parent;
				}
				else if(newparent->rightchild==NULL)
				{	
					newparent->rightchild=parent;
				}
				newparent->parent=parent->parent;
				parent->parent=newparent;
				if(parent->leftchild==node)
				{
					parent->leftchild=NULL;
				}
				else if(parent->rightchild==node)
				{
					parent->rightchild=NULL;
				}
			}
			else if(parentlevel==newparentlevel)
			{
				if(parent->leftchild!=node)
				{
					if(newparent->leftchild==NULL)
					{
						newparent->leftchild=parent->leftchild;
					}
					else if(newparent->rightchild==NULL)
					{	
						newparent->rightchild=parent->leftchild;
					}
				}
				else if(parent->rightchild!=node)
				{
					if(newparent->leftchild==NULL)
					{
						newparent->leftchild==parent->rightchild;
					}
					else if(newparent->rightchild==NULL)
					{
						newparent->rightchild=parent->rightchild;
					}
				}
			}
		}
		insertTree(newparent);
	//}
	return 0;
}


const char*
GeneTreeManager::getSpeciesName(GeneTree* node)
{
	const char* queryLeft=NULL;
	const char* queryRight=NULL;
	if(node->leftchild!=NULL)	
	{
		queryLeft=getSpeciesName(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		queryRight=getSpeciesName(node->rightchild);
	}
	if(queryLeft==NULL && queryRight==NULL)
	{
		SpeciesDistManager::Species* speciesnode=sdMgr->getSpecies(node->species);
		if(speciesnode==NULL)
		{
			cout <<"No species by name " << node->species << endl;
			return 0;
		}
		SpeciesDistManager::Species* speciesparent=speciesnode->parent;
		const char* parentname=speciesparent->name.c_str();
		return parentname;
	}
	if(queryLeft==NULL)
	{
		node->species.clear();
		node->species.append(queryRight);
		node->name.clear();
		node->name.append(queryRight);
	}
	else if(queryRight==NULL)
	{
		node->species.clear();
		node->species.append(queryLeft);
		node->name.clear();
		node->name.append(queryLeft);
	}
	else
	{
		int leftLevel=sdMgr->getLevelFromRoot(queryLeft);
		int rightLevel=sdMgr->getLevelFromRoot(queryRight);
		node->species.clear();
		node->name.clear();
		//Check is this is a duplicate node. This happens if both children have the same names
		//OR if one of child's vote for ancestor is the same as the other child's name 
		if(strcmp(queryLeft,node->rightchild->name.c_str())==0) 
		{
			node->species.append(queryLeft);
		}
		else if(strcmp(node->leftchild->name.c_str(),queryRight)==0)
		{
			node->species.append(queryRight);
		}
		//if(leftLevel>rightLevel)
		else
		{
			if(leftLevel<rightLevel)
			{
				node->species.append(queryLeft);
			}	
			else if(leftLevel>rightLevel)
			{
				node->species.append(queryRight);
			}
			else if(leftLevel==rightLevel)
			{
				if(strcmp(node->leftchild->species.c_str(),node->rightchild->species.c_str())==0)
				{
					node->species.append(node->leftchild->species);
				}
				else
				{
					node->species.append(queryLeft);
				}
			}
		}
		node->name.append(node->species.c_str());
	}
	SpeciesDistManager::Species* speciesnode=sdMgr->getSpecies(node->species);
	if(speciesnode==NULL)
	{
		cout <<"No species by name " << node->species << endl;
		return 0;
	}
	SpeciesDistManager::Species* speciesparent=speciesnode->parent;
	if(speciesparent==NULL)
	{
		return NULL;
	}
	const char* parentname=speciesparent->name.c_str();
	return parentname;
}


int
GeneTreeManager::getSpeciesName(GeneTree* node,map<string,int>& ancestors) 
{
	map<string,int> queryLeft;
	map<string,int> queryRight;
	if(node->leftchild!=NULL)	
	{
		getSpeciesName(node->leftchild,queryLeft);
	}
	if(node->rightchild!=NULL)
	{
		getSpeciesName(node->rightchild,queryRight);
	}
	if(queryLeft.size()==0 && queryRight.size()==0)
	{
		SpeciesDistManager::Species* speciesnode=sdMgr->getSpecies(node->species);
		if(speciesnode==NULL)
		{
			//cout <<"No species by name " << node->species << endl;
			return 0;
		}
		//speciesnode=speciesnode->parent;
		//while(speciesnode!=NULL)
		//{
		//	ancestors[speciesnode->name]=sdMgr->getLevelFromRoot(speciesnode->name.c_str());
		//	speciesnode=speciesnode->parent;
		//}
		//queryLeft.clear();
		//queryRight.clear();
		//return 0;
	}
	else if(queryLeft.size()==0)
	{
		int lowestAncestor=-1;
		for(map<string,int>::iterator aIter=queryRight.begin();aIter!=queryRight.end();aIter++)
		{
			if(aIter->second>lowestAncestor)
			{
				lowestAncestor=aIter->second;
				node->species.clear();
				node->species.append(aIter->first.c_str());
			}
		}
		node->name.clear();
		node->name.append(node->species.c_str());
	}
	else if(queryRight.size()==0)
	{
		int lowestAncestor=-1;
		for(map<string,int>::iterator aIter=queryLeft.begin();aIter!=queryLeft.end();aIter++)
		{
			if(aIter->second>lowestAncestor)
			{
				lowestAncestor=aIter->second;
				node->species.clear();
				node->species.append(aIter->first.c_str());
			}
		}
		node->name.clear();
		node->name.append(node->species.c_str());
	}
	
	else
	{
		int deepestCommonLevel=-1;
		string ancName;
		for(map<string,int>::iterator aIter=queryLeft.begin();aIter!=queryLeft.end();aIter++)
		{
			if(queryRight.find(aIter->first)==queryRight.end())
			{
				continue;
			}
			if(aIter->second>=deepestCommonLevel)
			{
				ancName.clear();
				deepestCommonLevel=aIter->second;
				ancName.append(aIter->first);
			}
		}
		
		node->species.clear();
		node->name.clear();
		if(deepestCommonLevel==-1)
		{
			cout <<"No common ancestors for " << node->leftchild->name << "\t" << node->rightchild->name << endl;
			exit(0);
		}
		//Check is this is a duplicate node. This happens if both children have the same names
		//OR if one of child's vote for ancestor is the same as the other child's name 
		node->species.append(ancName.c_str());
		node->name.append(node->species.c_str());
	}
	queryLeft.clear();
	queryRight.clear();
	SpeciesDistManager::Species* speciesnode=sdMgr->getSpecies(node->species);
	if(speciesnode==NULL)
	{
		//cout <<"No species by name " << node->species << endl;
		return 0;
	}
	if(node->parent!=NULL)
	{
		//add self as a possible ancestor only if we know that there is a possible duplicate
		if(node->parent->leftchild!=NULL && node->parent->rightchild!=NULL)
		{
			ancestors[speciesnode->name]=sdMgr->getLevelFromRoot(speciesnode->name.c_str());
		}
	}
	speciesnode=speciesnode->parent;
	while(speciesnode!=NULL)
	{
		ancestors[speciesnode->name]=sdMgr->getLevelFromRoot(speciesnode->name.c_str());
		speciesnode=speciesnode->parent;
	}

	return 0;
}


int
GeneTreeManager::setNodeType(GeneTree* node)
{
	if(node->leftchild!=NULL && node->rightchild!=NULL)
	{
		if(strcmp(node->leftchild->species.c_str(),node->rightchild->species.c_str())==0)
		{
			node->nodeType=2;
		}
		else
		{
			node->nodeType=1;
		}
	}
	if(node->leftchild!=NULL)
	{
		setNodeType(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		setNodeType(node->rightchild);
	}
	return 0;
}
*/
GeneTree* 
GeneTreeManager::getGeneTree(int ogid)
{
	GeneTree* gtree=NULL;
	if(genetreeSet.find(ogid)==genetreeSet.end())
	{
		return gtree;
	}
	gtree=genetreeSet[ogid];
	return gtree;
}	


GeneTree* 
GeneTreeManager::getGeneTree(MappedOrthogroup* mor)
{
	GeneTree* gtree=NULL;
	int ogid=mor->getID();
	if(genetreeSet.find(ogid)!=genetreeSet.end())
	{
		gtree=genetreeSet[ogid];
	}
	else
	{
		gtree=generateTree(mor);
		//gtree=generateTreeFromFile(mor);
	}
	return gtree;
}	
/*
int
GeneTreeManager::addKpolGene(GeneTree* gtree,string& geneName)
{
	string kpolname("Kpol");
	SpeciesDistManager::Species* kpol=sdMgr->getSpecies(kpolname);
	GeneTree* kpolanc=isKpolAnc(gtree,kpol->parent->name);
	if(kpolanc==NULL)
	{
		cout <<"No Kpol node in tree! " << endl;
		exit(0);
	}
	GeneTree* kpolnode=new GeneTree;
	kpolnode->species.append("Kpol");
	kpolnode->name.append(geneName);
	kpolnode->parent=kpolanc;
	if(kpolanc->rightchild!=NULL)
	{
		cout <<kpolanc->name << " already has a child!" << endl;
		exit(0);
	}
	kpolanc->rightchild=kpolnode;
	return 0;
}

GeneTree*
GeneTreeManager::isKpolAnc(GeneTree* gtree,string& parentName)
{
	
    GeneTree* hit=NULL;
	if(strcmp(gtree->species.c_str(),parentName.c_str())==0)
	{
		hit=gtree;
	}
	else
	{
		if(gtree->leftchild!=NULL)
		{
			hit=isKpolAnc(gtree->leftchild,parentName);
			if(hit!=NULL)
			{
				return hit;
			}
		}
		if(gtree->rightchild!=NULL)
		{
			hit=isKpolAnc(gtree->rightchild,parentName);
			if(hit!=NULL)
			{
				return hit;
			}
		}
	}
     
	return hit;
}
***************************/

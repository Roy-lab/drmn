
/*
CMINT: An algorithm to cluster functional genomics data from multiple cell types
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
#include <iostream>

#include <string.h>
#include "GeneTree.H"
#include "Matrix.H"
#include "Gamma.H"

Gamma::Gamma()
{
	dupAncName.clear();	
	nodeID=0;
}

Gamma::~Gamma()
{
}


int 
Gamma::setGeneClusterID(string& geneName,string& specName,int clusterID)
{
	if(geneMap.find(geneName)==geneMap.end())
	{
		cout <<"No gene " << geneName <<  " in species " << specName << endl;
		// no, everything is not fine. 
		exit(0);
	}
	Gamma::Node* gene=geneMap[geneName];
	gene->hasData=true;
	Matrix* gamma=gene->gamma;
/*	for(int r=0;r<gamma->getRowCnt();r++)
	{
		gamma->setValue(1.0,r,clusterID);
	}*/
	gamma->setValue(1.0,0,clusterID);
	return 0;
}

//Here we will assume that there is a single duplication even at most. In the future, we want
//to initialize this directly from the gene tree.
int 
Gamma::initUsingGeneTree(GeneTree* gtree)
{
	//initSubtree(gtree);
	initSubtree(gtree,gtree->species);
	return 0;
}

Gamma::Node*
Gamma::initSubtree(GeneTree* gtreenode)
{
	Gamma::Node* node=new Gamma::Node;
	node->hasData=false;
	node->dataLL=NULL;
	node->parent=NULL;
	node->excludeFlag=gtreenode->exclude;
	node->species.append(gtreenode->species);
	char nodeName[1024];
	if(strstr(node->species.c_str(),"Anc")!=NULL)
	{
		sprintf(nodeName,"%s_%d",gtreenode->name.c_str(),nodeID);
		nodeID++;
	}
	else
	{
		sprintf(nodeName,"%s",gtreenode->name.c_str());
	}
	node->name.append(nodeName);
	//if(gtreenode->nodeType==1)
	//{
		if(gtreenode->parent==NULL)
		{
			root=node;
			node->gamma=new Matrix(1,maxClusterCnt);
		}
		else
		{
			node->gamma=new Matrix(maxClusterCnt,maxClusterCnt);
		}
		node->gamma->setAllValues(0);
		node->normTerm=NULL;
	//}
	geneMap[node->name]=node;
	node->nodeType=gtreenode->nodeType;
	if(node->nodeType==2)
	{
		if(dupAncName.length()>0)
		{
			cout <<"Gamma already has a dup ancesor set! " << dupAncName << endl;
		}
		else
		{
			dupAncName.append(node->species);
		}
	}
    for (int i=0;i<gtreenode->children.size();i++) {
        Gamma::Node* child=initSubtree(gtreenode->children[i]);
        node->children.push_back(child);
	child->parent=node;
    }
    
	/*if(gtreenode->leftchild!=NULL)
	{
		Gamma::Node* leftchild=initSubtree(gtreenode->leftchild);
		//Gamma::Node* leftchild=geneMap[gtreenode->leftchild->name];
		node->leftchild=leftchild;
		leftchild->parent=node;
	}
	else
	{
		node->leftchild=NULL;
	}
	if(gtreenode->rightchild!=NULL)
	{
		Gamma::Node* rightchild=initSubtree(gtreenode->rightchild);
		//Gamma::Node* rightchild=geneMap[gtreenode->rightchild->name];
		node->rightchild=rightchild;
		rightchild->parent=node;
	}
	else
	{
		node->rightchild=NULL;
	}
     */
	return node;
}


Gamma::Node*
Gamma::initSubtree(GeneTree* gtreenode,string& rootname)
{
	Gamma::Node* node=new Gamma::Node;
	node->hasData=false;
	node->dataLL=NULL;
	node->parent=NULL;
	node->excludeFlag=gtreenode->exclude;
	node->species.append(gtreenode->species);
	char nodeName[1024];
	if(strstr(node->species.c_str(),"Anc")!=NULL)
	{
		sprintf(nodeName,"%s_%d",gtreenode->name.c_str(),nodeID);
		nodeID++;
	}
	else
	{
		sprintf(nodeName,"%s",gtreenode->name.c_str());
	}
	node->name.append(nodeName);
	//if(gtreenode->nodeType==1)
	//{
		//if(gtreenode->parent==NULL)
		if(strcmp(gtreenode->species.c_str(),rootname.c_str())==0)
		{
			root=node;
			//node->parent=NULL;
			node->gamma=new Matrix(1,maxClusterCnt);
		}
		else
		{
			node->gamma=new Matrix(maxClusterCnt,maxClusterCnt);
		}
		node->gamma->setAllValues(0);
		node->normTerm=NULL;
	//}
	geneMap[node->name]=node;
	node->nodeType=gtreenode->nodeType;
	if(node->nodeType==2)
	{
		if(dupAncName.length()>0)
		{
			cout <<"Gamma already has a dup ancesor set! " << dupAncName << endl;
		}
		else
		{
			dupAncName.append(node->species);
		}
	}
    
    for (int i=0;i<gtreenode->children.size();i++) {
        Gamma::Node* child=initSubtree(gtreenode->children[i],rootname);
        node->children.push_back(child);
        if(strcmp(child->species.c_str(),rootname.c_str())!=0)
		{
			child->parent=node;
		}

    }

    
    /*
	if(gtreenode->leftchild!=NULL)
	{
		Gamma::Node* leftchild=initSubtree(gtreenode->leftchild,rootname);
		//Gamma::Node* leftchild=geneMap[gtreenode->leftchild->name];
		node->leftchild=leftchild;
		if(strcmp(leftchild->species.c_str(),rootname.c_str())!=0)
		{
			leftchild->parent=node;
		}
	}
	else
	{
		node->leftchild=NULL;
	}
	if(gtreenode->rightchild!=NULL)
	{
		Gamma::Node* rightchild=initSubtree(gtreenode->rightchild,rootname);
		//Gamma::Node* rightchild=geneMap[gtreenode->rightchild->name];
		node->rightchild=rightchild;
		if(strcmp(rightchild->species.c_str(),rootname.c_str())!=0)
		{
			rightchild->parent=node;
		}
	}
	else
	{
		node->rightchild=NULL;
	}
     */
	return node;
}

Matrix* 
Gamma::getGamma(string& geneName,string& speciesName)
{
	if(geneMap.find(geneName)==geneMap.end())
	{
		return NULL;
	}
	Gamma::Node* node=geneMap[geneName];
	return node->gamma;
}

Matrix* 
Gamma::getNormTerm(string& geneName,string& speciesName)
{
	if(geneMap.find(geneName)==geneMap.end())
	{
		return NULL;
	}
	Gamma::Node* node=geneMap[geneName];
	return node->normTerm;
}


Gamma::Node* 
Gamma::getGeneNode(string& geneName,string& speciesName)
{
	if(geneMap.find(geneName)==geneMap.end())
	{
		return NULL;
	}
	Gamma::Node* node=geneMap[geneName];
	return node;
}

Gamma::Node* 
Gamma::getRoot()
{
	return root;
}


int 
Gamma::showTree()
{
	showTree(root);
	return 0;
}


map<int,double>*
Gamma::getDataLL(string& geneName)
{
	if(datall.find(geneName)==datall.end())
	{
		Gamma::Node* node=geneMap[geneName];
		if(node->hasData)
		{
			map<int,double>* pvals=new map<int,double>;
			datall[geneName]=pvals;
			node->dataLL=pvals;
			return pvals;
		}
		cout <<"No gene with name " << geneName << endl;
		return NULL;
	}
	return datall[geneName];
}

string&
Gamma::getDupAncestor()
{
	return dupAncName;
}

int
Gamma::showTree(Gamma::Node* node)
{
	//if(node->leftchild==NULL && node->rightchild==NULL)
    if(node->children.empty())
	{
		cout << node->name <<"-> null"<< endl;
        return 0;
	}
    
    for (int i=0;i<node->children.size();i++) {
        cout << node->name <<"-> "<< node->children[i]->name << endl;
        showTree(node->children[i]);
    }
    
	/*if(node->leftchild!=NULL)
	{
		cout << node->name <<"->left "<< node->leftchild->name << endl;
	}
	else
	{
		cout << node->name <<"->left null"<< endl;
	}
	if(node->rightchild!=NULL)
	{
		cout << node->name <<"->right "<< node->rightchild->name << endl;
	}
	else
	{
		cout << node->name <<"->right null"<< endl;
	}
	if(node->leftchild!=NULL)
	{
		showTree(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		showTree(node->rightchild);
	}
     */

	return 0;
}

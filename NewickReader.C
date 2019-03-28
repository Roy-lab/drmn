
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
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

#include "GeneTree.H"
#include "NewickReader.H"

NewickReader::NewickReader()
{
	root=NULL;
}

NewickReader::~NewickReader()
{
}

GeneTree*
NewickReader::readTree(const char* aFName)
{
	ifstream inFile(aFName);
	string strbuff;
	getline(inFile,strbuff);
	int len=strbuff.length();
	char* buffer=new char[len+1];
	strcpy(buffer,strbuff.c_str());
	inFile.close();
	GeneTree* gtree=parseNewickFormat(buffer);
	delete[] buffer;
	return gtree;	
}

GeneTree*
NewickReader::parseNewickFormat(char* newickStr)
{
	/*int depthCnt=0;
	int index=0;
	char nodeName[256];
	char distance[256];
	int nodeindex=0;
	int distindex=0;
	GeneTree* currNode=root;
	bool popName=true;
	while(newickStr[index]!='\0')
	{
		char c=newickStr[index];
		switch (c)
		{
			case '(':
			{
				GeneTree* newNode=new GeneTree;
				char aname[25];
				sprintf(aname,"Internal%d",depthCnt);
				newNode->name.append(aname);
				newNode->species.append(aname);
				if(currNode!=NULL)
				{
					//cout<<"Found child for "<< currNode->getName() << endl;
                    
					if((currNode->leftchild==NULL) && (strstr(currNode->name.c_str(),"Internal")!=NULL))
					{
						newNode->parent=currNode;
						currNode->leftchild=newNode;
					}
					else
					{
						GeneTree* parent=currNode->parent;
						if(parent==NULL)
						{
							cout <<"Bad format! No parent for " << currNode->name << endl;
							exit(0);
						}
						if(parent->rightchild!=NULL)
						{
							cout <<"Bad format! Right child of parent " << parent->name << " already set to " << parent->rightchild->name << endl;
							exit(0);
						}
						parent->rightchild=newNode;
						newNode->parent=parent;
					}
				}
				else
				{
					root=newNode;	
				}	
				currNode=newNode;
				depthCnt++;
				break;
			}
			case ')':
			{
				if(nodeindex>0)
				{
					nodeName[nodeindex]='\0';
					if(strcmp(nodeName,"Agos")==0)
					{
						cout <<"Stop here" << endl;
					}
					GeneTree* newChild=new GeneTree;
					char* pos=strchr(nodeName,'|');
					if(pos==NULL)
					{
						cout <<"No species name found in " << nodeName << endl;
						exit(0);
					}
					*pos='\0';
					newChild->name.append(pos+1);
					newChild->species.append(nodeName);
					newChild->parent=currNode->parent;
					GeneTree* parent=currNode->parent;
					parent->rightchild=newChild;
					nodeindex=0;
					distindex=0;
				}
				//cout <<"Moving up from " << currNode->getName() << " to " << currNode->getParent()->getName()<< endl;
				//Need to move up the tree
				currNode=currNode->parent;
				break;
			}
			case ',':
			{
				popName=true;
				if(nodeindex>0)
				{
					nodeName[nodeindex]='\0';
					GeneTree* newChild=new GeneTree;
					char* pos=strchr(nodeName,'|');
					if(pos==NULL)
					{
						cout <<"Bad format. No species name in " << nodeName << endl;
						exit(0);
					}
					*pos='\0';
					newChild->name.append(pos+1);
					newChild->species.append(nodeName);
					newChild->parent=currNode;
					currNode->leftchild=newChild;
					currNode=newChild;
				}
				
				distance[distindex]='\0';
				nodeindex=0;
				distindex=0;
				break;
			}
			case ';':
			{
				if(nodeindex>0)
				{
					nodeName[nodeindex]=0;
				}
				break;
			}
			case ':':
			{
				popName=false;	
				break;
			}
			default:
			{
				if(popName)
				{
					nodeName[nodeindex]=c;
					nodeindex++;
				}
				else
				{
					distance[distindex]=c;
					distindex++;
				}
			}
		}
		index++;
	}
	return root;*/return NULL;
}

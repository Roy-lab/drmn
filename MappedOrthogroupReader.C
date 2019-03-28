
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
#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"

MappedOrthogroupReader::MappedOrthogroupReader()
{
}

MappedOrthogroupReader::~MappedOrthogroupReader()
{
}

//This function reads the orthogroup mapping produced in any one of the PROPER, FILLED and UNFILLED formats.
//We are essentialy interested in the first two columns
int 
MappedOrthogroupReader::readFile(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string buffstr;
	int lineCnt=0;
	int oldogid=-2;
	MappedOrthogroup* ogrp=NULL;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(lineCnt==0)
		{
			lineCnt++;
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(strstr(buffer,"orf19.993")!=NULL)
		{
			cout <<"Stop here" << endl;
		}
		//Get the orthogroup and the genes
		int ogid=-1;
		string geneSet;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				// here we're splitting tok by _
				char* pos=strchr(tok,'_');
				if(pos==NULL)
				{
					cout <<"Bad OG name " << tok<<  endl;
					exit(0);
				}
				*pos='\0';

				// now we need to get just the numeric part of the token
				//ogid=atoi(tok+2); // assumes format is "OG1_0"
				string oString(tok);
				string r1="OG";
				string r2="ID";
				size_t tpos=oString.find(r1);
				if (tpos != string::npos) oString.erase(tpos,r1.length());
				tpos=oString.find(r2);
				if (tpos != string::npos) oString.erase(tpos,r2.length());

				// now convert to integer
				ogid=atoi(oString.c_str());

				if(oldogid!=ogid)
				{
					ogrp=new MappedOrthogroup;
					ogrp->setID(ogid);
					orthogroupSet[ogrp->getID()]=ogrp;
					oldogid=ogid;
				}
				ogrp->incrCnt();
				
			}
			else if(tokCnt==1)
			{
				addMembers(tok,ogrp);
			}
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
		lineCnt++;
	}
	inFile.close();
	generateGeneOrthoMap();	
	return 0;
}

/*
 * Returns number of orthogroups
 */ 
int 
MappedOrthogroupReader::size() 
{
	return orthogroupSet.size();
}

/*
 * Reads the species order file
 */
int
MappedOrthogroupReader::readSpeciesMapping(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string speciesName(buffer);
		speciesIDNameMap[lineCnt]=speciesName;
		lineCnt++;
	}
	inFile.close();

	// DC add -- problem reading species order file. no entries.
	if (lineCnt==0)
	{
		return 1;
	}

	return 0;
}

/*
 * Populates a vector with the species in the order we read them
 */
int
MappedOrthogroupReader::getSpeciesOrder(vector<string>& order)
{
	for (map<int,string>::iterator sIter=speciesIDNameMap.begin(); sIter!=speciesIDNameMap.end(); sIter++)
	{
		order.push_back(sIter->second);
	}

	return 0;
}


int 
MappedOrthogroupReader::getMappedOrthogroupID(const char* geneName, const char* species)
{
	string speciesKey(species);
	if(geneOrthoMap.find(speciesKey)==geneOrthoMap.end())
	{
		return -1;
	}
	string geneKey(geneName);
	STRINTMAP* geneSet=geneOrthoMap[speciesKey]; 
	if(geneSet->find(geneKey)==geneSet->end())
	{
		return -1;
	}
	int orthoID=(*geneSet)[geneKey];
	return orthoID;
}


map<int,MappedOrthogroup*>& 
MappedOrthogroupReader::getMappedOrthogroups()
{
	return orthogroupSet;
}

MappedOrthogroup*
MappedOrthogroupReader::getMappedOrthogroup(const char* geneName, const char* species)
{
	int orthoID=getMappedOrthogroupID(geneName,species);
	MappedOrthogroup* og=orthogroupSet[orthoID];
	return og;
}


STRINTMAP* 
MappedOrthogroupReader::getOrtholog(const char* srcSpecName, const char* geneName, const char* targetSpecName)
{
	int ogid=getMappedOrthogroupID(geneName,srcSpecName);
	if(ogid==-1)
	{
		return NULL;
	}
	MappedOrthogroup* mgrp=orthogroupSet[ogid];
	STRINTMAP* orthohits=mgrp->getSpeciesHitsForGene(srcSpecName,targetSpecName,geneName);
	return orthohits;
}

int
MappedOrthogroupReader::addMembers(char* abuffer, MappedOrthogroup* ogrp)
{
	char* tok=abuffer;
	int tokCnt=0;
	map<string,string> specGeneMap;
	while(tok!=NULL)
	{
		char* end=strchr(tok,',');
		if(end!=NULL)
		{
			*end='\0';
		}
		if(speciesIDNameMap.find(tokCnt)==speciesIDNameMap.end())
		{
			cout <<"No species with id " << tokCnt << endl;
			exit(0);
		}
		string speciesName=speciesIDNameMap[tokCnt];
		string geneName(tok);
		if(strcmp(geneName.c_str(),"NONE")!=0)
		{
			specGeneMap[speciesName]=geneName;
		}
		tokCnt++;
		if(end!=NULL)
		{
			tok=end+1;
		}	
		else
		{
			tok=end;
		}
	}
	ogrp->setMembers(specGeneMap);
	specGeneMap.clear();
	return 0;
}

int 
MappedOrthogroupReader::generateGeneOrthoMap()
{
	for(map<int,MappedOrthogroup*>::iterator oIter=orthogroupSet.begin();oIter!=orthogroupSet.end();oIter++)
	{
		MappedOrthogroup* og=oIter->second;
		map<string,GeneMap*>& members=og->getOrthoMembers();
		for(map<string,GeneMap*>::iterator sIter=members.begin();sIter!=members.end();sIter++)
		{
			map<string,map<string,STRINTMAP*>*>& genes=sIter->second->getGeneSet();
			STRINTMAP* allgenesForSpecies=NULL;
			if(geneOrthoMap.find(sIter->first)==geneOrthoMap.end())
			{
				allgenesForSpecies=new STRINTMAP;
				geneOrthoMap[sIter->first]=allgenesForSpecies;
			}
			else
			{
				allgenesForSpecies=geneOrthoMap[sIter->first];
			}
			for(map<string,map<string,STRINTMAP*>*>::iterator gIter=genes.begin();gIter!=genes.end();gIter++)
			{
				if(strcmp(gIter->first.c_str(),"orf19.993")==0)
				{
					cout <<"Found " << gIter->first.c_str() << endl;
				}
				(*allgenesForSpecies)[gIter->first]=oIter->first; 
			}
		}
	}
	return 0;
}

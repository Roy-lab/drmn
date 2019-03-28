
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
#include "GeneMap.H"
#include "MappedOrthogroup.H"


MappedOrthogroup::MappedOrthogroup()
{
	cnts=0;
}

MappedOrthogroup::~MappedOrthogroup()
{
}

int 
MappedOrthogroup::setMembers(map<string,string>& speciesGeneMap)
{
	map<string,string>* geneSet=new map<string,string>;
	int setid=geneSets.size();
	geneSets[setid]=geneSet;
	for(map<string,string>::iterator sIter=speciesGeneMap.begin();sIter!=speciesGeneMap.end();sIter++)
	{
		(*geneSet)[sIter->first]=sIter->second;
		GeneMap* speciesHit_s=NULL;
		if(orthoMembers.find(sIter->first)==orthoMembers.end())
		{
			speciesHit_s=new GeneMap;
			orthoMembers[sIter->first]=speciesHit_s;
		}
		else
		{
			speciesHit_s=orthoMembers[sIter->first];
		}
		if(speciesGeneMap.size()==1)
		{
			speciesHit_s->addPair(sIter->second,"","");
		}
		
		map<string,string>::iterator tIter=sIter;
		tIter++;
		for(;tIter!=speciesGeneMap.end();tIter++)
		{
			GeneMap* speciesHit_t=NULL;
			speciesHit_s->addPair(sIter->second,tIter->first,tIter->second);
			if(orthoMembers.find(tIter->first)==orthoMembers.end())
			{
				speciesHit_t=new GeneMap;
				orthoMembers[tIter->first]=speciesHit_t;
			}	
			else
			{
				speciesHit_t=orthoMembers[tIter->first];
			}
			speciesHit_t->addPair(tIter->second,sIter->first,sIter->second);
		}
	}
	return 0;
}

	
int 
MappedOrthogroup::setID(int aid)
{
	oid=aid;
	return 0;
}

int 
MappedOrthogroup::getID()
{
	return oid;
}


map<string,GeneMap*>&
MappedOrthogroup::getOrthoMembers()
{
	return orthoMembers;
}

int 
MappedOrthogroup::incrCnt()
{
	cnts=cnts+1;
	return 0;
}

int 
MappedOrthogroup::getCnt()
{
	return cnts;
}

GeneMap* 
MappedOrthogroup::getSpeciesHits(const char* specName)
{
	string key(specName);
	if(orthoMembers.find(key)==orthoMembers.end())
	{
		return NULL;
	}
	return orthoMembers[key];
}

STRINTMAP* 
MappedOrthogroup::getSpeciesHitsForGene(const char* srcSpecName, const char* targetSpecName, const char* geneName)
{
	string key(srcSpecName);
	if(orthoMembers.find(key)==orthoMembers.end())
	{
		cout <<"No hits for species  " << srcSpecName << endl;
		return NULL;
	}
	GeneMap* geneMap=orthoMembers[key];
	STRINTMAP* hits=geneMap->getHits(geneName,targetSpecName);
	return hits;
}
map<int,map<string,string>*>& 
MappedOrthogroup::getGeneSets()
{
	return geneSets;
}


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

GeneMap::GeneMap()
{
}

GeneMap::~GeneMap()
{
}

int 
GeneMap::addPair(const string& srcGene, const string& targetSpecies, const string& targetGene)
{
	map<string,STRINTMAP*>* hitsForGene=NULL;
	if(geneSet.find(srcGene)==geneSet.end())
	{
		hitsForGene=new map<string,STRINTMAP*>;
		geneSet[srcGene]=hitsForGene;
	}
	else
	{
		hitsForGene=geneSet[srcGene];
	}
	if(targetSpecies.length()==0)
	{
		return 0;
	}
	STRINTMAP* hitsInSpecies=NULL;
	if(hitsForGene->find(targetSpecies)==hitsForGene->end())
	{
		hitsInSpecies=new STRINTMAP;
		(*hitsForGene)[targetSpecies]=hitsInSpecies;	
	}
	else
	{
		hitsInSpecies=(*hitsForGene)[targetSpecies];	
	}
	(*hitsInSpecies)[targetGene]=0;
	return 0;
}

STRINTMAP* 
GeneMap::getHits(const char* srcGene, const char* targetSpecies)
{
	string sKey(srcGene);
	if(geneSet.find(sKey)==geneSet.end())
	{
		cout <<"Weird! No  gene " << sKey << " in its orthogroups"<< endl;
		return NULL;
	}
	map<string,STRINTMAP*>* hits=geneSet[sKey];
	string tKey(targetSpecies);
	if(hits->find(tKey)==hits->end())
	{
		return NULL;
	} 
	return (*hits)[tKey];
}

map<string,map<string,STRINTMAP*>*>&
GeneMap::getGeneSet()
{
	return geneSet;
}

#ifndef _GENETREEMGR_
#define _GENETREEMGR_
#include <map>
#define GENETREEPATH "data/TREES"
using namespace std;
class GeneTree;
class SpeciesDistManager;
class MappedOrthogroup;
class GeneTreeManager
{
	public:
	GeneTreeManager();
	~GeneTreeManager();
	int setSpeciesTree(SpeciesDistManager::Species*);
	int setSpeciesDistManager(SpeciesDistManager*);
	int setGeneTreeDir(const char*);
	GeneTree* getGeneTree(int);	
	GeneTree* getGeneTree(MappedOrthogroup*);	
	GeneTree* generateTree(MappedOrthogroup*);
	GeneTree* generateTreeFromFile(MappedOrthogroup*);
	inline int getTotalAdditionalParents()
	{
		return totalAdditionalParents;
	}
	int showTotalAdditionalParentDist();
	private:
	int showTree(GeneTree*);
	int pruneTree(GeneTree*);
	int insertTree(GeneTree* node);
	int setNodeType(GeneTree* node);
	int identifyAncestorWithDuplicates(SpeciesDistManager::Species* speciesnode,MappedOrthogroup* mor,SpeciesDistManager::Species** duplicateSpecies);
	const char* getSpeciesName(GeneTree* node);
	int getSpeciesName(GeneTree* node,map<string,int>&);
	GeneTree* createTree(SpeciesDistManager::Species* duplicateSpecies,SpeciesDistManager::Species* node, MappedOrthogroup* mor, int copyno);
	GeneTree* readTreeFromFile(const char* treePath);
	int addKpolGene(GeneTree*,string&);
	GeneTree* isKpolAnc(GeneTree* gtree,string& parentName);
	map<int,GeneTree*> genetreeSet;
	SpeciesDistManager::Species* speciestreeRoot;
	SpeciesDistManager* sdMgr;
	char dirName[1024];
	int totalAdditionalParents;
	map<string,int> addedParents;
};
#endif

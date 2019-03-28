#ifndef _GENETREE_
#define _GENETREE_
#include <string>
#include <vector>
using namespace std;
class GeneTree
{
	public:
	GeneTree();
	~GeneTree();
	//GeneTree* leftchild;
	//GeneTree* rightchild;
	vector<GeneTree*> children;
        GeneTree* parent;
	string species;
	string name;
	int nodeType;
	bool exclude;
};
#endif


#ifndef _GAMMA_
#define _GAMMA_
#include <map>
#include <string>
using namespace std;
class Matrix;
class GeneTree;
class Gamma
{
	public:
	Gamma();
	~Gamma();
	inline int setMaxClusterCnt(int k)
	{
		maxClusterCnt=k;
	}
	int addGene(string& geneName,string& specName,int maxClusterCnt);
	int setGeneClusterID(string&,string&,int);
	int initUsingGeneTree(GeneTree*);
	Matrix* getGamma(string&,string&);
	Matrix* getNormTerm(string&,string&);
	double getScore();
	struct Node;
	struct Node
	{
		string name;
		string species;
		//Node* leftchild;
		//Node* rightchild;
        vector<Node *> children;
		Node* parent;
		//This is 2 for duplication and 1 for speciation
		int nodeType;
		Matrix* normTerm;
		Matrix* gamma;
		bool hasData;
		map<int,double>* dataLL;
		map<int,double> mixProbs;
        bool excludeFlag;
	};
	Node* getRoot();
	Node* getGeneNode(string&,string&);
	Node* root;
	int showTree();
	map<int,double>* getDataLL(string&);
	string& getDupAncestor();
	private:
	Gamma::Node* initSubtree(GeneTree* gtreenode);
	Gamma::Node* initSubtree(GeneTree* gtreenode,string& rootname);
	int showTree(Node*);
	map<string,Node*> geneMap;
	map<string,Matrix*> gammas;
	map<string,map<int,double>*> datall;
	int maxClusterCnt;
	int nodeID;
	string dupAncName;
	
};
#endif

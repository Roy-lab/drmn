#ifndef _GENEEXPMGR_
#define _GENEEXPMGR_

#include <map>
#include <string>
#include <vector>

using namespace std;

class GeneExpManager
{
	public:
	GeneExpManager();
	~GeneExpManager();
	int readExpression(const char* aFName);
	int readExpression_Withheader(const char* aFName);
	vector<double>* getExp(const string&);
	map<string,vector<double>*>& getGeneSet();
	vector<string>& getColNames();
	
	private:
	int readColNames(char*);
	map<string,vector<double>*> exprSet;
	vector<string> colNames;
};
#endif

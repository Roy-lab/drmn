#ifndef _GENENAME_MAPPER
#define _GENENAME_MAPPER
#include <map>
#include <string>

using namespace std;

class GeneNameMapper
{
	public:
	GeneNameMapper();
	~GeneNameMapper();
	int readGeneNames();
	const char* getCommonName(const char*);

	private:
	map<string,string> orfToCommon;
};
#endif

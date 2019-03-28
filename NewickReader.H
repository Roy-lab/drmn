#ifndef _NEWICK_READER
#define _NEWICK_READER
class GeneTree;
class NewickReader
{
	public:
	NewickReader();
	~NewickReader();
	GeneTree* readTree(const char*);
	GeneTree* parseNewickFormat(char*);

	private:
	GeneTree* root;
};
#endif

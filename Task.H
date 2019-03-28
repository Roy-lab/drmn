#ifndef _TASK
#define _TASK

#include <map>
#include <string>
#include "Matrix.H"

using namespace std;

class Task_T
{
public:
	Task_T();
	~Task_T();
	Matrix* X;
	Matrix* Y;
	map<string,int>* YNames;
	string name;
};

#endif //_TASK

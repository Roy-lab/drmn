#ifndef _LEASTDIRTY
#define _LEASTDIRTY

#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
#include "Task.H"
#include "GenericLearner.H"

using namespace std;

class LeastDirty: public GenericLearner
{
	public:
		LeastDirty(double l,double b):GenericLearner(l,b)
		{
		}
		~LeastDirty()
		{
		}
		int doCV(int,int,int,vector<Task_T*>*,map<string,int>,const char*);
		int doAGM(vector<Task_T*>*,vector<Matrix*>&,vector<Matrix*>&,vector<Matrix*>&);
	private:
		int sum2vec(vector<Matrix*>,vector<Matrix*>,vector<Matrix*>&);
		int sub2vec(vector<Matrix*>,vector<Matrix*>,vector<Matrix*>&);
		double getError(vector<Matrix*> allXT, vector<Matrix*> allW, vector<Matrix*> allY);
		double getL1infnorm(vector<Matrix*>);
		double getL11norm(vector<Matrix*>);
		int proximalL1infnorm(vector<Matrix*>,vector<Matrix*>,double,vector<Matrix*>&);
		int proximalL11norm(vector<Matrix*>,vector<Matrix*>,double,vector<Matrix*>&);
		bool calcCond(vector<Matrix*>,vector<Matrix*>,vector<Matrix*>,double);
		int getPn(vector<Matrix*>,vector<Matrix*>,double,double,vector<Matrix*>&);
};

#endif //_LEASTDIRTY

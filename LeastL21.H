#ifndef _LEASTL21
#define _LEASTL21

#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
#include "Task.H"
#include "GenericLearner.H"

using namespace std;

class LeastL21: public GenericLearner
{
	public:
		LeastL21(double l,double b):GenericLearner(l,b)
		{
		}
		~LeastL21()
		{
		}
		int doAGM(vector<Task_T*>*,vector<Matrix*>&);
	private:
		double lossFunc(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>);
		int getGrad(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>&);
		int getFGLasso(vector<Matrix*>,double,vector<Matrix*>&);
		double getZPGamma(double,double,vector<Matrix*>,vector<Matrix*>);
		double getNonSmooth(vector<Matrix*>);
};

#endif //_LEASTL21

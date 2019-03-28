#ifndef _LEASTLASSO
#define _LEASTLASSO

#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
#include "Task.H"
#include "GenericLearner.H"

using namespace std;

class LeastLasso: public GenericLearner
{
	public:
		LeastLasso(double l,double b):GenericLearner(l,b)
		{
		}
		~LeastLasso()
		{
		}
		int doAGM(vector<Task_T*>*,vector<Matrix*>&);
	private:
		double lossFunc(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>);
		int getGrad(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>&);
		int getL1Proj(vector<Matrix*>,double,vector<Matrix*>&,double&);
		double getZPGamma(double,double,vector<Matrix*>,vector<Matrix*>);
};

#endif //_LEASTLASSO

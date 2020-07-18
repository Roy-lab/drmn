#ifndef _LEASTCFGLASSO
#define _LEASTCFGLASSO

#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
#include "Task.H"
#include "GenericLearner.H"

using namespace std;

class LeastCFGLasso: public GenericLearner
{
	public:
		LeastCFGLasso(double r1, double r2, double r3);
		~LeastCFGLasso()
		{
		}
		int doAGM(vector<Task_T*>*,vector<Matrix*>&,map<int,vector<int>>&);
	private:
		int getFGLasso(vector<Matrix*>&,double,double,double,vector<Matrix*>&);
		int supportSet(double*,double*,double*,double*,int*,double,int);
		double dualityGap2(double*,double*,double*,double*,double,int);
		int generateSolution(double*,double*,double*,double*,double*,double*,double*,int*,double,int);
		int sfa_one(double*,double*,int*,double*,double*,double*,double,int,int,double*,double*,double,int);
		double Thomas(double*,double*,int);
		int flsa(double*,double*,double*,double,double,int,int,double,int,int);
		double getNonSmooth(vector<Matrix*>,map<int,vector<int>>&);
		int getGrad(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>&);
		double lossFunc(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>);

		double rho1;
		double rho2;
		double rho3;
		double delta;
};

#endif //_LEASTCFGLASSO

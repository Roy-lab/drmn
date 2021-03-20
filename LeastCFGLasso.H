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
		//int getFGLasso(vector<Matrix*>&,double,double,double,vector<Matrix*>&,Matrix*,Matrix*);
		int getFGLasso_generic(vector<Matrix*>&,double,double,double,vector<Matrix*>&,Matrix*,Matrix*,Matrix*,double);
		int supportSet(double*,double*,double*,double*,int*,double,int);
		int supportSet(double*,double*,double*,double*,int*,double,int,Matrix*);
		double dualityGap2(double*,double*,double*,double*,double,int);
		double dualityGap2_generic(double*,double*,double*,double*,double,int);
		int generateSolution(double*,double*,double*,double*,double*,double*,double*,int*,double,int);
		int sfa_one(double*,double*,int*,double*,double*,double*,double,int,int,double*,double*,double,int,Matrix*,Matrix*);
		int sfa_one_SFAg(double*,double*,int*,double*,double*,double*,double,int,int,double*,double*,double,int,Matrix*,Matrix*,Matrix*,double);
		int sfa_one_SFAg_take2(double*,double*,int*,double*,double*,double*,double,int,int,double*,double*,double,int,Matrix*,Matrix*,Matrix*,double);
		double Thomas(double*,double*,int,Matrix*,Matrix*);
		//int flsa(double*,double*,double*,double,double,int,int,double,int,int,Matrix*,Matrix*);
		int flsa_generic(double*,double*,double*,double,double,int,int,double,int,int,Matrix*,Matrix*,Matrix*,double);
		double getNonSmooth(vector<Matrix*>,map<int,vector<int>>&);
		int getGrad(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>,vector<Matrix*>&);
		double lossFunc(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>);


		//int showAllWeights(int,int,vector<Matrix*>&);

		double rho1;
		double rho2;
		double rho3;
		double delta;

		Matrix* getTree(map<int,vector<int>>&,int);
};

#endif //_LEASTCFGLASSO

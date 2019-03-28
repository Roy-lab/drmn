#ifndef _GENERICLEARNER
#define _GENERICLEARNER

#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
#include "Task.H"

using namespace std;

class GenericLearner
{
	public:
		GenericLearner(double,double);
		~GenericLearner();
		virtual int doCV(int,int,int,vector<Task_T*>*,map<string,int>,const char*);
		virtual int doAGM(vector<Task_T*>*,vector<Matrix*>&);
		virtual int doAGM(vector<Task_T*>*,vector<Matrix*>&,vector<Matrix*>&,vector<Matrix*>&);
	protected:
		int maxItr;
		double beta;
		double lambda;
		double tol;
		int populateRandIntegers(gsl_rng*,int*,int);
		int split(int,int,int,Task_T*,Task_T*&,Task_T*&);

		int clearMatVector(vector<Matrix*>&);
		int copyMatVectorA2B(vector<Matrix*>,vector<Matrix*>&);

		int predict(vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>&);
		int writeAllPred(const char*,vector<Task_T*>*,vector<Task_T*>*,vector<Matrix*>,vector<Matrix*>);
		int writePred(const char*,map<string,int>*,Matrix*,Matrix*);
		int writeRegs(const char*,const char*,vector<Task_T*>*,vector<Matrix*>,map<string,int>);

		int showMe(const char*,vector<Matrix*>);
};

#endif //_GENERICLEARNER

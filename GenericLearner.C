// ************************************************************
// This class just process the data. For actual algorithms see:
// LeastLasso, LeastL21, and LeastDirty
// ************************************************************

#include "GenericLearner.H"
#include <algorithm>
#include <fstream>

using namespace std;

GenericLearner::GenericLearner()
{
	lambda = 0;
	beta   = 0;
	maxItr = 1000;
	tol    = 1e-4;
}

GenericLearner::GenericLearner(double l, double b)
{
	lambda = l;
	beta   = b;
	maxItr = 1000;
	tol    = 1e-4;
}

GenericLearner::~GenericLearner()
{
}

int
GenericLearner::split(int folds, int myFold, int foldSeed, Task_T* onet, Task_T*& traint, Task_T*& testt)
{
	int n  = onet->X->getColCnt();

	// this is where we populate randInds with the order of examples
	gsl_rng* randnum;
	randnum = gsl_rng_alloc(gsl_rng_default);  
	gsl_rng_set(randnum, foldSeed);
	int* randInds = new int[n];
	populateRandIntegers(randnum,randInds,n);

	int foldSize=(int)(n/folds); 
	int testBegin=myFold*foldSize;
	int testEnd=testBegin+foldSize;

	vector<int> trainIDs;
	vector<int> testIDs;
	for(int i=0;i<n;i++)
	{
		int randid=randInds[i];
		if(i>=testBegin && i<testEnd)
		{
			testIDs.push_back(randid);
		}
		else
		{
			trainIDs.push_back(randid);
		}
	}
	Matrix* Xtrain = onet->X->copyMe();Xtrain->removeCols(testIDs);
	Matrix* Ytrain = onet->Y->copyMe();Ytrain->removeRows(testIDs);
	Matrix* Xtest  = onet->X->copyMe();Xtest->removeCols(trainIDs);
	Matrix* Ytest  = onet->Y->copyMe();Ytest->removeRows(trainIDs);
	map<string,int>* YNamestrain = new map<string,int>;
	map<string,int>* YNamestest  = new map<string,int>;
	map<int,string>  id2name;
	for (map<string,int>::iterator siter = onet->YNames->begin(); siter != onet->YNames->end(); siter++)
	{
		id2name[siter->second] = siter->first;
	}
	int trainc=0;
	int testc=0;
	for (int i=0;i<n;i++)
	{
		if (find(trainIDs.begin(),trainIDs.end(),i) != trainIDs.end())
		{
			(*YNamestrain)[id2name[i]] = trainc;
			trainc++;
		}
		else
		{
			(*YNamestest)[id2name[i]] = testc;
			testc++;
		}
	}
	traint = new Task_T;
	traint->X = Xtrain;
	traint->Y = Ytrain;
	traint->YNames = YNamestrain;
	traint->name = onet->name;
	testt  = new Task_T;
	testt->X = Xtest;
	testt->Y = Ytest;
	testt->YNames = YNamestest;
	testt->name = onet->name;
	return 0;
}

int
GenericLearner::doCV(int folds, int myFold, int foldSeed, vector<Task_T*>* allt,map<string,int>regs, const char* outdir)
{
	//standardize
	for (int i=0;i<allt->size();i++)
	{
		Task_T* t = allt->at(i);
		t->X->rowStandardize();
		t->Y->colStandardize();
	}
	vector<Task_T*>* alltrain = new vector<Task_T*>;
	vector<Task_T*>* alltest  = new vector<Task_T*>;
	for (int i=0;i<allt->size();i++)
	{
		Task_T* onet = allt->at(i);
		Task_T* traint;
		Task_T* testt;
		split(folds, myFold, foldSeed, onet, traint, testt);
		alltrain->push_back(traint);
		alltest->push_back(testt);
	}
	vector<Matrix*> predTrain;
	vector<Matrix*> predTest;
	vector<Matrix*> allW;
	doAGM(alltrain, allW);
	for (int i=0;i<allW.size();i++)
	{
		int cnt =0;
		Matrix* W = allW[i];
		int l = W->getRowCnt();
		for (int j=0;j<l;j++)
		{
			double v = W->getValue(j,0);
			if (fabs(v)>0.001)
			{
				cnt++;
			}
		}
		cout << "Task" << i << " " << cnt << "/" << l << ",";
	}
	cout << endl;
	cout << "predict training..." << endl;
	predict(alltrain,allW,predTrain);
	cout << "predict test..." << endl;
	predict(alltest,allW,predTest);
	cout << "write prediction..." << endl;
	writeAllPred(outdir, alltrain, alltest, predTrain, predTest);
	cout << "write regulators..." << endl;
	writeRegs(outdir,"regs",alltrain, allW, regs);
	return 0;
}

int
GenericLearner::writeRegs(const char* outdir, const char* suff, vector<Task_T*>* alltrain, vector<Matrix*> allW, map<string,int> regs)
{
	char outname[1024];
	for (int i=0;i<alltrain->size();i++)
	{
		Task_T* traint = alltrain->at(i);
		Matrix* W = allW[i];
		sprintf(outname,"%s/%s.%s.txt",outdir,traint->name.c_str(),suff);
		ofstream outFile(outname);
		int l = W->getRowCnt();
		for (map<string,int>::iterator ritr=regs.begin();ritr!=regs.end();ritr++)
		{
			string r = ritr->first;
			int j = ritr->second;
			double v = W->getValue(j,0);
			if (fabs(v)>0.001)
			{
				outFile << r << "\t" << v << endl;
			}
		}
		outFile.close();
	}
}

int 
GenericLearner::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}

int
GenericLearner::writeAllPred(const char* outdir, vector<Task_T*>* alltrain, vector<Task_T*>* alltest, vector<Matrix*> predTrain, vector<Matrix*> predTest)
{
	char outname[1024];
	for (int i=0;i<alltrain->size();i++)
	{
		Task_T* traint = alltrain->at(i);
		Task_T* testt  = alltest->at(i);
		Matrix* YTrain = traint->Y;
		Matrix* YTest  = testt->Y;
		Matrix* PTrain = predTrain[i];
		Matrix* PTest  = predTest[i];
		map<string,int>* YTrainNames = traint->YNames;
		map<string,int>* YTestNames  = testt->YNames;
		sprintf(outname,"%s/%s.train.txt",outdir,traint->name.c_str());
		writePred(outname,YTrainNames, YTrain, PTrain);
		sprintf(outname,"%s/%s.test.txt",outdir,testt->name.c_str());
		writePred(outname,YTestNames, YTest, PTest);
	}
	return 0;
}

int
GenericLearner::writePred(const char* outname,map<string,int>* YNames, Matrix* Y, Matrix* P)
{
	ofstream outFile(outname);
	for (map<string,int>::iterator gitr=YNames->begin();gitr!=YNames->end();gitr++)
	{
		string gene = gitr->first;
		int index   = gitr->second;
		double y = Y->getValue(index,0);
		double p = P->getValue(index,0);
		outFile << gene << "\t" << y << "\t" << p << endl;
	}
	outFile.close();
	return 0;
}

int
GenericLearner::clearMatVector(vector<Matrix*>& vec)
{
	for(int i=0;i<vec.size();i++)
	{
		Matrix* M = vec[i];
		delete M;
	}
	vec.clear();
	return 0;
}

int
GenericLearner::copyMatVectorA2B(vector<Matrix*> allA, vector<Matrix*>& allB)
{
	clearMatVector(allB);
	for (int i=0;i<allA.size();i++)
	{
		Matrix* A = allA[i];
		Matrix* B = A->copyMe();
		allB.push_back(B);
	}
	return 0;
}

int
GenericLearner::predict(vector<Task_T*>* allt, vector<Matrix*> allW, vector<Matrix*>& allPY)
{
	clearMatVector(allPY);
	for (int i=0;i<allt->size();i++)
	{
		Task_T* t = allt->at(i);
		Matrix* X = t->X;
		Matrix* T = X->transMatrix();
		Matrix* W = allW[i];
		Matrix* P = T->multiplyMatrix(W);
		allPY.push_back(P);
		delete T;
	}
	return 0;
}

int
GenericLearner::doAGM(vector<Task_T*>* allt, vector<Matrix*>& allW, map<int,vector<int>>& tree)
{
	return 0;
}

int
GenericLearner::doAGM(vector<Task_T*>* allt, vector<Matrix*>& allW, vector<Matrix*>& allP, vector<Matrix*>& allQ)
{
	return 0;
}

int
GenericLearner::doAGM(vector<Task_T*>* allt, vector<Matrix*>& allW)
{
	return 0;
}

int
GenericLearner::showMe(const char* pref, vector<Matrix*> allA)
{
	cout << pref << endl;
	int row=allA[0]->getRowCnt();
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<allA.size();j++)
		{
			double v = allA[j]->getValue(i,0);
			printf("%10.4f",v);
		}
		cout << endl;
	}
	return 0;
}

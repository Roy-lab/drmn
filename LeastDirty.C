// ******************************************************
// This is c++ implementation of Least_Dirty fomr MALSAR
// https://github.com/jiayuzhou/MALSAR/blob/master/MALSAR/functions/dirty/Least_Dirty.m
// It is using accelerated gradient method (AGM) to solve 
// the optimization
// ******************************************************

#include "LeastDirty.H"
#include <algorithm>
#include <fstream>
#include <stdio.h>

using namespace std;

//I am redefining this to pass P and Q to doAGM. 
//Is there a better way to do this?
int
LeastDirty::doCV(int folds, int myFold, int foldSeed, vector<Task_T*>* allt,map<string,int>regs, const char* outdir)
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
	vector<Matrix*> allP;
	vector<Matrix*> allQ;
	doAGM(alltrain, allW, allP, allQ);
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
	writeRegs(outdir,"regs",alltrain,allW,regs);
	writeRegs(outdir,"Pregs",alltrain,allP,regs);
	writeRegs(outdir,"Qregs",alltrain,allQ,regs);
	return 0;
}

int
LeastDirty::doAGM(vector<Task_T*>* allt, vector<Matrix*>& allW, vector<Matrix*>& allP, vector<Matrix*>& allQ)
{
	vector<Matrix*> allXY;
	vector<Matrix*> allXT;
	vector<Matrix*> allXX;
	vector<Matrix*> allP_old;
	vector<Matrix*> allQ_old;
	vector<Matrix*> allPn;
	vector<Matrix*> allQn;
	vector<Matrix*> allDP;
	vector<Matrix*> allDQ;
	vector<Matrix*> allGW;
	vector<Matrix*> allY; // just to keep them outside Task_T
	int task_cnt = allt->size();
	int gene_cnt;
	int reg_cnt;
	for(int i=0;i<allt->size();i++)
	{
		Task_T* t  = allt->at(i);
		Matrix* X  = t->X;
		Matrix* Y  = t->Y;
		Matrix* XY = X->multiplyMatrix(Y);
		Matrix* XT = X->transMatrix();
		Matrix* XX = X->multiplyMatrix(XT);
		allXY.push_back(XY);
		allXT.push_back(XT);
		allXX.push_back(XX);
		int row    = X->getRowCnt();
		Matrix* P0 = new Matrix(row,1);
		allP.push_back(P0);

		allY.push_back(Y);

		gene_cnt = X->getColCnt();
		reg_cnt  = X->getRowCnt();
	}
	copyMatVectorA2B(allP,allQ);
	copyMatVectorA2B(allP,allPn);
	copyMatVectorA2B(allP,allQn);

	int t_new = 1;
	int t_old = 1;
	//L1norm = max(sum(abs(X),1));   //sum columns of correct X, rows of my X
	//Linfnorm = max(sum(abs(X),2)); //sum rows of correct X, cols of my X
	double L1norm = 0;
	double Linfnorm = 0;
	double cmax = 0;
	double rmax = 0;
	for(int t=0;t<allt->size();t++)
	{
		Task_T* T  = allt->at(t);
		Matrix* X  = T->X;
		int row = X->getRowCnt();
		int col = X->getColCnt();
		for (int i=0;i<row;i++)
		{
			double rsum = 0;
			for (int j=0;j<col;j++)
			{
				double v = X->getValue(i,j);
				rsum += fabs(v);
			}
			if (rsum>rmax)
			{
				rmax = rsum;
			}
		}
		for (int j=0;j<col;j++)
		{
			double csum = 0;
			for (int i=0;i<row;i++)
			{
				double v = X->getValue(i,j);
				csum += fabs(v);
			}
			if (csum>cmax)
			{
				cmax = csum;
			}
		}
	}
	L1norm = rmax;
	Linfnorm = cmax;
	//L = 2*max(L1norm*L1norm/size(X,1),Linfnorm*Linfnorm/size(X,2));
	double L1 = (L1norm*L1norm)/(gene_cnt*task_cnt);
	double L2 = (Linfnorm*Linfnorm)/(reg_cnt*task_cnt);
	double L  = L1;
	if (L<L2)
	{
		L=L2;
	}
	L = L*2;
	
	//funcVal = cat(1, funcVal, norm(X*(P(:) + Q(:)) - y)^2 + lambda1*L1infnorm(P) + lambda2*L11norm(Q));
	double fVal = 0;
	double fVal_old = 0;
	sum2vec(allP,allQ,allW);
	double e = getError(allXT, allW, allY);
	L1 = getL1infnorm(allP);
	L2 = getL11norm(allQ);
	fVal = e+lambda*L1+beta*L2;

	for(int itr=0;itr<maxItr;itr++)
	{
		//if (itr>3)
		//{
		//	exit(0);
		//}
		//showMe("P0",allP);
		//showMe("Q0",allQ);
		//cout << "P0" << endl;
		//allP[0]->showMatrix();
		//cout << "Q0" << endl;
		//allQ[0]->showMatrix();
		//cout << endl;
		copyMatVectorA2B(allP,allP_old);
		copyMatVectorA2B(allQ,allQ_old);
		t_old = t_new;

		cout << "itr:" << itr << ":" ; // << endl;
		//gradvec = 2*(XtX*(Pn(:)+Qn(:)) - Xty);
		//gradmat = reshape(gradvec,dimension,task_num);
		clearMatVector(allGW);
		for (int t=0;t<allXX.size();t++)
		{
			Matrix* XX = allXX[t];
			Matrix* Pn = allPn[t];
			Matrix* Qn = allQn[t];
			Matrix* XY = allXY[t];
			Matrix* W  = Pn->addMatrix(Qn);
			Matrix* GW = XX->multiplyMatrix(W);
			GW->subtractWithMatrix(XY);
			GW->multiplyScalar(2);
			allGW.push_back(GW);
			delete W;
		}
		//showMe("GW",allGW);

		for (int initr=0;initr<20;initr++)
		{
			cout << ".";
			//P = proximalL1infnorm(Pn - gradmat/L, lambda1/L);
			//Q = proximalL11norm(Qn - gradmat/L, lambda2/L);
			proximalL1infnorm(allPn, allGW, L, allP);
			proximalL11norm(allQn, allGW, L, allQ);
			//dP = P - Pn;  dQ = Q - Qn;
			sub2vec(allP,allPn,allDP);
			sub2vec(allQ,allQn,allDQ);
			//if 2*((dP(:) + dQ(:))'*XtX*(dP(:) + dQ(:))) <= L*sum(sum((dP.*dP + dQ.*dQ)))
			//	break;
			//else
			//	L = L*2;
			//end
			//showMe("P:",allP);
			//showMe("Q:",allQ);
			bool flag = calcCond(allDP, allDQ, allXX, L);
			if (flag)
			{
				break;
			}
			else
			{
				L = L*2;
			}
		}
		//cout << endl;

		fVal_old = fVal;
		sum2vec(allP,allQ,allW);
		e = getError(allXT, allW, allY);
		L1 = getL1infnorm(allP);
		L2 = getL11norm(allQ);
		fVal = e+lambda*L1+beta*L2;
		
		if (itr>=2)
		{
			if ( fabs(fVal-fVal_old) <= tol*fVal_old )
			{
				break;
			}
		}
		//t_new = (1+sqrt(1+4*t_old^2))/2;
		t_new = (1+sqrt(1+4*t_old*t_old))/2;

		//Pn = P + (t_old-1)/t_new*(P - P_old);
		getPn(allP, allP_old, t_old, t_new, allPn);
		//Qn = Q + (t_old-1)/t_new*(Q - Q_old);
		getPn(allQ, allQ_old, t_old, t_new, allQn);
	}
	cout << endl;
	
	sum2vec(allP,allQ,allW);

	clearMatVector(allPn);
	clearMatVector(allQn);
	clearMatVector(allP_old);
	clearMatVector(allQ_old);
	clearMatVector(allXT);
	clearMatVector(allXY);
	allY.clear();//don't delete Y Matrices
	return 0;
}

bool
LeastDirty::calcCond(vector<Matrix*> allDP, vector<Matrix*> allDQ, vector<Matrix*> allXX, double L)
{
	double dsum = 0;
	for (int t=0;t<allDP.size();t++)
	{
		Matrix* DP = allDP[t];
		Matrix* DQ = allDQ[t];
		int row = DP->getRowCnt();
		for (int i=0;i<row;i++)
		{
			double dp = DP->getValue(i,0);
			double dq = DQ->getValue(i,0);
			dsum += (dp*dp+dq*dq);
		}
	}
	double xsum = 0;
	for (int t=0;t<allDQ.size();t++)
	{
		Matrix* XX = allXX[t];
		Matrix* DP = allDP[t];
		Matrix* DQ = allDQ[t];
		Matrix* W  = DP->addMatrix(DQ);
		Matrix* WT = W->transMatrix();
		Matrix* R  = WT->multiplyMatrix(XX);
		R->multiplyWithMatrix(W);
		double v = R->getValue(0,0);
		xsum += v;
		delete W;
		delete WT;
		delete R;
	}
	//2*((dP(:) + dQ(:))'*XtX*(dP(:) + dQ(:))) <= L*sum(sum((dP.*dP + dQ.*dQ)))
	if (2*xsum <= L*dsum)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int
LeastDirty::proximalL1infnorm(vector<Matrix*> allPn, vector<Matrix*> allGW, double L, vector<Matrix*>& allP)
{
	clearMatVector(allP);
	//P = proximalL1infnorm(Pn - gradmat/L, lambda1/L);
	//function [X] = proximalL1infnorm(D, tau)
	//[m,n]=size(D);
	//[mu,~,~]=prf_lbm(D,m,n,tau);
	//X = D - mu;
	double tau = lambda/L;
	int i, j, gnum, rho = 0;
	double theta = 0;
	double s = 0;
	double iter_step=0; 

	double ltol = 1e-8;

	int m = allPn[0]->getRowCnt();
	int n = allPn.size();

	Matrix* D = new Matrix(m,n);
	Matrix* X = new Matrix(m,n);
	for (int t=0;t<allPn.size();t++)
	{
		Matrix* Pn = allPn[t];
		Matrix* GW = allGW[t];
		for (int i=0;i<m;i++)
		{
			double pn = Pn->getValue(i,0);
			double gw = GW->getValue(i,0);
			D->setValue(pn - gw/L,i,t);
		}
	}

	for (int j=0;j<m;j++)
	{
		gnum = 0; theta = 0; s = 0;
		iter_step = 0; rho = 0;
		for(int i=0;i<n;i++)
		{
			//x[gnum*m+j] = fabs(c[i*m+j]);
			double v = D->getValue(j,i);
			v = fabs(v);
			X->setValue(v,j,gnum);
			//s += x[gnum*m+j];
			s += v;
			gnum++;
		}
		/* If ||c||_1 <= t, then c is the solution  */
		if (s <= tau)
		{
			theta=0;
			for(int i=0;i<n;i++)
			{
				//x[i*m+j]=c[i*m+j];
				double v = D->getValue(j,i);
				X->setValue(v,j,i);
			}
			continue;
		}
		/*while loops*/
		while (fabs(s - tau - gnum*theta) > ltol)
		{
			iter_step++;
			theta = (s-tau)/gnum;
			s=0; rho = 0;
			for (int i=0;i<gnum;i++)
			{
				double v = X->getValue(j,i);
				//if (x[i*m+j] >= theta)
				if (v >= theta)
				{
					//x[rho*m+j] = x[i*m+j];
					X->setValue(v,j,rho);
					//s+=x[i*m+j]; 
					s+=v; 
					rho++;
				}
			}
			gnum = rho;
		}
		/*end of while*/
		/*projection result*/
		for(int i=0;i<n;i++)
		{
			//x[i*m+j] = (c[i*m+j] > theta)?(c[i*m+j]-theta):((c[i*m+j]< -theta)?(c[i*m+j]+theta):0);
			double v = D->getValue(j,i);
			if (v>theta)
			{
				X->setValue(v-theta,j,i);
			}
			else
			{
				if (v < -theta)
				{
					X->setValue(v+theta,j,i);
				}
				else
				{
					X->setValue(0,j,i);
				}
			}
		}
	}
	//X = D - mu;
	for (int t=0;t<allPn.size();t++)
	{
		Matrix* P = new Matrix(m,1);
		for (int i=0;i<m;i++)
		{
			double d = D->getValue(i,t);
			double x = X->getValue(i,t);
			P->setValue(d-x,i,0);
		}
		allP.push_back(P);
	}
	delete X;
	delete D;
	return 0;
}

int
LeastDirty::proximalL11norm(vector<Matrix*> allQn, vector<Matrix*> allGW, double L, vector<Matrix*>& allQ)
{
	//Q = proximalL11norm(Qn - gradmat/L, lambda2/L);
	//function [X] = proximalL11norm(D, tau)
	//X = sign(D).*max(0,abs(D)-tau);
	clearMatVector(allQ);
	double tau = beta/L;
	for (int t=0;t<allQn.size();t++)
	{
		Matrix* Qn = allQn[t];
		Matrix* GW = allGW[t];
		Matrix* D  = GW->copyMe();
		D->multiplyScalar(-1/L);
		D->addWithMatrix(Qn);
		int row = D->getRowCnt();
		Matrix* Q = new Matrix(row,1);
		for (int i=0;i<row;i++)
		{
			double di = D->getValue(i,0);
			int s = 1;
			if (di<0)
			{
				s=-1;
			}
			di = fabs(di);
			di = di-tau;
			if (di<0)
			{
				di = 0;
			}
			Q->setValue(s*di,i,0);
		}
		allQ.push_back(Q);
	}
	return 0;
}

int
LeastDirty::getPn(vector<Matrix*> allP, vector<Matrix*> allP_old, double t_old, double t_new, vector<Matrix*>& allPn)
{
	//Pn = P + (t_old-1)/t_new*(P - P_old);
	clearMatVector(allPn);
	for (int i=0;i<allP.size();i++)
	{
		Matrix* P = allP[i];
		Matrix* P_old = allP_old[i];
		Matrix* Pn = P->copyMe();
		Pn->subtractWithMatrix(P_old);
		Pn->multiplyScalar((t_old-1)/t_new);
		Pn->addWithMatrix(P);
		allPn.push_back(Pn);
	}
	return 0;
}

double
LeastDirty::getL1infnorm(vector<Matrix*> allP)
{
	//Xnorm = sum(max(abs(X),[],2)); //max of rows, sum of that
	double res = 0;
	int row = allP[0]->getRowCnt();
	for (int i=0;i<row;i++)
	{
		double rmax = 0;
		for (int t=0;t<allP.size();t++)
		{
			Matrix* P = allP[t];
			double  v = P->getValue(i,0);
			v = fabs(v);
			if (v>rmax)
			{
				rmax = v;
			}
		}
		res = res+rmax;
	}
	return res;
}

double
LeastDirty::getL11norm(vector<Matrix*> allQ)
{
	//Xnorm = sum(sum(abs(X)));
	double res = 0;
	for (int t=0;t<allQ.size();t++)
	{
		Matrix* Q = allQ[t];
		int row = Q->getRowCnt();
		int col = Q->getColCnt();
		for (int i=0;i<row;i++)
		{
			for (int j=0;j<col;j++)
			{
				double v = Q->getValue(i,j);
				res += fabs(v);
			}
		}
	}
	return res;
}

int
LeastDirty::sub2vec(vector<Matrix*> allA, vector<Matrix*> allB, vector<Matrix*>& allC)
{
	clearMatVector(allC);
	for (int i=0;i<allA.size();i++)
	{
		Matrix* A = allA[i];
		Matrix* B = allB[i];
		Matrix* C = A->subtractMatrix(B);
		allC.push_back(C);
	}
	return 0;
}

int
LeastDirty::sum2vec(vector<Matrix*> allA, vector<Matrix*> allB, vector<Matrix*>& allC)
{
	clearMatVector(allC);
	for (int i=0;i<allA.size();i++)
	{
		Matrix* A = allA[i];
		Matrix* B = allB[i];
		Matrix* C = A->addMatrix(B);
		allC.push_back(C);
	}
	return 0;
}

double
LeastDirty::getError(vector<Matrix*> allXT, vector<Matrix*> allW, vector<Matrix*> allY)
{
	//norm(X*(P(:) + Q(:)) - y)^2
	double e=0;
	for (int i=0;i<allXT.size();i++)
	{
		Matrix* XT = allXT[i];
		Matrix* W  = allW[i];
		Matrix* Y  = allY[i];
		Matrix* E  = XT->multiplyMatrix(W);
		E->subtractWithMatrix(Y);
		double en  = E->getFNorm();
		e += (en*en);
	}
	return e;
}

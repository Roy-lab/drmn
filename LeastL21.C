// ******************************************************
// This is c++ implementation of Least_L21 fomr MALSAR
// https://github.com/jiayuzhou/MALSAR/blob/master/MALSAR/functions/joint_feature_learning/Least_L21.m
// It is using accelerated gradient method (AGM) to solve 
// the optimization
// ******************************************************

#include "LeastL21.H"
#include <algorithm>
#include <fstream>

using namespace std;

int
LeastL21::doAGM(vector<Task_T*>* allt, vector<Matrix*>& allW)
{
	vector<Matrix*> allXY;
	vector<Matrix*> allXT;
	vector<Matrix*> allV;
	vector<Matrix*> allGW;
	vector<Matrix*> allWs;
	vector<Matrix*> allWz;
	vector<Matrix*> allWz_old;
	vector<Matrix*> allWzp;
	vector<Matrix*> allDelta;
	for(int i=0;i<allt->size();i++)
	{
		Task_T* t  = allt->at(i);
		Matrix* X  = t->X;
		Matrix* Y  = t->Y;
		Matrix* XY = X->multiplyMatrix(Y);
		Matrix* XT = X->transMatrix();
		allXY.push_back(XY);
		allXT.push_back(XT);
		int row    = X->getRowCnt();
		Matrix* W0;
		W0 = new Matrix(row,1);
		allWz.push_back(W0);
		W0 = new Matrix(row,1);
		allWz_old.push_back(W0);
	}
	double t = 1;
	double t_old = 0;
	double gamma = 1;
	double gamma_inc = 2;

	double l1c_wzp=0;

	double fVal = 0;
	double fVal_old = 0;
	int bFlag=0;

	for(int itr=0;itr<maxItr;itr++)
	{
		cout << "itr:" << itr << ":" ; // << endl;
		double alpha = (t_old - 1) /t;
		//Ws   = (1 + alpha) * Wz - alpha * Wz_old;
		clearMatVector(allWs);
		for(int i=0;i<allWz.size();i++)
		{
			Matrix* Wz = allWz[i];
			Matrix* Wz_old = allWz_old[i];
			Matrix* Temp = Wz_old->copyMe();
			Temp->multiplyScalar(-1*alpha);
			Matrix* Ws = Wz->copyMe();
			Ws->multiplyScalar((1+alpha));
			Ws->addWithMatrix(Temp);
			delete Temp;
			allWs.push_back(Ws);
			if (i==0)
			{
				//cout << "itr" << itr << ", Ws0" << endl;
				//Ws->showMatrix();
			}
		}
		//gWs  = gradVal_eval(Ws);
		clearMatVector(allGW);
		getGrad(allt, allXT, allXY, allWs, allGW);
		if(true)
		{
			Matrix* GW = allGW[0];
			//cout << "itr" << itr << ", GW0" << endl;
			//GW->showMatrix();
		}
		//Fs   = funVal_eval  (Ws);
		double Fs = lossFunc(allt, allXT, allWs);
		double Fzp;
		while(true)
		{
			cout << ".";
			//Wzp = FGLasso_projection(Ws - gWs/gamma, rho1 / gamma);
			clearMatVector(allV);
			for(int i=0;i<allWs.size();i++)
			{
				Matrix* Ws = allWs[i];
				Matrix* GW = allGW[i];
				Matrix* V  = GW->copyMe();
				V->multiplyScalar(-1/gamma);
				V->addWithMatrix(Ws);
				allV.push_back(V);
			}
			clearMatVector(allWzp);
			getFGLasso(allV, lambda / gamma, allWzp);
			//Fzp = funVal_eval  (Wzp);
			Fzp = lossFunc(allt, allXT, allWzp);
			//delta_Wzp = Wzp - Ws;
			double r_sum = 0;
			clearMatVector(allDelta);
			for (int i=0;i<allWzp.size();i++)
			{
				Matrix* Wzp = allWzp[i];
				Matrix* Ws  = allWs[i];
				Matrix* DL  = Wzp->subtractMatrix(Ws);
				allDelta.push_back(DL);
				double dl = DL->getFNorm();
				r_sum += (dl*dl);
			}
			if (r_sum <= 1e-20)
			{
				bFlag = 1;
				break;
			}
			//Fzp_gamma = Fs + sum(sum(delta_Wzp .* gWs)) + gamma/2 * sum(sum(delta_Wzp.*delta_Wzp));
			//showVec("Delta",allDelta);
			//showVec("GW", allGW);
			double Fzp_gamma = getZPGamma(Fs, gamma, allDelta, allGW);
			if (Fzp <= Fzp_gamma)
			{
				break;
			}
			else
			{
				gamma = gamma * gamma_inc;
			}
		}
		//cout << endl;
		//Wz_old = Wz;
		copyMatVectorA2B(allWz,allWz_old);
		//Wz = Wzp;
		copyMatVectorA2B(allWzp,allWz);

		//funcVal = cat(1, funcVal, Fzp + rho1 * l1c_wzp);
		double ns = getNonSmooth(allWz);
		fVal = Fzp + ns;

		if (bFlag == 1)
		{
			break;
		}
		//if iter>=2
		//	if (abs( funcVal(end) - funcVal(end-1) ) <= opts.tol* funcVal(end-1))
		//		break;
		//	end
		//end
		if (itr>=2)
		{
			if ( fabs(fVal-fVal_old) <= tol*fVal_old )
			{
				break;
			}
		}

		t_old = t;
		t = 0.5 * (1 + sqrt(1+ 4 * t*t));
		fVal_old = fVal;
	}
	cout << endl;

	copyMatVectorA2B(allWzp,allW);
	clearMatVector(allXT);
	clearMatVector(allXY);
	clearMatVector(allWzp);
	clearMatVector(allWz);
	clearMatVector(allWz_old);
	clearMatVector(allV);
	clearMatVector(allGW);
	clearMatVector(allWs);
	clearMatVector(allDelta);
	return 0;
}

double
LeastL21::getZPGamma(double Fs, double gamma, vector<Matrix*> allDelta, vector<Matrix*> allGW)
{
	//Fzp_gamma = Fs + sum(sum(delta_Wzp .* gWs)) + gamma/2 * sum(sum(delta_Wzp.*delta_Wzp));
	double Fzp_gamma = 0;
	Fzp_gamma += Fs;
	double s1 = 0;
	double s2 = 0;
	for (int i=0;i<allDelta.size();i++)
	{
		Matrix* DL = allDelta[i];
		Matrix* GW = allGW[i];
		int l = DL->getRowCnt();
		for (int j=0;j<l;j++)
		{
			double dwz = DL->getValue(j,0);
			double gw  = GW->getValue(j,0);
			Fzp_gamma += ((dwz*gw) + (gamma/2)*(dwz*dwz));
			s1 = s1+(dwz*gw);
			s2 = s2+(dwz*dwz);
		}
	}
	return Fzp_gamma;
}

int
LeastL21::getFGLasso(vector<Matrix*> allV, double v, vector<Matrix*>& allZ)
{
	//X = repmat(max(0, 1 - lambda./sqrt(sum(D.^2,2))),1,size(D,2)).*D;
	int l=allV[0]->getRowCnt();
	Matrix* T = new Matrix(l,1);
	for (int j=0;j<l;j++)
	{
		double tj = 0;
		for (int i=0;i<allV.size();i++)
		{
			Matrix* V = allV[i];
			double vj = V->getValue(j,0);
			tj += (vj*vj);
		}
		double res = 1 - v/sqrt(tj);
		if (res<0)
		{
			res = 0;
		}
		T->setValue(res,j,0);
	}
	allZ.clear();
	for (int i=0;i<allV.size();i++)
	{
		Matrix* V = allV[i];
		int l = V->getRowCnt();
		Matrix* Z = new Matrix(l,1);
		for (int j=0;j<l;j++)
		{
			double vj = V->getValue(j,0);
			double tj = T->getValue(j,0);
			Z->setValue(vj*tj,j,0);
		}
		allZ.push_back(Z);
	}
	return 0;
}

int
LeastL21::getGrad(vector<Task_T*>* allt, vector<Matrix*> allXT, vector<Matrix*> allXY, vector<Matrix*> allW, vector<Matrix*>& allGW)
{
	//grad_W = [];
	//for i = 1:task_num
	//	grad_W = cat(2, grad_W, X{i}*(X{i}' * W(:,i)-Y{i}) );
	//	grad_W = cat(2, grad_W, X{i}*(X{i}' * W(:,i)) - X{i}*Y{i} );
	//end
	//grad_W = grad_W+ rho_L2 * 2 * W;


	allGW.clear();
	for (int i=0;i<allt->size();i++)
	{
		Task_T* t  = allt->at(i);
		Matrix* W  = allW[i];
		Matrix* X  = t->X;
		//Matrix* Y  = t->Y;
		//Matrix* XT = X->transMatrix();
		Matrix* XT = allXT[i];
		Matrix* XW = XT->multiplyMatrix(W);
		Matrix* R  = X->multiplyMatrix(XW);
		//Matrix* XY = X->multiplyMatrix(Y);
		Matrix* XY = allXY[i];
		Matrix* GW = R->subtractMatrix(XY);
		Matrix* S  = W->copyMe();
		S->multiplyScalar(2*beta);
		GW->addWithMatrix(S);
		allGW.push_back(GW);
		delete XW;
		delete R;
		delete S;
		//delete XY;
		//delete XT;
	}
	return 0;
}

double
LeastL21::lossFunc(vector<Task_T*>* allt, vector<Matrix*> allXT, vector<Matrix*> allW)
{
	double funcVal = 0;
	for (int i=0;i<allt->size();i++)
	{
		Task_T* t  = allt->at(i);
		Matrix* W  = allW[i];
		Matrix* X  = t->X;
		Matrix* Y  = t->Y;
		//Matrix* XT = X->transMatrix();
		Matrix* XT = allXT[i];
		Matrix* R  = XT->multiplyMatrix(W);
		Matrix* S  = Y->subtractMatrix(R);
		double  f  = S->getFNorm();
		double  wn = W->getFNorm();
		funcVal = funcVal + (0.5*f*f) + (beta*wn*wn);
		//delete XT;
		delete R;
		delete S;
	}
	return funcVal;
}

double 
LeastL21::getNonSmooth(vector<Matrix*> allW)
{
	double res=0;
	for(int i=0;i<allW.size();i++)
	{
		Matrix* W = allW[i];
		double wn = W->getFNorm();
		res += wn*lambda;
	}
	return res;
}

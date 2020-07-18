//See https://github.com/jiayuzhou/MALSAR/blob/master/MALSAR/functions/progression_model/cFSGL/Least_CFGLasso.m
#include "LeastCFGLasso.H"

LeastCFGLasso::LeastCFGLasso(double r1, double r2, double r3)
{
	rho1 = r1;
	rho2 = r2;
	rho3 = r3;

	maxItr = 1000;
	tol    = 1e-4;
	delta  = 1e-10;
}

int 
LeastCFGLasso::doAGM(vector<Task_T*>* allt, vector<Matrix*>& allW, map<int,vector<int>>& tree)
{
	vector<Matrix*> allXY;//cleared
	vector<Matrix*> allXT;//cleared
	vector<Matrix*> allV;//cleared
	vector<Matrix*> allGW;//cleared
	vector<Matrix*> allWs;//cleared
	vector<Matrix*> allWz;//cleared
	vector<Matrix*> allWz_old;//cleared
	vector<Matrix*> allWzp;//cleared
	vector<Matrix*> allDelta;//cleared
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
    	//alpha = (t_old - 1) /t;
		double alpha = (t_old - 1) /t;
		//Ws   = (1 + alpha) * Wz - alpha * Wz_old;
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
		}
		//gWs  = gradVal_eval(Ws);
		clearMatVector(allGW);
		getGrad(allt, allXT, allXY, allWs, allGW);
		//Fs   = funVal_eval  (Ws);
		double Fs = lossFunc(allt, allXT, allWs);
		double Fzp;
		while(true)
		{
			cout << ".";
        	//Wzp = FGLasso_projection(Ws - gWs/gamma, rho1 / gamma, rho2 / gamma, rho3 / gamma);
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
			getFGLasso(allV, rho1/gamma, rho2/gamma, rho3/gamma, allWzp);

			//Fzp = funVal_eval  (Wzp);
			Fzp = lossFunc(allt, allXT, allWzp);
			//delta_Wzp = Wzp - Ws;
        	//nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;
        	//Fzp_gamma = Fs + sum(sum(delta_Wzp .* gWs)) + gamma/2 * nrm_delta_Wzp;
			double nrm_delta_Wzp = 0;
			double Fzp_gamma = 0;
			clearMatVector(allDelta);
			for (int i=0;i<allWzp.size();i++)
			{
				Matrix* Wzp = allWzp[i];
				Matrix* Ws  = allWs[i];
				Matrix* DL  = Wzp->subtractMatrix(Ws);
				allDelta.push_back(DL);
				double dl = DL->getFNorm();
				nrm_delta_Wzp += (dl*dl);

        		//Fzp_gamma = Fs + sum(sum(delta_Wzp .* gWs)) + gamma/2 * nrm_delta_Wzp;
				Matrix* GW = allGW[i];
				for (int j=0;j<GW->getRowCnt();j++)
				{
					double gwj = GW->getValue(j,0);
					double dlj = DL->getValue(j,0);
					Fzp_gamma += (gwj*dlj);
				}
			}
        	//Fzp_gamma = Fs + sum(sum(delta_Wzp .* gWs)) + gamma/2 * nrm_delta_Wzp;
			Fzp_gamma += Fs;
			Fzp_gamma += ((gamma/2)*nrm_delta_Wzp);
        	//r_sum = nrm_delta_Wzp;
			double r_sum = nrm_delta_Wzp;
			if (r_sum <= 1e-20)
			{
				bFlag = 1;
				break;
			}
        	//if (Fzp <= Fzp_gamma)
            //	break;
        	//else
            //	gamma = gamma * gamma_inc;
        	//end
			if (Fzp <= Fzp_gamma)
			{
				break;
			}
			else
			{
				gamma = gamma * gamma_inc;
			}
		}
		//Wz_old = Wz;
		copyMatVectorA2B(allWz,allWz_old);
		//Wz = Wzp;
		copyMatVectorA2B(allWzp,allWz);

    	//funcVal = cat(1, funcVal, Fzp + nonsmooth_eval(Wz, rho1, rho2, rho3));
		double ns = getNonSmooth(allWz,tree);
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

	clearMatVector(allXY);
	clearMatVector(allXT);
	clearMatVector(allV);
	clearMatVector(allGW);
	clearMatVector(allWs);
	clearMatVector(allWz);
	clearMatVector(allWz_old);
	clearMatVector(allWzp);
	clearMatVector(allDelta);

	return 0;
}

int
LeastCFGLasso::getFGLasso(vector<Matrix*>& allV, double v1, double v2, double v3, vector<Matrix*>& allZ)
{
	int row = allV[0]->getRowCnt();
	int n   = allV.size();
	allZ.clear();
	for (int i=0;i<n;i++)
	{
		Matrix* Z = new Matrix(row,1);
		allZ.push_back(Z);
	}
	double* v  = new double [n];
	double* w1 = new double [n];
	double* w0 = new double [n-1];
	for (int j=0;j<row;j++)
	{
		for (int i=0;i<n;i++)
		{
			v[i] = allV[i]->getValue(j,0);
			w1[i] = 0;
		}
		for (int i=0;i<n-1;i++)
		{
			w0[i] = 0;
		}
		//v = W(i, :);
		//w0 = zeros(length(v)-1, 1);
		//w_1 = flsa(v, w0,  lambda_1, lambda_2, length(v), 1000, 1e-9, 1, 6);
		flsa(w1,v,w0,v1,v2,n,1000,1e-9,1,6);
		//nm = norm(w_1, 2);
		double nm = 0;
		for (int i=0;i<n;i++)
		{
			nm += w1[i]*w1[i];
		}
		nm = sqrt(nm);
		//if nm == 0
		//	w_2 = zeros(size(w_1));
		//else
		//	w_2 = max(nm - lambda_3, 0)/nm * w_1;
		//end
		if (nm > 0)
		{
			//w_2 = max(nm - lambda_3, 0)/nm * w_1;
			double t = nm-v3;
			if (t<0)
			{
				t = 0;
			}
			for (int i=0;i<n;i++)
			{
				w1[i] *= (t/nm);
			}
		}
		//Wp(i, :) = w';
		for (int i=0;i<n;i++)
		{
			allZ[i]->setValue(w1[i],j,0);
		}
	}
	delete [] v;
	delete [] w0;
	delete [] w1;
	
	return 0;
}

int 
LeastCFGLasso::supportSet(double *x, double *v, double *z, double *g, int * S, double lambda, int nn)
{
	int i, j, n=nn+1, numS=0;
	double temp;

	/*
	we first scan z and g to obtain the support set S
	*/

	/*numS: number of the elements in the support set S*/
	for(i=0;i<nn; i++)
	{
		if ( ( (z[i]==lambda) && (g[i] < delta) ) || ( (z[i]==-lambda) && (g[i] >delta) ) )
		{
			S[numS]=i;
			numS++;
		}
	}
	

	if (numS==0) /*this shows that S is empty*/
	{
		temp=0;
		for (i=0;i<n;i++)
		{
			temp+=v[i];
		}

		temp=temp/n;
		for(i=0;i<n;i++)
		{
			x[i]=temp;
		}
		return numS;
	}


    /*
	 Next, we deal with numS >=1
     */

	/*process the first block
	   j=0
	*/
	temp=0;
	for (i=0;i<=S[0]; i++)
	{
		temp+=v[i];
	}
	/*temp =sum (v [0: s[0] ]*/
	temp=( temp + z[ S[0] ] ) / (S[0] +1);
	for (i=0;i<=S[0]; i++)
	{
		x[i]=temp;
	}

	/*process the middle blocks
	  If numS=1, it belongs the last block
	*/
	for (j=1; j < numS; j++)
	{
		temp=0;
		for (i= S[j-1] +1; i<= S[j]; i++)
		{
			temp+=v[i];
		}

		/*temp =sum (v [ S[j-1] +1: s[j] ]*/

		temp=(temp - z[ S[j-1] ] + z[ S[j] ])/ (S[j]- S[j-1]);

		for (i= S[j-1] +1; i<= S[j]; i++)
		{
			x[i]=temp;
		}
	}

	/*process the last block
	j=numS-1;
	*/
	temp=0;
	for (i=S[numS-1] +1 ;i< n; i++)
	{
		temp+=v[i];
	}
	/*temp =sum (v [  (S[numS-1] +1): (n-1) ]*/

	temp=( temp - z[ S[numS-1] ] ) / (nn - S[numS-1]); /*S[numS-1] <= nn-1*/

	for (i=S[numS-1] +1 ;i< n; i++)
	{
		x[i]=temp;
	}
	return numS;
}

double 
LeastCFGLasso::dualityGap2(double *z, double *g, double *s, double *Av, double lambda, int nn)
{
	int i, m;
	double temp;

	/*
	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++){
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}	
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];
    */
	
	for (i=0;i<nn;i++)
	{
		if (g[i]>0)
		{
			s[i]=lambda + z[i];
		}
		else
		{
			s[i]=-lambda + z[i];
		}
	}	

	temp=0;
	/* :'(
	 *
	m=nn%5;
	
	if (m!=0){
		for(i=0;i<m;i++)
			temp+=s[i]*g[i];
	}
	
	for(i=m;i<nn;i+=5)
		temp=temp + s[i]  *g[i]
		          + s[i+1]*g[i+1]
		          + s[i+2]*g[i+2]
		          + s[i+3]*g[i+3]
		          + s[i+4]*g[i+4];
	*/
	for(i=0;i<nn;i++)
	{
		temp += s[i]*g[i];
	}
	return temp;
	//*gap=temp;
}


int 
LeastCFGLasso::generateSolution(double *x, double *z, double *gap, double *v, double *Av, double *g, double *s, int *S, double lambda, int nn)
{
/*
  generate the solution x based on the information of z and g 
  (!!!!we assume that g has been computed as the gradient of z!!!!)

*/

	int i, m, numS, n=nn+1;
	double temp, funVal1, funVal2;
	    
	/*
	z is the appropriate solution,
	and g contains its gradient
	*/


   /*
		We assume that n>=3, and thus nn>=2
		
		  We have two ways for recovering x. 
		  The first way is x = v - A^T z
		  The second way is x =omega(z)
  */

	temp=0;
	/*
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=z[i]*(g[i] + Av[i]);
    }
    for (i=m;i<nn;i+=5)
        temp=temp + z[i]  *(g[i]   + Av[i])
                  + z[i+1]*(g[i+1] + Av[i+1])
                  + z[i+2]*(g[i+2] + Av[i+2])
                  + z[i+3]*(g[i+3] + Av[i+3])
                  + z[i+4]*(g[i+4] + Av[i+4]);
	*/
	for (i=0;i<nn;i++)
	{
		temp += z[i]*(g[i]+Av[i]);
	}
    funVal1=temp/2;
    
    temp=0;
	/*
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=fabs(g[i]);
    }
    for (i=m;i<nn;i+=5)
        temp=temp + fabs(g[i])
        + fabs(g[i+1])
        + fabs(g[i+2])
        + fabs(g[i+3])
        + fabs(g[i+4]);
	*/
	for (i=0;i<nn;i++)
	{
		temp += fabs(g[i]);
	}
    funVal1 = funVal1 + temp*lambda;
    
    /*
           we compute the solution by the second way
    */

    numS= supportSet(x, v, z, g, S, lambda, nn);
    
	/*
        we compute the objective function of x computed in the second way
    */
    
    temp=0;
	/*
    m=n%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=(x[i]-v[i]) * (x[i]-v[i]);
    }
    for (i=m;i<n;i+=5)
        temp=temp + (x[i]-  v[i]) * (  x[i]-  v[i])
                  + (x[i+1]-v[i+1]) * (x[i+1]-v[i+1])
                  + (x[i+2]-v[i+2]) * (x[i+2]-v[i+2])
                  + (x[i+3]-v[i+3]) * (x[i+3]-v[i+3])
                  + (x[i+4]-v[i+4]) * (x[i+4]-v[i+4]);
	*/
	for (i=0;i<nn;i++)
	{
		temp += ( (x[i] - v[i]) * (x[i] - v[i]) );
	}
    funVal2=temp/2;
    
    temp=0;
	/*
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=fabs( x[i+1]-x[i] );
    }
    for (i=m;i<nn;i+=5)
        temp=temp + fabs( x[i+1]-x[i] )
                  + fabs( x[i+2]-x[i+1] )
                  + fabs( x[i+3]-x[i+2] )
                  + fabs( x[i+4]-x[i+3] )
                  + fabs( x[i+5]-x[i+4] );
	*/
	for (i=0;i<nn;i++)
	{
		temp += fabs(x[i+1]-x[i]);
	}
    funVal2 = funVal2 + lambda*temp;

    if (funVal2 > funVal1)
	{
		/*
		  we compute the solution by the first way
		*/
        x[0]=v[0] + z[0];
        for(i=1;i<n-1;i++)
		{
            x[i]= v[i] - z[i-1] + z[i];
		}
        x[n-1]=v[n-1] - z[n-2];
    }
    else
	{    
        /*
        the solution x is computed in the second way
        the gap can be further reduced
        (note that, there might be numerical error)
         */
        
        *gap = *gap - (funVal1- funVal2);
        if (*gap <0)
		{
            *gap=0;
		}
	}
	return (numS);
}

int 
LeastCFGLasso::sfa_one(	double *x,     double *gap, int * activeS,
						double *z,     double * v,   double * Av, 
						double lambda, int nn,       int maxStep,
						double *s,     double *g,
						double tol,    int tau)
{
	int i, iterStep, m, tFlag=0, n=nn+1;
	double temp;
	//int* S=(int *) malloc(sizeof(int)*nn);
	int *S = new int[nn];
	double gapp=-1, gappp=-2;	/*gapp denotes the previous gap*/
	int numS=-100, numSp=-200, numSpp=-300;    
	/*
	numS denotes the number of elements in the Support Set S
	numSp denotes the number of elements in the previous Support Set S
	*/

	*gap=-1; /*initialize *gap a value*/

	/*
	The main algorithm by Nesterov's method

     B is an nn x nn tridiagonal matrix.

     The nn eigenvalues of B are 2- 2 cos (i * PI/ n), i=1, 2, ..., nn
	*/

	/*
	we first do a gradient step based on z
	*/


	/*
		---------------------------------------------------
		  A gradient step  begins
	*/
		g[0]=z[0] + z[0] - z[1] - Av[0];
		for (i=1;i<nn-1;i++)
		{
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];

		
		/* 
		do a gradient step based on z to get the new z
		*/
		/*
		 * Again, why?!
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=z[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = z[i]   -  g[i]  /4;
			z[i+1] = z[i+1] -  g[i+1]/4;
			z[i+2] = z[i+2] -  g[i+2]/4;
			z[i+3] = z[i+3] -  g[i+3]/4;
			z[i+4] = z[i+4] -  g[i+4]/4;
		}
		*/
		for (i=0;i<nn;i++)
		{
			z[i] = z[i] - g[i]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda

        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++)
		{
			if (z[i]>lambda)
			{
				z[i]=lambda;
			}
			else
			{
				if (z[i]<-lambda)
				{
					z[i]=-lambda;
				}
			}
		}

		/*
		---------------------------------------------------
		  A gradient descent step ends
		*/


	/*compute the gradient at z*/

	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++)
	{
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}	
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];

	for (iterStep=1; iterStep<=maxStep; iterStep++)
	{
		/*
		---------------------------------------------------
		restart the algorithm with x=omega(z)
		*/
				
		numSpp=numSp;
		numSp=numS; /*record the previous numS*/
		numS = supportSet(x, v, z, g, S, lambda, nn);
		
		/*With x, we compute z via
		AA^T z = Av - Ax
		 */
	
		/*
		compute s= Av -Ax
		*/
		
		for (i=0;i<nn; i++)
		{
			s[i]=Av[i] - x[i+1] + x[i];
		}
		/*
		Apply Rose Algorithm for solving z
		*/
		temp = Thomas(z, s, nn);
	    /*
	    project z to [-lambda, lambda]
	    */
		for(i=0;i<nn;i++)
		{
			if (z[i]>lambda)
			{
				z[i]=lambda;
			}
			else
			{
				if (z[i]<-lambda)
				{
					z[i]=-lambda;
				}
			}
		}
		/*
		---------------------------------------------------
		restart the algorithm with x=omega(z)
        we have computed a new z, based on the above relationship
		*/

		/*
		---------------------------------------------------
		  A gradient step  begins
		*/
		g[0]=z[0] + z[0] - z[1] - Av[0];
		for (i=1;i<nn-1;i++)
		{
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];

		/* 
		do a gradient step based on z to get the new z
		*/
		/* -_-
		 *
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=z[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = z[i]   -  g[i]  /4;
			z[i+1] = z[i+1] -  g[i+1]/4;
			z[i+2] = z[i+2] -  g[i+2]/4;
			z[i+3] = z[i+3] -  g[i+3]/4;
			z[i+4] = z[i+4] -  g[i+4]/4;
		}
		*/
		for (i=0;i<nn;i++)
		{
			z[i] = z[i] - g[i]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda

        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++)
		{
			if (z[i]>lambda)
			{
				z[i]=lambda;
			}
			else
			{
				if (z[i]<-lambda)
				{
					z[i]=-lambda;
				}
			}
		}

		/*
		---------------------------------------------------
		  A gradient descent step ends
		*/

		/*compute the gradient at z*/

		g[0]=z[0] + z[0] - z[1] - Av[0];
		for (i=1;i<nn-1;i++)
		{
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}	
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];


		if (iterStep % tau==0)
		{
			gappp=gapp;
			gapp=*gap;  /*record the previous gap*/

			(*gap) = dualityGap2(z, g, s, Av, lambda, nn);
			/*g, the gradient of z should be computed before calling this function*/

			if (*gap <=tol)
			{
				tFlag=1;
				break;
			}

			m=1;
			if (nn > 1000000)
			{
				m=5;
			}
			else
			{
				if (nn > 100000)
				{
					m=3;
				}
			}

			if ( abs( numS-numSp) <m )
			{
				m=generateSolution(x, z, gap, v, Av, g, s, S, lambda, nn);
				/*g, the gradient of z should be computed before calling this function*/

				if (*gap < tol)
				{
					numS=m;
					tFlag=2;					
					break;
				}

				if ( (*gap ==gappp) && (numS==numSpp) )
				{
					tFlag=2;
					break;
							
				}
				
            } /*end of if*/

		}/*end of if tau*/
	} /*end of for*/

	if (tFlag!=2)
	{
		numS=generateSolution(x, z, gap, v, Av, g, s, S, lambda, nn);
       /*g, the gradient of z should be computed before calling this function*/
    }

	//free(S);
	delete [] S;

	*activeS=numS;
	return iterStep;
}


double
LeastCFGLasso::Thomas(double* z0, double* Av, int nn)
{
	//
	// We are supposed to solve
	// A A^T z0 = Av
	// Right now, it only works if A is linear tree 0->1->2->3 -_-
	// So, if A is
	// [-1     1     0     0     0
	//   0    -1     1     0     0
	//   0     0    -1     1     0
	//   0     0     0    -1     1 ]
	// B = A A^T is
	// [ 2    -1     0     0
	//  -1     2    -1     0
	//   0    -1     2    -1
	//   0     0    -1     2] 
	// I guess if we want to fix this, we should use a different solution than Thomas :/
	//

	int i;
	double tt, z_max;

	z0[0] = Av[0]/2;
	for (i=1;i < nn; i++)
	{
		tt = Av[i] + z0[i-1];
		z0[i] = tt - tt/(i+2);
	}

	z_max = fabs(z0[nn-1]);

	for (i=nn-2; i>=0; i--)
	{
		z0[i] += z0[i+1] - z0[i+1]/(i+2);
		tt = fabs(z0[i]);
		if (tt > z_max)
		{
			z_max=tt;
		}
	}
	return z_max;
}

//w_1 = flsa(v, w0,  lambda_1, lambda_2, length(v), 1000, 1e-9, 1, 6);
int
LeastCFGLasso::flsa(double* x, double* v, double* z0, double lambda1, double lambda2, int n, int maxStep, double tol, int tau, int flag)
{
	//see https://github.com/jiayuzhou/MALSAR/blob/master/MALSAR/c_files/flsa/flsa.h
	
	//double* x = new double[n];
	double* z = new double[n-1];
	//double* infor = new double[4];

	int i, nn=n-1, m;
	double zMax, temp;
	double *Av, *g, *s;
	int iterStep, numS;
	double gap;

	//Av=(double *) malloc(sizeof(double)*nn);
    Av = new double[nn];

	/*
	Compute Av= A*v(n=4, nn=3)
     A= [ -1   1   0   0
           0  -1   1   0
	       0   0  -1   1 ]
	*/
	//It assumes a linear tree
	for (i=0;i<nn; i++)
	{
		Av[i]=v[i+1]-v[i];
	}

	/*
	Sovlve the linear system via Thomas's algorithm or Rose's algorithm
     B * z0 = Av
	*/
	// B = is A A^T
	// Again, assumes a linear tree
    zMax = Thomas(z, Av, nn);

	/*
	We consider two cases:
	   1) lambda2 >= zMax, which leads to a solution with same entry values
	   2) lambda2 < zMax, which needs to first run sfa, and then perform soft thresholding
	*/

	/*
	First case: lambda2 >= zMax
	*/
	if (lambda2 >= zMax)
	{
		temp=0;
		/* I'm not sure why they were doing it like this
		 * Speed?!
		m = n%5; //Why 5?! :/
		if (m!=0)
		{
			for (i=0;i<m;i++)
			{
				temp+=v[i];
			}
		}		
		for (i=m;i<n;i+=5)
		{
			temp += v[i] + v[i+1] + v[i+2] + v[i+3] + v[i+4];
		}
		temp/=n; 
		*/
		/* temp is the mean value of v*/
		for (i=0;i<n;i++)
		{
			temp += v[i];
		}
		temp/=n; 

		/*
		soft thresholding by lambda1
		*/
		if (temp> lambda1)
		{
			temp= temp-lambda1;
		}
		else
		{
			if (temp < -lambda1)
			{
				temp= temp+lambda1;
			}
			else
			{
				temp=0;
			}
		}

		/* Again, not sure why they were doing it like this :/
		m=n%7;
		if (m!=0){
			for (i=0;i<m;i++)
				x[i]=temp;
		}
		for (i=m;i<n;i+=7){
			x[i]   =temp;
			x[i+1] =temp;
			x[i+2] =temp;
			x[i+3] =temp;
			x[i+4] =temp;
			x[i+5] =temp;
			x[i+6] =temp;
		}
		*/
		for (i=0;i<n;i++)
		{
			x[i] = temp;
		}
		
		gap=0;

		delete[] Av;
		//free(Av);

		//We don't need the info
		//infor[0]= gap;
		//infor[1]= 0;
		//infor[2]=zMax;
		//infor[3]=0;

		return 0;
	}

	/*
	Second case: lambda2 < zMax

    We need to call sfa for computing x, and then do soft thresholding

    Before calling sfa, we need to allocate memory for g and s, 
	           and initialize z and z0.
	*/


	/*
	Allocate memory for g and s
	*/

	//g    =(double *) malloc(sizeof(double)*nn),
	//s    =(double *) malloc(sizeof(double)*nn);
	g = new double [nn],
	s = new double [nn];

	m = flag /10;
	// flag was 6 in FGLasso_projection_rowise, so I don't go into else
	//if (m==0)
	//{
	for (i=0;i<nn;i++)
	{
		if (z0[i] > lambda2)
		{
			z[i]=lambda2;
		}
		else
		{
			if (z0[i]<-lambda2)
			{
				z[i]=-lambda2;
			}
			else
			{
				z[i]=z0[i];
			}
		}
	}
	//}
	
	//flag was 6, I go to sfa_one
	//flag=flag %10;
	
	iterStep=sfa_one(x, &gap, &numS, z, v, Av, lambda2, nn, maxStep, s, g, tol, tau);
	/*
	soft thresholding by lambda1
	*/

	for(i=0;i<n;i++)
	{
		if (x[i] > lambda1)
		{
			x[i]-=lambda1;
		}
		else
		{
			if (x[i]<-lambda1)
			{
				x[i]+=lambda1;
			}
			else
			{
				x[i]=0;
			}
		}
	}
	//free(Av);
	//free(g);
	//free(s);
	delete [] Av;
	delete [] g;
	delete [] s;

	delete [] z;

	//we don't need infor
	//infor[0]=gap;
	//infor[1]=iterStep;
	//infor[2]=zMax;
	//infor[3]=numS;
	
	return 0;
}

double 
LeastCFGLasso::getNonSmooth(vector<Matrix*> allW,map<int,vector<int>>& tree)
{
    // for i = 1 : size(W, 1)
    //    w = W(i, :);
    //    non_smooth_value = non_smooth_value ...
    //        + rho_1 * norm(w, 1) + rho_2 * norm(R * w', 1) ...
    //        + rho_3 * norm(w, 2);
    // end
	double res=0;
	int row = allW[0]->getRowCnt();
	for(int j=0;j<row;j++)
	{
		double wn=0;//norm2
		double w1=0;//norm1
		double wr=0;//pairwise
		for(int i=0;i<allW.size();i++)
		{
			Matrix* W = allW[i];
			double  v = W->getValue(j,0);
			wn += (v*v);
			w1 += fabs(v);
		}
		wn = sqrt(wn);
		for (map<int,vector<int>>::iterator itr=tree.begin();itr!=tree.end();itr++)
		{
			int t1 = itr->first;
			Matrix* W1 = allW[t1];
			double  v1 = W1->getValue(j,0);

			vector<int>& es = itr->second;
			for(int i=0;i<es.size();i++)
			{
				int t2 = es[i];
				Matrix* W2 = allW[t2];
				double  v2 = W2->getValue(j,0);

				wr += fabs(v1-v2);
			}
		}
		res = res + rho1*w1 + rho2*wr + rho3*wn;
	}
	return res;
}

int
LeastCFGLasso::getGrad(vector<Task_T*>* allt, vector<Matrix*> allXT, vector<Matrix*> allXY, vector<Matrix*> allW, vector<Matrix*>& allGW)
{
	// grad_W = [];
	// for i = 1:task_num
	//   grad_W = cat(2, grad_W, X{i}*(X{i}' * W(:,i)-Y{i}) );
	// end
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
		allGW.push_back(GW);
		delete XW;
		delete R;
		//delete XY;
		//delete XT;
	}
	return 0;
}

double
LeastCFGLasso::lossFunc(vector<Task_T*>* allt, vector<Matrix*> allXT, vector<Matrix*> allW)
{
	// for i = 1: task_num
	//	funcVal = funcVal + 0.5 * norm (Y{i} - X{i}' * W(:, i))^2;
	// end
	double funcVal = 0;
	for (int i=0;i<allt->size();i++)
	{
		Task_T* t  = allt->at(i);
		Matrix* W  = allW[i];
		Matrix* Y  = t->Y;
		//Matrix* XT = X->transMatrix();
		Matrix* XT = allXT[i];
		Matrix* R  = XT->multiplyMatrix(W);
		Matrix* S  = Y->subtractMatrix(R);
		double  f  = S->getFNorm();
		funcVal = funcVal + (0.5*f*f);
		//delete XT;
		delete R;
		delete S;
	}
	return funcVal;
}

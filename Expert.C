/*
CMINT: An algorithm to cluster functional genomics data from multiple celltypes
    Copyright (C) 2016 Sushmita Roy sushroy@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <math.h>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
//#include "Evidence.H"
#include "DRMNPotential.H"
#include "Expert.H"


Expert::Expert()
{
	covariance=NULL;
	mean=NULL;
	invCovariance=NULL;
	dPot=NULL;
}

Expert::~Expert()
{
	if(covariance!=NULL)
	{
		delete covariance;
	}
	if(mean!=NULL)
	{
		delete mean;
	}
	if(invCovariance!=NULL)
	{
		delete invCovariance;
	}
	if(dPot!=NULL)
	{
		delete dPot;
	}
}


int 
Expert::setMean(Matrix* m)
{
	if(mean!=NULL)
	{
		delete mean;
	}
	mean=m;
	return 0;
}

int 
Expert::setCovariance(Matrix* c)
{
	if(covariance!=NULL)
	{
		delete covariance;
	}
	covariance=c;
	//c->showMatrix();
	invCovariance=covariance->invMatrix();
	
	double det=1;
	double n=((double)c->getRowCnt());
	for(int i=0;i<n;i++)
	{
		//det=det+log(c->getValue(i,i));
		det=det*(c->getValue(i,i));
	}
	
	normFactor=pow(2*PI,n)*det;
	//normFactor=(n*log(2*PI)) + det;
	normFactor=sqrt(normFactor);
	//normFactor=normFactor;
	return 0;
}



int
Expert::updateCovariance()
{
	if(invCovariance!=NULL)
	{
		delete invCovariance;
	}
	invCovariance=covariance->invMatrix();
	/*double det=covariance->detMatrix();
	double n=((double)covariance->getRowCnt());
	normFactor=pow(2*PI,n)*det;
	normFactor=sqrt(normFactor);*/
	double det=1;
	double n=((double)covariance->getRowCnt());
	for(int i=0;i<n;i++)
	{
		//det=det+log(covariance->getValue(i,i));
		det=det*covariance->getValue(i,i);
	}
	
	normFactor=pow(2*PI,n)*det;
	//normFactor=(n*log(2*PI)) + det;
	normFactor=sqrt(normFactor);
	//normFactor=normFactor/2;
	return 0;
}


double 
Expert::getOutputPDF(vector<double>* y)
{
	double pdf=0;
	Matrix diffMat(1,y->size());
	int colCnt=y->size();
	for(int i=0;i<colCnt;i++)
	{
		double diff=(*y)[i]-mean->getValue(0,i);
		diffMat.setValue(diff,0,i);
	}
	Matrix* p1=diffMat.multiplyMatrix(invCovariance);
	double sum=0;
	for(int i=0;i<colCnt;i++)
	{
		double s=p1->getValue(0,i)*diffMat.getValue(0,i);
		if(isinf(s))
		{
			cout <<"Nan: at "<< i << " diffmat " << diffMat.getValue(0,i) <<endl; 
		}
		sum=sum+s;
	}
	pdf=exp(-0.5*sum)/normFactor;
	if(isnan(pdf))
	{
		cout << "PDF is nan for sum=" <<sum << " normfactor=" << normFactor << endl;
	}
	//pdf=(-0.5*sum)-normFactor;
	//pdf=exp(pdf);
	if(pdf<1e-80)
	{
		if(pdf<minpdf)
		{	
			minpdf=pdf;
		}
		pdf=1e-80;
		clipCnt++;
	}
	delete p1;
	return pdf;
}


double 
Expert::getOutputPDF_Nocov(vector<double>* y)
{
	double pdf=0;
	Matrix diffMat(1,y->size());
	int colCnt=y->size();
	double lpdf=0;
	for(int i=0;i<colCnt;i++)
	{
		double mean_i=mean->getValue(0,i);
		double var_i=covariance->getValue(i,i);
		double val_i=(*y)[i];
		double gpdf=gsl_ran_gaussian_pdf(val_i-mean_i,sqrt(var_i));
		if(gpdf<1e-30)
		{
			gpdf=1e-30;
		}
		lpdf=lpdf+log(gpdf);
	}
	pdf=exp(lpdf);
	if(pdf<1e-40)
	{
		pdf=1e-40;
	}
	return pdf;
}


int 
Expert::setPrior(double p)
{
	priorProb=p;
	return 0;
}

double 
Expert::getPrior()
{
	return priorProb;
}

Matrix*
Expert::getMean()
{
	return mean;
}

Matrix*
Expert::getCovariance()
{
	return covariance;
}



int
Expert::assignGeneToExpert(const char* geneName)
{
	string key(geneName);
	geneSet[key]=0;
	return 0;
}

map<string,int>&
Expert::getGeneSet()
{
	return geneSet;
}

int
Expert::resetAssignedGenes()
{
	geneSet.clear();
}

double
Expert::getEntropy()
{
	double determinant=covariance->detMatrix();
	//cout <<"Determinant: " << determinant << endl;
	double n=((double)covariance->getRowCnt())/2.0;
	double commFact=1+log(2*3.1472);
	double jointEntropy=0.5*((n*commFact) + log(determinant));
	//cout <<"Joint entropy: " << jointEntropy << endl;
	return jointEntropy;
}

int
Expert::generateSample(gsl_rng* r, vector<double>& sample)
{
	int dim=mean->getColCnt();
	for(int d=0;d<dim;d++)
	{
		double meanval=mean->getValue(0,d);
		double sdval=covariance->getValue(d,d);
		double aval=gsl_ran_gaussian(r,sqrt(sdval));
		aval=aval+meanval;
		sample.push_back(aval);
	}
	return 0;
}


int 
Expert::resetClip()
{
	clipCnt=0;
	return 0;
}

int 
Expert::getClip()
{
	return clipCnt;
}

int
Expert::addRegulator(string& reg)
{
	regulators[reg]=0;
}

map<string,int>& 
Expert::getCurrentRegSet()
{
	return regulators;
}

int 
Expert::setLLScore(double score)
{
	llScore=score;
	return 0;
}

double 
Expert::getLLScore()
{
	return llScore;
}

int
Expert::setDRMNPotential(DRMNPotential* ptr)
{
	dPot=ptr;
	return 0;
}

DRMNPotential*
Expert::getDRMNPotential()
{
	return dPot;
}

double
//Expert::getDRMNProb(map<int,Evidence*>* evidSet)
Expert::getDRMNProb(map<int,double>* evidSet)
{
	double prob=dPot->getCondPotValueFor(evidSet);
	return prob;
}

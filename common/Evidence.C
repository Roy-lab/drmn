#include "Evidence.H"

Evidence::Evidence()
{
}

Evidence::~Evidence()
{
}

int
Evidence::assocVariable(int id)
{
	varID=id;
	return 0;
}

int
Evidence::getAssocVariable()
{
	return varID;
}

int 
Evidence::setEvidVal(double aval)
{
	contValue=aval;
	return 0;
}

double 
Evidence::getEvidVal()
{
	return contValue;
}

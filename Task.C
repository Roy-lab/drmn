#include "Task.H"

Task_T::Task_T()
{
	X=NULL;
	Y=NULL;
	YNames=NULL;
}

Task_T::~Task_T()
{
	if (X!=NULL)
	{
		delete X;
	}
	if (Y!=NULL)
	{
		delete Y;
	}
	if (YNames!=NULL)
	{
		YNames->clear();
		delete YNames;
	}
}

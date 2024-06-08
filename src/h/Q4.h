# pragma once

# include "Element.h"

using namespace std;

//! Q4 element class
class CQ4 : public CElement
{
public:

//! Constructor
	CQ4();

//! Deconstructor
	~CQ4();

//! Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//! Write element data to stream
	virtual void Write(COutputter& output);

//! Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//! Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};

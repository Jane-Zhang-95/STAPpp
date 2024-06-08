#pragma once

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CGauss
{
public:

	unsigned int ngp;	//!< Number of gauss point

public:

//! Virtual deconstructor
    virtual ~CGauss() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//for T3
class CT3Gauss : public CGauss
{
public:
	double* csi; //gauss point
    double* w;//weight

public:
//!	Constructor
	CT3Gauss();

//!	Desconstructor
    ~CT3Gauss();

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);

//
    virtual void T3_GPoint();

};


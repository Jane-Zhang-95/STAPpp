/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//! Material class for 2D plate element
class C2DMaterial : public CMaterial
{
public:

	double mu, thickness;	//!< Poisson ratio, Thickness

	bool plane_stress;	//!< Plane stress = True, Plane strain = False 

	double* D;

	double Area;

public:

	C2DMaterial(): D(nullptr){}

	~C2DMaterial(){
		if(D)
			delete[] D;
	}

//! Read material data from stream Input
	virtual bool Read(ifstream& Input);

//! Write material data to Stream
	virtual void Write(COutputter& output);
};

//! Material class for 3D plate element
class C3DMaterial : public CMaterial
{
public:

	double mu;	//!< Poisson ratio, Thickness

	double* D;

	double Area;

public:

	C3DMaterial(): D(nullptr){}

	~C3DMaterial(){
		if(D)
			delete[] D;
	}

//! Read material data from stream Input
	virtual bool Read(ifstream& Input);

//! Write material data to Stream
	virtual void Write(COutputter& output);
};
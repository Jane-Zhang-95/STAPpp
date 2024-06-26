/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include <ctime>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
	
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();
	
	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT                                      BODY" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES                                      FORCE" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, Q4 ELEMENTS" << endl
			  << "     EQ.3, T3 ELEMENTS" << endl
			  << "     EQ.4, H8 ELEMENTS" << endl
			  << "     EQ.5, Tet4 ELEMENTS" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::T3: // T3 element
				OutputT3Elements(EleGrp);
				break;
			case ElementTypes::Q4: // Q4 element
				OutputQ4Elements(EleGrp);
				break;
			case ElementTypes::H8: // H8 element
				OutputH8Elements(EleGrp);
				break;	
			case ElementTypes::Tet4: // Q4 element
				OutputTet4Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}
//	Output bar element data
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	Output T3 element data
void COutputter::OutputT3Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S         POISSON          THICK            PLANE" << endl
		<< " NUMBER     MODULUS          RATIO           -NESS           STRESS" << endl
		<< "               E               MU              T              FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      A        B        C        SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

void COutputter::OutputQ4Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND PLATE MATERIAL CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S         POISSON          THICK            PLANE" << endl
		<< " NUMBER     MODULUS          RATIO           -NESS           STRESS" << endl
		<< "               E               MU              T              FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE     NODE     NODE     NODE       MATERIAL" << endl
		<< " NUMBER-N      A        B        C        D       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}

void COutputter::OutputTet4Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND PLATE MATERIAL CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S         POISSON" << endl
		<< " NUMBER     MODULUS          RATIO" << endl
		<< "               E               MU" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE     NODE     NODE     NODE       MATERIAL" << endl
		<< " NUMBER-N      A        B        C        D       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}


//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S         POISSON          THICK            PLANE" << endl
		<< " NUMBER     MODULUS          RATIO           -NESS           STRESS" << endl
		<< "               E               MU              T              FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      A        B        C        SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

void COutputter::OutputH8Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND PLATE MATERIAL CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S         POISSON" << endl
		<< " NUMBER     MODULUS          RATIO" << endl
		<< "               E               MU" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE     NODE     NODE     NODE      NODE      NODE     NODE      NODE     MATERIAL" << endl
		<< " NUMBER-N      A        B        C        D         E         F        G         H     SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl;
		*this << "     NUMBER OF NATURAL BC . =" << setw(6) << LoadData->nnbc << endl
			  << endl;
		*this << "    NODE          DIRECTION         LOAD            ELEMENT" << endl
			  << "   NUMBER                         MAGNITUDE         NUMBER" << endl;

		LoadData->Write(*this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()  
{  
	CDomain* FEMData = CDomain::GetInstance();  
	  
	CNode* NodeList = FEMData->GetNodeList();  
	  
	double* Displacement = FEMData->GetDisplacement();  
  
	*this << setiosflags(ios::scientific);  
  
	*this << " D I S P L A C E M E N T S" << endl  
		  << endl;  
	   
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;  
  
	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)  
	{   
		NodeList[np].WriteNodalDisplacement(*this, Displacement);  
	}  
  
	*this << endl;  
}

//	Calculate stresses
void COutputter::OutputElementStress()  
{  
	CDomain* FEMData = CDomain::GetInstance();  
  
	double* Displacement = FEMData->GetDisplacement();  
  
	unsigned int NUMEG = FEMData->GetNUMEG();  
  
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)  
	{  
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)  
			  << EleGrpIndex + 1 << endl  
			  << endl;  
  
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];  
  
		unsigned int NUME = EleGrp.GetNUME();  
  
		ElementTypes ElementType = EleGrp.GetElementType();  
		
		double stress_Bar = 0;
		double* stress_T3 = new double[3];
  
		switch (ElementType)  
		{  
			case ElementTypes::Bar:  
				*this << "  ELEMENT             FORCE            STRESS" << endl  
					<< "  NUMBER" << endl;   
  
				// for all element 
				for (unsigned int Ele = 0; Ele < NUME; Ele++)  
				{  
					// get element now
					CElement& Element = EleGrp[Ele];  
  
					Element.ElementStress(&stress_Bar, Displacement);  
  
					// to CBarmaterial  
					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());  
   
					*this << setw(5) << Ele + 1 << setw(22) << stress_Bar * material.Area << setw(18)  
						<< stress_Bar << endl;  
				}  
  
				*this << endl;  
  
				break;  

			case ElementTypes::T3: 
				*this << "  ELEMENT            STRESS" << endl  
					<< "  NUMBER" << endl;  
  
				for (unsigned int Ele = 0; Ele < NUME; Ele++)  
				{   
					CElement& Element = EleGrp[Ele];  
  
					Element.ElementStress(stress_T3, Displacement);  
  
					//CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());  
    
					*this << setw(5) << Ele + 1 << setw(22) << stress_T3[0] << setw(18)
						 << stress_T3[1] << setw(18) << stress_T3[2] << endl;  
				}  
  
				*this << endl; 
  
				break; 

			case ElementTypes::Q4: // Q4 element
				*this << "  ELEMENT(GAUSS)                                      STRESS COMPONENTS" << endl
					  << "     NUMBER          S11             S22             S33             S12             S13             S23" << endl;
				
				double Stress[24];
        
//	Print nodal displacement
void COutputter::OutputNodalDisplacement()  
{  
	CDomain* FEMData = CDomain::GetInstance();  
	  
	CNode* NodeList = FEMData->GetNodeList();  
	  
	double* Displacement = FEMData->GetDisplacement();  
  
	*this << setiosflags(ios::scientific);  
  
	*this << " D I S P L A C E M E N T S" << endl  
		  << endl;  
	   
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;  
  
	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)  
	{   
		NodeList[np].WriteNodalDisplacement(*this, Displacement);  
	}  
  
	*this << endl;  
}

//	Calculate stresses
void COutputter::OutputElementStress()  
{  
	CDomain* FEMData = CDomain::GetInstance();  
  
	double* Displacement = FEMData->GetDisplacement();  
  
	unsigned int NUMEG = FEMData->GetNUMEG();  
  
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)  
	{  
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)  
			  << EleGrpIndex + 1 << endl  
			  << endl;  
  
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];  
  
		unsigned int NUME = EleGrp.GetNUME();  
  
		ElementTypes ElementType = EleGrp.GetElementType();  
		
		double stress_Bar = 0;
		double* stress_T3 = new double[3];
		double* L2_error = new double[1];
		L2_error[0] = 0;
  
		switch (ElementType)  
		{  
			case ElementTypes::Bar:  
				*this << "  ELEMENT             FORCE            STRESS" << endl  
					<< "  NUMBER" << endl;   
  
				// for all element 
				for (unsigned int Ele = 0; Ele < NUME; Ele++)  
				{  
					// get element now
					CElement& Element = EleGrp[Ele];  
  
					Element.ElementStress(&stress_Bar, Displacement);  
  
					// to CBarmaterial  
					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());  
   
					*this << setw(5) << Ele + 1 << setw(22) << stress_Bar * material.Area << setw(18)  
						<< stress_Bar << endl;  
				}  
  
				*this << endl;  
  
				break;  

			case ElementTypes::T3: 
				*this << "  ELEMENT                          STRESS" << endl  
					<< "  NUMBER" << endl;  
  

				for (unsigned int Ele = 0; Ele < NUME; Ele++)  
				{   
					CElement& Element = EleGrp[Ele];  
  
					Element.ElementStress(stress_T3, Displacement);

					CT3& T3Element = dynamic_cast<CT3&>(Element);

					T3Element.GetError(stress_T3, Displacement, L2_error);
  
					//CT3Material& material = *dynamic_cast<CT3Material*>(Element.GetElementMaterial());  
    
					*this << setw(5) << Ele + 1 << setw(22) << stress_T3[0] << setw(18)
						 << stress_T3[1] << setw(18) << stress_T3[2] << endl;  
				}  
				std::cout<<endl;
				std::cout<< "L2_error is: "<<L2_error[0] << endl;
				*this  << endl; 
  
				break; 

			case ElementTypes::Q4: // Q4 element
				*this << "  ELEMENT(GAUSS)                                      STRESS COMPONENTS" << endl
					  << "     NUMBER          S11             S22             S33             S12             S13             S23" << endl;
				
				double Stress[24];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(Stress, Displacement);

					C2DMaterial& material = *dynamic_cast<C2DMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << " (-,-)" << setw(18) << Stress[0] << setw(16) << Stress[1] << setw(16) << Stress[2]
						<< setw(16) << Stress[3] << setw(16) << Stress[4] << setw(16) << Stress[5] << endl;
					*this << setw(5) << Ele + 1 << " (-,+)" << setw(18) << Stress[6] << setw(16) << Stress[7] << setw(16) << Stress[8]
						<< setw(16) << Stress[9] << setw(16) << Stress[10] << setw(16) << Stress[11] << endl;
					*this << setw(5) << Ele + 1 << " (+,-)" << setw(18) << Stress[12] << setw(16) << Stress[13] << setw(16) << Stress[14]
						<< setw(16) << Stress[15] << setw(16) << Stress[16] << setw(16) << Stress[17] << endl;
					*this << setw(5) << Ele + 1 << " (+,+)" << setw(18) << Stress[18] << setw(16) << Stress[19] << setw(16) << Stress[20]
						<< setw(16) << Stress[21] << setw(16) << Stress[22] << setw(16) << Stress[23] << endl;
				}

				*this << endl;

				break;

			case ElementTypes::Tet4: // Tet4 element
				*this << "  ELEMENT                                             STRESS COMPONENTS" << endl
					  << "  NUMBER          S11             S22             S33             S12             S13             S23" << endl;
				
				double stress[6];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
          				Element.ElementStress(stress, Displacement);

					C3DMaterial& material = *dynamic_cast<C3DMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(18) << stress[0] << setw(16) << stress[1] << setw(16) << stress[2]
						<< setw(16) << stress[3] << setw(16) << stress[4] << setw(16) << stress[5] << endl;
        			}
        			*this << endl;
				break;

			case ElementTypes::H8: // H8 element
				*this << "  ELEMENT(GAUSS)                                      STRESS COMPONENTS" << endl
					  << "     NUMBER          S11             S22             S33             S12             S13             S23" << endl;
				
				double Stress_H8[48];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(Stress_H8, Displacement);

					C3DMaterial& material = *dynamic_cast<C3DMaterial*>(Element.GetElementMaterial());
					
					*this << setw(5) << Ele + 1 << " (-,-,-)" << setw(18) << Stress_H8[0] << setw(16) << Stress_H8[1] << setw(16) << Stress_H8[2]
						<< setw(16) << Stress_H8[3] << setw(16) << Stress_H8[4] << setw(16) << Stress_H8[5] << endl;
					*this << setw(5) << Ele + 1 << " (-,-,+)" << setw(18) << Stress_H8[6] << setw(16) << Stress_H8[7] << setw(16) << Stress_H8[8]
						<< setw(16) << Stress_H8[9] << setw(16) << Stress_H8[10] << setw(16) << Stress_H8[11] << endl;
					*this << setw(5) << Ele + 1 << " (-,+,-)" << setw(18) << Stress_H8[12] << setw(16) << Stress_H8[13] << setw(16) << Stress_H8[14]
						<< setw(16) << Stress_H8[15] << setw(16) << Stress_H8[16] << setw(16) << Stress_H8[17] << endl;
					*this << setw(5) << Ele + 1 << " (-,+,+)" << setw(18) << Stress_H8[18] << setw(16) << Stress_H8[19] << setw(16) << Stress_H8[20]
						<< setw(16) << Stress_H8[21] << setw(16) << Stress_H8[22] << setw(16) << Stress_H8[23] << endl;
					*this << setw(5) << Ele + 1 << " (+,-,-)" << setw(18) << Stress_H8[24] << setw(16) << Stress_H8[25] << setw(16) << Stress_H8[26]
						<< setw(16) << Stress_H8[27] << setw(16) << Stress_H8[28] << setw(16) << Stress_H8[29] << endl;
					*this << setw(5) << Ele + 1 << " (+,-,+)" << setw(18) << Stress_H8[30] << setw(16) << Stress_H8[31] << setw(16) << Stress_H8[32]
						<< setw(16) << Stress_H8[33] << setw(16) << Stress_H8[34] << setw(16) << Stress_H8[35] << endl;
					*this << setw(5) << Ele + 1 << " (+,+,-)" << setw(18) << Stress_H8[36] << setw(16) << Stress_H8[37] << setw(16) << Stress_H8[38]
						<< setw(16) << Stress_H8[39] << setw(16) << Stress_H8[40] << setw(16) << Stress_H8[41] << endl;
					*this << setw(5) << Ele + 1 << " (+,+,+)" << setw(18) << Stress_H8[42] << setw(16) << Stress_H8[43] << setw(16) << Stress_H8[44]
						<< setw(16) << Stress_H8[45] << setw(16) << Stress_H8[46] << setw(16) << Stress_H8[47] << endl;

       				 	}
          				*this << endl;
					break;

			default:   
				cerr << "*** Error *** Elment type " << ElementType  
					<< " has not been implemented.\n\n";  
		}  
		
				delete[] stress_T3; 

				delete[] L2_error;
	}  
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

// 调试
#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif

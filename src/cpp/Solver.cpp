/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>

using namespace std;

// LDLT facterization
void CLDLTSolver::LDLT()  
{  
	unsigned int N = K.dim(); 
	//std::cout<<N<<endl;
      
    unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Hights  
  
	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)  
	{  
		//std::cout<<"^^ "<<K(j,j)<<endl;
		
		unsigned int mj = j - ColumnHeights[j-1];  
           
		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)  
		{   
			unsigned int mi = i - ColumnHeights[i-1];  
  
			double C = 0.0;  
              
			for (unsigned int r = max(mi, mj); r <= i-1; r++)  
				C += K(r,i) * K(r,j);		// C += L_ri * U_rj  
  
			K(i,j) -= C;	// U_ij = K_ij - C  
		}  
  
		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)  
		{  
			double Lrj = K(r,j) / K(r,r);	// L_rj = U_rj / D_rr  
  
			K(j,j) -= Lrj * K(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			//std::cout<<K(j,j)<<endl;
			K(r,j) = Lrj;  
			//std::cout<<K(j,j)<<endl;
		}  
  
        if (fabs(K(j,j)) <= FLT_MIN)  
        {   
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl  
            	 << "    Euqation no = " << j << endl  
            	 << "    Pivot = " << K(j,j) << endl;  
              
            exit(4);  
        }  
    }  
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)  
{  
	unsigned int N = K.dim();    
    unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Heights  
    
	// Reduce right-hand-side load vector (LV = R)  
	for (unsigned int i = 2; i <= N; i++)	
	{   
        unsigned int mi = i - ColumnHeights[i-1];  
  
		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1  
			Force[i-1] -= K(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)  
	}  
  
	// Back substitute (Vbar = D^(-1) V, L^T a = Vbar)  
	for (unsigned int i = 1; i <= N; i++)	
		Force[i-1] /= K(i,i);	// Vbar = D^(-1) V  
   
	for (unsigned int j = N; j >= 2; j--)	 
	{  
        unsigned int mj = j - ColumnHeights[j-1];  
  
		for (unsigned int i = mj; i <= j-1; i++)	
			Force[i-1] -= K(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)  
	}  

};

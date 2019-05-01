/* Written by Leonardo Silvestri 2007-2013 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <string> 
#include <sstream>
#include <complex>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <string>
#include <fstream>
#include <iomanip>
#include <OPutils.h>
#include <N_Functions.h> 
#include <Transfer_matrix.h>

namespace ublas = boost::numeric::ublas;
using namespace std;

#define hc 0.197 // eV  micron
#ifndef Pi 
#define Pi 3.141592654
#endif
#ifndef ZERO_TOL
#define ZERO_TOL 1e-10  // Numbers smaller than this are considered equal to 0
#endif

typedef std::complex<double>							dComplex;
typedef ublas::vector<double>							dVector;
typedef ublas::vector<dComplex>							dComplexVector;
typedef ublas::matrix<double,	ublas::column_major>	dMatrix;
typedef ublas::matrix<dComplex,	ublas::column_major>	dComplexMatrix;

// Note: if there is an extra (empty) line in the data file this routine will put 1 wrong extra term
// DataFile must have 13 columns: energy (in eV), Re[eps_11], Im[eps_11], Re[eps_12], Im[eps_12],
// Re[eps_13], Im[eps_13], Re[eps_22], Im[eps_22], Re[eps_23], Im[eps_23], Re[eps_33], Im[eps_33]
// Epsilon is then assumed to be symmetric
void read_epsilon_data( string DataFile,
						vector<double>& energy,
						vector<vector<dComplex> >& epsilon_Data	)
{

	dComplex I(0,1);
	double dummy_r, dummy_i;
	int i;
	string line;
	ifstream infile(DataFile.c_str());
	if  (!infile.is_open())	
		nrerror("read_epsilon_data: File not found!");

	epsilon_Data.clear();
	epsilon_Data.resize(6);
	for (i=0; i<6; i++)
		epsilon_Data[i].resize(0);
	energy.resize(0);
	
	getline(infile, line);	// this reads the first line which only contains text
	infile >>  dummy_r;		// this reads the lowest energy 
	
	do {
		energy.push_back(dummy_r);
		for (i=0; i<6; i++)
		{
				infile >>  dummy_r >> dummy_i; 
				epsilon_Data[i].push_back(dummy_r+I*dummy_i);
		}
		infile >>  dummy_r; 
	} while ( !infile.eof() );
	infile.close();
	return;
}

// Given the energy, this routine finds the correct epsilon by interpolating data 
void set_epsilon( double hw, dComplexMatrix& epsilon, vector<double>& energy, vector<vector<dComplex> >&  epsilon_Data)
{
	if ( hw<energy[0] || hw>energy.back() )
		nrerror("set_epsilon: not enough data for the requested energy range!");
	mylinint(energy,  epsilon_Data[0], energy.size(), hw, &epsilon(0,0));
	mylinint(energy,  epsilon_Data[1], energy.size(), hw, &epsilon(0,1));
	mylinint(energy,  epsilon_Data[2], energy.size(), hw, &epsilon(0,2));
	mylinint(energy,  epsilon_Data[3], energy.size(), hw, &epsilon(1,1));
	mylinint(energy,  epsilon_Data[4], energy.size(), hw, &epsilon(1,2));
	mylinint(energy,  epsilon_Data[5], energy.size(), hw, &epsilon(2,2));
	epsilon(2,0)=epsilon(0,2);
	epsilon(1,0)=epsilon(0,1);
	epsilon(2,1)=epsilon(1,2);
	return;
}

// A routine to manually set a diagonal epsilon
void set_diag_epsilon( dComplex e11, dComplex e22, dComplex e33, dComplexMatrix&  eps)
{
	eps(0,0)=e11;
	eps(0,1)=0.0;
	eps(0,2)=0.0;
	eps(1,0)=0.0;
	eps(2,0)=0.0;
	eps(1,1)=e22;
	eps(1,2)=0.0;
	eps(2,2)=e33;
	eps(2,1)=0.0;
	return;
}

// A routine to manually set an isotropic epsilon
void set_iso_epsilon( dComplex e11, dComplexMatrix& eps)
{
	eps(0,0)=e11;
	eps(1,1)=e11;
	eps(2,2)=e11;
	eps(0,1)=0.0;
	eps(0,2)=0.0;
	eps(1,0)=0.0;
	eps(2,0)=0.0;
	eps(1,2)=0.0;
	eps(2,1)=0.0;
	return;
}

// Checks if the used epsilons are symmetric
int check_symmetry(dComplexMatrix&  eps)
{
	double emax=0.0;
	for (int i=0; i<3; i++) 
		for (int j=i+1; j<3; j++)
			if (abs(eps(i,j)-eps(j,i))>emax)   emax=abs(eps(i,j)-eps(j,i));
	if (emax < ZERO_TOL ) 
		return(true);
	else
		return(false);
}

// Checks if epsilon corresponds to an isotropic material
int check_isotropy(dComplexMatrix&  eps)
{
	double emax=0.0;
	for (int i=0; i<3; i++) 
		for (int j=0; j<3; j++)
			if (i!=j && abs(eps(i,j))>emax)   emax=abs(eps(i,j));
	if (emax+abs(eps(0,0)-eps(1,1))+abs(eps(0,0)-eps(2,2)) < ZERO_TOL ) 
		return(true);
	else
		return(false);
}

// This routine computes the Tp matrix (see Schubert) 
// It also stores all the properties of the 4 waves propagating in the material
// kz (kz[i] contains the value for wave i=0,..,3), 
// polarization (pol[i] is a vector containing {Ex, Ey, Hx, Hy}, for each wave i=0,..,3)
// Poynting vector direction (the vectors s[i], for each wave i=0,..,3) 
void set_Tp(  
			  dComplexMatrix&  Tp,
			  dComplexMatrix&  eps,
			  double kx,
			  double hw,
			  double d,
			  vector<dComplex>& kz,
			  vector<dComplexVector>& pol,
			  vector<dComplexVector>& s,
			  int info_flag       // if different from 0 data about the waves are printed on a log file
								)
{
	dComplex I(0,1);
	int i, j, m;

	if (Tp.size1() != 4 || Tp.size2() != 4 )
		nrerror("set_Tp: Tp matrix size is wrong!");

	if (pol.size() != 4 || s.size() != 4)
		nrerror("set_Tp: pol or s matrix size is wrong!");
	else
		for (i=0; i<4; i++)
			if (pol[i].size() != 4 || s[i].size() != 3)
				nrerror("set_Tp: pol or s matrix size is wrong!");

	if (kz.size() != 4) 
		nrerror("set_Tp: kz vector size is wrong!");

	dComplexMatrix Delta(4,4);
	double k0=hw/hc;

	dComplex e11=eps(0,0);
	dComplex e12=eps(0,1);
	dComplex e13=eps(0,2);
	dComplex e22=eps(1,1);
	dComplex e23=eps(1,2);
	dComplex e33=eps(2,2);
	dComplex e21=eps(1,0);
	dComplex e31=eps(2,0);
	dComplex e32=eps(2,1);

	Delta(0,0)=-kx*e31/e33;
	Delta(0,1)=-kx*e32/e33;
	Delta(0,2)=0.0;
	Delta(0,3)=1.0-kx*kx/e33;

	Delta(1,0)=0.0;
	Delta(1,1)=0.0;
	Delta(1,2)=-1.0;
	Delta(1,3)=0.0;

	Delta(2,0)=e23*e31/e33-e21;
	Delta(2,1)=kx*kx-e22+e23*e32/e33;
	Delta(2,2)=0.0;
	Delta(2,3)=kx*e23/e33;

	Delta(3,0)=e11-e13*e31/e33;
	Delta(3,1)=e12-e13*e32/e33;
	Delta(3,2)=0.0;
	Delta(3,3)=-kx*e13/e33;

//	cout << "Delta" << Delta << endl; 
	
	dComplexMatrix CMtemp(4,4);
	dComplexMatrix EigVecs(4,4);
	dComplexMatrix* EigVecsPtr=&EigVecs;
	CMtemp=Delta;
	ComplexMatrixEigendecomposition(CMtemp, kz, EigVecsPtr);

//	Eig_Delta(Delta, kz);

//	cout << "kz: " << kz << endl;
//	cout << "Delta: " << Delta;

// Finds the direction of propagation of energy (vector s)

	e11+=1.0*ZERO_TOL;		// This is made to obtain the correct limit for isotropic materials
	e22+=2.0*ZERO_TOL;
	e33+=3.0*ZERO_TOL;

	for (i=0; i<4; i++)
	{
		s[i][0]=(e12*e12+e13*e13-e11*e22-e11*e33)*kx+2.0*e11*kx*kx*kx
		         +kz[i]*(-e13*e22+e12*e23)+3.0*e13*kx*kx*kz[i]
				 +(e11+e33)*kx*kz[i]*kz[i]+e13*kz[i]*kz[i]*kz[i];
	
		s[i][1]=kx*(e13*e23-e12*e33)+e12*kx*kx+e12*kx*kx*kx
				+kz[i]*(e12*e13-e11*e23)+e23*kx*kx*kz[i]+e12*kx*kz[i]*kz[i]+e23*kz[i]*kz[i]*kz[i];

		s[i][2]=kx*(-e13*e22+e12*e23)+kx*kx*kx*e13+kz[i]*(e13*e13+e23*e23-e11*e33-e22*e33)
				+(e11+e33)*kx*kx*kz[i]+3.0*e13*kx*kz[i]*kz[i]+2.0*e33*kz[i]*kz[i]*kz[i];
	}

	dComplex temp;
	for (i=0; i<4; i++)
	{
		temp=kx*s[i][0]+kz[i]*s[i][2];
		for (j=0; j<3; j++)		s[i][j]/=temp;
	}

	e11-=1e-10;   // Here we restore the correct values
	e22-=2e-10;
	e33-=3e-10;

// waves propagating toward z>0 (i.e. with sz>0) are placed in the first two places of kz vector

	if (s[0][2].real() < 0.0) 
	{
		if (s[2][2].real() > 0.0) 
		{
			swap(kz[0],kz[2]);
			swap(s[0], s[2]);
		}
		else
		{			
			swap(kz[0],kz[3]);
			swap(s[0],s[3]);
		}
	}

	if (s[1][2].real() < 0.0) 
	{
		if (s[3][2].real() > 0.0) 
		{
			swap(kz[1],kz[3]);
			swap(s[1],s[3]);
		}
		else
		{
			swap(kz[1],kz[2]);
			swap(s[1],s[2]);
		}
	}
	
	int rank;
	vector<double>	SingularValues(4);
	dComplexMatrix	Eigensystem(4,4);
	dComplexMatrix	range(4,4);
	dComplexMatrix	null_space(4,4);
	dComplexMatrix	debug(4,4);

	for (i=0; i<2; i++)
	{
		for (j=0; j<4; j++) for (m=0; m<4; m++)
			if (j==m) Eigensystem(j,m)=-kz[i]; else Eigensystem(j,m)=0.0;
		Eigensystem+=Delta;
		for (j=0; j<4; j++) for (m=0; m<4; m++) debug(j,m)=Eigensystem(j,m);
		rank=CSVD(Eigensystem, SingularValues, range, null_space);
		if (rank<2 || rank>3)
		{	
			cout << "hw: " << hw << endl;
			cout << "Epsilon: " << eps << endl;
			cout << "Delta: " << Delta << endl;
			cout << "kz: " << kz << endl;
			cout << "eigen1: " << debug;
			cout << "Sing: " << SingularValues << endl;
			cout << "range: " << range << endl;
			cout << "null: " << null_space << endl;
			cout << "rank: " << rank<< endl;
			nrerror("set_Tp: wrong rank!");
		}

		if ( rank == 3 )
		{
			for (j=0; j<4; j++)
				pol[i][j]=null_space(j,3);
		}
		if ( rank == 2 )
		{
			for (j=0; j<4; j++)
			{
				pol[0][j]=null_space(j,2);
				pol[1][j]=null_space(j,3);
			}
			break;
		}
	}
	
	for (i=2; i<4; i++)
	{
			for (j=0; j<4; j++) for (m=0; m<4; m++)
				if (j==m) Eigensystem(j,m)=-kz[i]; else Eigensystem(j,m)=0.0;
			Eigensystem+=Delta;
			rank=CSVD(Eigensystem, SingularValues, range, null_space);
			if (rank<2 || rank>3)
			{	
				cout << "hw: " << hw << endl;
				cout << "Epsilon: " << eps << endl;
				cout << "Delta: " << Delta << endl;
				cout << "kz: " << kz << endl;
				cout << "eigen2: " << Eigensystem;
				cout << "Sing: " << SingularValues << endl;
				cout << "range: " << range << endl;
				cout << "null: " << null_space << endl;
				cout << "rank: " << rank<< endl;
				nrerror("SEt_Tp: wrong rank!");
			}

		if ( rank == 3 )
		{
			for (j=0; j<4; j++)
				pol[i][j]=null_space(j,3);
		}
		if ( rank == 2 )
		{
			for (j=0; j<4; j++)
			{
				pol[2][j]=null_space(j,2);
				pol[3][j]=null_space(j,3);
			}
			break;
		}
		
	} 

	for (i=0; i<4; i++) 
	{
		if (kx < ZERO_TOL)
			temp=pow(pol[i][0],2)+pow(pol[i][1],2);
		else
			temp=pow(pol[i][0],2)*(1.0+pow(kz[i]/kx,2))-2.0*pol[i][0]*pol[i][3]*kz[i]/kx/kx
			      +pow(pol[i][3]/kx,2)+pow(pol[i][1],2);
		for (j=0; j<4; j++)		pol[i][j]/=sqrt(temp);
	}

// We put first the polarization which is mainly along y (i.e. the s polarization in the isotropic case)
			
	if (abs(pol[0][0]) > abs(pol[0][1]) ) 
	{
		swap(pol[0],pol[1]);
		swap(kz[0],kz[1]);
		swap(s[0],s[1]);
	}

// We put first the polarization which is mainly along y (i.e. the s polarization in the isotropic case)

	if (abs(pol[2][0]) > abs(pol[2][1]) ) 
	{
		swap(pol[2],pol[3]);
		swap(kz[2],kz[3]);
		swap(s[2],s[3]);
	}
		
// Following Schubert, waves travelling toward z>0 are placed in the 1st and 3rd position 

	swap(pol[1],pol[2]);
	swap(kz[1],kz[2]);
	swap(s[1],s[2]);
	
// In order to keep the same convention as in Schubert (look at the formulas for LA and Lf not the figure)
// we make sure that s waves have a positive value of Ey,
// forward propagating p waves have a positive value of Ex, 
// backward propagating p waves have a negative value of Ex
	
	if (pol[0][1].real()<0)
		for (j=0; j<4; j++)		pol[0][j]*=-1.0;
	if (pol[1][1].real()<0)
		for (j=0; j<4; j++)		pol[1][j]*=-1.0;
	if (pol[2][0].real()<0)
		for (j=0; j<4; j++)		pol[2][j]*=-1.0;
	if (pol[3][0].real()>0)
		for (j=0; j<4; j++)		pol[3][j]*=-1.0;
	
// This alert usually appears when all the waves are evanescent: just check them
	if (s[0][2].real() < 0.0 ||  s[2][2].real() < 0.0)
	{
		cout <<  "set_Tp alert: waves in 1st and 3rd position are not propagating toward z>0" << endl;
		cout << "(This alert usually appears when all the waves are evanescent: just check them)" << endl;
	}

	double angle_k, angle_s;
	if (info_flag != 0) 
	{
		ofstream logfile("info.log", ios::out | ios::app);

		for (i=0; i<4;i++) 
		{
			logfile  << "Wave " << i ;

			angle_k=atan(kx/kz[i].real())*180.0/Pi;
			angle_s=atan(s[i][0].real()/s[i][2].real())*180.0/Pi;

			if ( fabs(angle_k-angle_s) < 1e-2)
				logfile << "  Ordinary ";
			else
				logfile << "  ExtraOrdinary ";

			if (i%2==0)
				logfile << "  z ----> ";
			else
				logfile << "  <---- z ";

			if ( fabs(kz[i].real()) < fabs(kz[i].imag()) )
				logfile << "Evanescent ";
			logfile << ":" << endl;

			logfile << "k: (" << kx << ", 0, " << kz[i] << "); angle = " 
			    <<  angle_k << endl;
			logfile << "s: (" << s[i][0] << ", " << s[i][1] << ", " << s[i][2] << "); angle = " 
			    <<  angle_s << endl;
			logfile << "ex: " << pol[i][0] << ";  ey: " << pol[i][1]  
				       << ";  hx: " << pol[i][2] << ";  hy: " << pol[i][3]<< endl;
		}
		logfile.close();
	}

// If the material is isotropic, equations for beta cannot be solved

	if (check_isotropy(eps)) 
	{
		for (i=0; i<Tp.size1(); i++) for (j=0; j<Tp.size2(); j++) 	Tp(i,j)=0.0;
		dComplex q=sqrt(e11-kx*kx);
		Tp(0,0)=cos(-k0*d*q);			// minus sign is because we need Tp(-d) 
		Tp(1,1)=cos(-k0*d*q);
		Tp(2,2)=cos(-k0*d*q);
		Tp(3,3)=cos(-k0*d*q);
		Tp(0,3)=I*q/e11*sin(-k0*d*q);
		Tp(1,2)=-I/q*sin(-k0*d*q);
		Tp(2,1)=-I*q*sin(-k0*d*q);
		Tp(3,0)=I*e11/q*sin(-k0*d*q);	
	}
	else
	{
		std::vector<dComplex> vec_4x4(4*4);
		std::vector<dComplex> Beta_Vec(4);
		ublas::identity_matrix<dComplex,	ublas::column_major> Id(4);

		for (i=0; i<4; i++) for (j=0;j<4;j++) vec_4x4[j*4+i]=pow(kz[i], j);

		for (j=0;j<4;j++) Beta_Vec[j]=exp(-I*k0*d*kz[j]);  // minus sign is because we need Tp(-d) 

		gesv_cpp(4, 1, vec_4x4, Beta_Vec);    	         
	
		Tp=Beta_Vec[0]*Id+Beta_Vec[1]*Delta+Beta_Vec[2]*prod(Delta,Delta)
		  +Beta_Vec[3]*prod(Delta,dComplexMatrix(prod(Delta,Delta)));
	}
	return;
}

// Computes matrix T (see Schubert)
void set_T(					   vector<dComplexMatrix >&  epsilon,
							   dComplexMatrix&  T,
							   dComplexMatrix& La,
							   dComplexMatrix& Lf,
							   double kx,
							   double hw,
							   vector<double>& d,
							   int info_flag
								)
{
	if (epsilon.size() != d.size() ) 
		nrerror("set_T: epsilon and d are specified for a different number of media");
	dComplexMatrix Tp(4,4);
	dComplexMatrix Mat_temp(4,4);
	vector<dComplex> kz(4);
	vector<dComplexVector> pol(4);
	vector<dComplexVector> s(4);
	
	if (La.size1() != 4 || La.size2() != 4 )
		nrerror("set_T: La matrix size is wrong!");
	
	if (Lf.size1() != 4 || Lf.size2() != 4 )
		nrerror("set_T: Lf matrix size is wrong!");
	
	if (T.size1() != 4 || T.size2() != 4 )
		nrerror("set_T: T matrix size is wrong!");

	int i, j,n;
	for (i=0; i<4; i++)
	{
		pol[i].resize(4);
		s[i].resize(3);
	}
		
	ofstream logfile("info.log", ios::out | ios::app);
		
	if ( info_flag != 0 )		
	{
		logfile << endl << "-------------------------------------------------------------------" << endl;
		logfile << "Energy : " << hw << " eV" << endl;
		logfile << endl << "Ambient =====================" << endl;
		logfile << "Epsilon: " << epsilon[0] << endl;
	}

	set_Tp(Tp, epsilon[0], kx, hw, d[0], kz, pol, s, info_flag);
	
	for (i=0; i<4; i++)	for (j=0; j<4; j++)
		La(i,j)=pol[j][i];

/*	for (i=0; i<4; i++)	
	{
		cout << "pol " << i << " : ";
		for (j=0; j<4; j++)
			cout << pol[i][j] << " \t";
		cout  << endl;
	}*/
	
	InvertComplexMatrix(La, Mat_temp);
	
//	cout << "Mat_temp: " << Mat_temp << endl;
	for (i=1; i<epsilon.size()-1; i++)
	{	
		if (info_flag != 0)		
		{
				logfile << endl << "Slab n. " << i << " =====================" << endl;
				logfile << "Epsilon: " << epsilon[i] << endl;
		}
		set_Tp(Tp, epsilon[i], kx, hw, d[i], kz, pol, s, info_flag);
//		cout << "Mat_temp:" << Mat_temp << endl;
//		cout << "Tp:" << Tp << endl;
		T=prod(Mat_temp,Tp);
//		cout << "T:" << T << endl;
		Mat_temp=T;
//		cout << "Mat_temp (" << i << "): " << Mat_temp << endl;
	}
	if (info_flag != 0)		
	{
		logfile << endl << "Substrate =====================" << endl;
		logfile << "Epsilon: " << epsilon.back() << endl;
	}
	set_Tp(Tp, epsilon.back(), kx, hw, d.back(), kz, pol, s, info_flag);

/*	Older (working) version
	for (i=0; i<4; i++)	for (j=0; j<4; j++)
		Lf(i,j)=pol[i][j];
	for (i=1; i<4; i+=2) for (j=0; j<4; j++) Lf(i,j)=0.0;
	Lf=trans(Lf);*/
	
	for (i=0; i<4; i++)	
	{
		Lf(i,0)=pol[0][i];
		Lf(i,1)=0.0;
		Lf(i,2)=pol[2][i];
		Lf(i,3)=0.0;
	}
	
/*	for (i=0; i<4; i++)	
	{
		cout << "Lf col " << i << " : ";
		for (j=0; j<4; j++)
			cout << Lf(j,i) << " \t";
		cout  << endl;
	}*/
	
	T=prod(Mat_temp,Lf);
	logfile.close();
	return;
}

// Rotates dMatrix M1 and puts the rotated matrix in M2
// alpha, beta , gamma are Euler angles in degrees
void rotate_matrix(	 dComplexMatrix&  M1, 
					 dComplexMatrix&  M2, 
					 double alpha,
					 double beta,
					 double gamma)
{
	dComplexMatrix RotMat(3,3);
	dComplexMatrix RotMatInv(3,3);
	double a=alpha*Pi/180.0;   // It is easier to have radians in this routine 
	double b=beta*Pi/180.0;
	double c=gamma*Pi/180.0;

	RotMat(0,0)=cos(a)*cos(c)-sin(a)*cos(b)*sin(c);
	RotMat(0,1)=-sin(a)*cos(c)-cos(a)*cos(b)*sin(c);
	RotMat(0,2)=sin(b)*sin(c);
	
	RotMat(1,0)=cos(a)*sin(c)+sin(a)*cos(b)*cos(c);
	RotMat(1,1)=-sin(a)*sin(c)+cos(a)*cos(b)*cos(c);
	RotMat(1,2)=-sin(b)*cos(c); 

	RotMat(2,0)=sin(a)*sin(b);
	RotMat(2,1)=cos(a)*sin(b);
	RotMat(2,2)=1.0*cos(b); 

	InvertComplexMatrix(RotMat,RotMatInv); 
	M2=prod(RotMat, dComplexMatrix(prod(M1, RotMatInv))); 

	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// These are the solutions for Bs, Bp, Cs, Cp (see Schubert) as a function of As, Ap

dComplex Bs(dComplexMatrix& T, 
								  dComplex As, 
								  dComplex Ap)
{
	return(

			( Ap*(-T(0,2)*T(1,0)+T(0,0)*T(1,2))+	As*(-T(1,2)*T(2,0)+T(1,0)*T(2,2)) )/
			(T(0,0)*T(2,2)-T(0,2)*T(2,0))

				);				
}

dComplex Bp(dComplexMatrix& T, 
								  dComplex As, 
								  dComplex Ap)
{
	return(

			( Ap*(-T(0,2)*T(3,0)+T(0,0)*T(3,2)) + As*(-T(3,2)*T(2,0)+T(3,0)*T(2,2)) )/
			(T(0,0)*T(2,2)-T(0,2)*T(2,0))

				);				
}

dComplex Cs(dComplexMatrix& T, 
								  dComplex As, 
								  dComplex Ap)
{
	return(

			( -Ap*T(0,2) +	As*T(2,2) )/
			(T(0,0)*T(2,2)-T(0,2)*T(2,0))

				);				
}

dComplex Cp(dComplexMatrix& T, 
								  dComplex As, 
								  dComplex Ap)
{
	return(

			( Ap*T(0,0) -	As*T(2,0) )/
			(T(0,0)*T(2,2)-T(0,2)*T(2,0))

				);				
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z component of incident, reflected and transmitted 
// time averaged Poynting vectors for each polarization (s or p)

double  SzI( dComplexMatrix& T, 
				dComplexMatrix& La,
				dComplexMatrix& Lf,
				dComplex As, 
				dComplex Ap)
{

	return(
		norm( As ) * real(  La(0,0)*conj(La(3,0)) -  La(1,0)*conj(La(2,0)) )
		+norm( Ap ) * real(  La(0,2)*conj(La(3,2)) -  La(1,2)*conj(La(2,2)) )
		);
}

double  SzRs( dComplexMatrix& T, 
						dComplexMatrix& La,
						dComplexMatrix& Lf,
						dComplex As, 
						dComplex Ap)
{

	return(
		norm( Bs(T, As, Ap) ) * real(  La(0,1)*conj(La(3,1)) -  La(1,1)*conj(La(2,1)) )		
		);

}

double  SzRp( dComplexMatrix& T, 
						dComplexMatrix& La,
						dComplexMatrix& Lf,
						dComplex As, 
						dComplex Ap)
{

	return(
		norm( Bp(T, As, Ap) ) * real(  La(0,3)*conj(La(3,3)) -  La(1,3)*conj(La(2,3)) )		
		);
}

double  SzTs( dComplexMatrix& T, 
						dComplexMatrix& La,
						dComplexMatrix& Lf,
						dComplex As, 
						dComplex Ap)
{

	return(
		norm( Cs(T, As, Ap) ) * real(  Lf(0,0)*conj(Lf(3,0)) -  Lf(1,0)*conj(Lf(2,0)) )	
		);
}

double  SzTp( dComplexMatrix& T, 
						dComplexMatrix& La,
						dComplexMatrix& Lf,
						dComplex As, 
						dComplex Ap)
{

	return(
		norm( Cp(T, As, Ap) ) * real(  Lf(0,2)*conj(Lf(3,2)) -  Lf(1,2)*conj(Lf(2,2)) )
		);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reflection and transmission coefficients for exiting polarizations s and p 
// as a function of As, Ap (which give the incoming polarization)

double  Rs_coeff(	  dComplexMatrix& T, 
										dComplexMatrix& La,
										dComplexMatrix& Lf,
										dComplex As, 
										dComplex Ap)
{

	return(
		
		-SzRs(T, La, Lf, As, Ap) /  SzI(T, La, Lf, As, Ap)

		);
}

double  Rp_coeff(	  dComplexMatrix& T, 
										dComplexMatrix& La,
										dComplexMatrix& Lf,
										dComplex As, 
										dComplex Ap)
{

	return(
		
		-SzRp(T, La, Lf, As, Ap) /  SzI(T, La, Lf, As, Ap)

		);
}

double  Ts_coeff(	  dComplexMatrix& T, 
										dComplexMatrix& La,
										dComplexMatrix& Lf,
										dComplex As, 
										dComplex Ap)
{
	return(
		
		SzTs(T, La, Lf, As, Ap) /  SzI(T, La, Lf, As, Ap)

		);
}

double  Tp_coeff(	  dComplexMatrix& T, 
										dComplexMatrix& La,
										dComplexMatrix& Lf,
										dComplex As, 
										dComplex Ap)
{
	return(
		
	 SzTp(T, La, Lf, As, Ap) /  SzI(T, La, Lf, As, Ap)

		);
}


//////////////////////////////////////////////////////////////////////////////////////
/// Ellipsometry parameters (in degrees) 

dComplex rho(dComplexMatrix& T,
		          dComplex As, 
				  dComplex Ap)
{
	if (abs(As) < ZERO_TOL || abs(Ap) < ZERO_TOL )
		nrerror("rho: incident light cannot be pure s or p !");
	
	return(

		( Bp(T, As, Ap) / Ap ) / (Bs(T, As, Ap) / As )

		);
}

double Phi(dComplexMatrix& T,
		          dComplex As, 
				  dComplex Ap)
{
	return(

		atan( abs( rho(T, As, Ap) ) ) *180.0 / Pi

		);
}

double Delta(dComplexMatrix& T,
		          dComplex As, 
				  dComplex Ap)
{
	return(

		arg( rho(T, As, Ap) ) *180.0 / Pi

		);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Routine that calculates T, La and Lf  
void compute_transfer_matrices(		   int dim,						// number of media (layers+2)
							   vector<dComplexMatrix>  epsilon_eff,	// dielectric tensor for each layer (including ambient and final)
							   vector<double> d,						// thickness of the n slabs (including ambient and final, which are not used)
							   double inc_angle,						// angle of incidence (degree)
							   double gamma,							// angle between the incident plane and xz plane (degree) 
							   double wl,								// wavelength (in microns)
							   vector<dComplexMatrix> Tp,				// Partial transfer matrices Tip (see Schubert's paper) 
							   dComplexMatrix& T,						// General transfer matrix T (see Schubert's paper) 
							   dComplexMatrix& La,						// Incident matrix La (see Schubert's paper)
							   dComplexMatrix& Lf						// Exit matrix Lf (see Schubert's paper)
							   )
{
	int i,j; 
	if (epsilon_eff.size()!=dim || d.size()!=dim)
		nrerror("compute_transfer_matrices: wrong input dim!");
	if (!check_isotropy(epsilon_eff[0])) 
		nrerror("compute_transfer_matrices: ambient is not isotropic!");
	if (!check_isotropy(epsilon_eff[dim-1])) 
		nrerror("compute_transfer_matrices: exit medium is not isotropic!");
	if (epsilon_eff[0](0,0).imag()>ZERO_TOL )
		nrerror("compute_transfer_matrices: ambient medium is absorbing!");
	double na=epsilon_eff[0](0,0).real();
	double kx=na*sin(inc_angle*Pi/180.0);	// in units of k0
	double hw=0.197*2.0*Pi/wl;				// in eV (energy is used in the set_T routine instead of wavelength)
	
	for (i=1; i<dim-1; i++) 
		rotate_matrix(epsilon_eff[i], epsilon_eff[i], 0.0, 0.0, -gamma);
	
	set_T(epsilon_eff, T, La, Lf, kx, hw, d, 0);  // last argument is just a flag for writing info.log
	
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example of a routine that computes absorption for a slab 
// placed between two isotropic media, as a function of energy. 

void spectrum_vs_energy(string DataFile,			// name of the file where epsilon data are
						string OutputFile,			// name of the file where the spectrum is written
						double na,                  // ambient n
						double nf,					// substrate n
						double d_slab,				// thickness of the first slab
						double inc_angle,			// angle of incidence (°)
						double gamma				// angle between the incident plane and xz plane 
						)
{
	dComplexMatrix  epsilon(3, 3);	
	dComplexMatrix  epsilon_rot(3, 3);				// it is initialized to 0
	vector<vector<dComplex> > epsilon_Data;
	vector<double> energy;
	double hw;
	dComplexMatrix  T(4, 4);
	//	dComplexMatrix  T_background(4, 4);
	dComplexMatrix  La(4, 4); 
	dComplexMatrix  Lf(4, 4);
	ofstream file_out;
	double kx=na*sin(inc_angle*Pi/180.0);
	int i, n; 
	int dim=3;	
	
	epsilon.clear();
	epsilon_rot.clear();
	file_out.open(OutputFile.c_str());
	
	std::vector<double> d(dim);						// vector of slab thicknesses in micron  
	d[1]=d_slab;									// (first and last are not used but must be set)											
	std::vector<dComplexMatrix>  epsilon_eff(dim);
	for (i=0; i<dim; i++)
	{
		epsilon_eff[i].resize(4,4);
		epsilon_eff[i].clear();			
	}		// it is initialized to 0
	
	set_iso_epsilon(na, epsilon_eff[0]);								// ambient medium
	set_iso_epsilon(nf, epsilon_eff[dim-1]);								// exit medium
	read_epsilon_data(DataFile, energy, epsilon_Data);			// Reads slab epsilon
	
	file_out <<		"#Energy(eV)\tAbs_s\tAbs_p\tRef_s\tRef_p" //	T_ss  T_sp   T_ps   T_pp   R_ss  R_sp   R_ps   R_pp" 
	<< endl;
	//	cout << "Background epsilon: " << epsilon_Data[0][0].real() << endl;
	for (n=0; n<epsilon_Data[0].size(); n++)
	{
		hw=energy[n]; 
		epsilon(0,0)=epsilon_Data[0][n];
		epsilon(0,1)=epsilon_Data[1][n];
		epsilon(0,2)=epsilon_Data[2][n];
		epsilon(1,1)=epsilon_Data[3][n];
		epsilon(1,2)=epsilon_Data[4][n];
		epsilon(2,2)=epsilon_Data[5][n];
		epsilon(2,0)=epsilon(0,2);
		epsilon(1,0)=epsilon(0,1);
		epsilon(2,1)=epsilon(1,2);
		rotate_matrix(epsilon, epsilon_eff[1], 0.0, 0.0, -gamma);
		
		if (!check_symmetry(epsilon_eff[1]))
			nrerror("spectrum_vs_energy: epsilon is not symmetric!");
		
		set_T(epsilon_eff, T, La, Lf, kx, hw, d, 0);  // last argument is a flag for writing info.log
		
		//		set_iso_epsilon(epsilon_Data[0][0].real(), epsilon_eff[1]);	// background epsilon
		//		set_T(epsilon_eff, T_background, La, Lf, kx, hw, d, 0);		// last argument is a flag for writing info.log
		
		file_out << hw << "\t"
		<< -log10(Ts_coeff(T, La, Lf, 1.0, 0.0)+Tp_coeff(T, La, Lf, 1.0, 0.0))/*+log10(Ts_coeff(T_background, La, Lf, 1.0, 0.0))*/ << "\t" 
		<< -log10(Tp_coeff(T, La, Lf, 0.0, 1.0)+Ts_coeff(T, La, Lf, 0.0, 1.0))/*+log10(Tp_coeff(T_background, La, Lf, 0.0, 1.0))*/ << "\t"
		<< Rs_coeff(T, La, Lf, 1.0, 0.0) + Rp_coeff(T, La, Lf, 1.0, 0.0) << "\t"
		<< Rs_coeff(T, La, Lf, 0.0, 1.0) + Rp_coeff(T, La, Lf, 0.0, 1.0) << "\t"
		//						<< Ts_coeff(T, La, Lf, 1.0, 0.0) << "  " 
		//						<< Tp_coeff(T, La, Lf, 1.0, 0.0) << "  "
		//                      << Tp_coeff(T, La, Lf, 0.0, 1.0) << "  "
		//						<< Ts_coeff(T, La, Lf, 0.0, 1.0) << "  " 
		//						<< Rs_coeff(T, La, Lf, 1.0, 0.0) << "  " 
		//						<< Rp_coeff(T, La, Lf, 1.0, 0.0) << "  "
		//						<< Rs_coeff(T, La, Lf, 0.0, 1.0) << "  " 
		//		                << Rp_coeff(T, La, Lf, 0.0, 1.0) << "  "
		<< endl;
	}
	file_out.close();
	return;
}

template<typename T>
void read_spectrum_input_data(string file_name, T& hw_min, T& hw_max, T& hw_step, T& epsilon_inf, T& gamma_coeff, T& d_slab, 
							  vector<T>& inc_angle, vector<T>& gamma_angle, T& na, T& nf)
{ 
	string line;
	istringstream iss, iss2;
	ifstream infile(file_name.c_str());
	if (!infile.is_open()) nrerror("read_spectrum_input_data: File not found!");
	
	getline(infile, line, ':');
	infile >> d_slab;
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<T>( iss ), istream_iterator<T>(), back_inserter( inc_angle ) );
	iss.clear();
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<T>( iss ), istream_iterator<T>(), back_inserter( gamma_angle ) );
	getline(infile, line, ':');
	infile >> na;
	getline(infile, line, ':');
	infile >> nf;
	//	get_next_from_file(infile, hw_min);
	//	get_next_from_file(infile, hw_max);
	get_next_from_file(infile, hw_step);
	get_next_from_file(infile, epsilon_inf);
	get_next_from_file(infile, gamma_coeff);
	infile.close();
	return;
};

template<typename T>
void read_spectrum_input_data_full(string file_name, T& hw_min, T& hw_max, T& hw_step, T& epsilon_inf, T& gamma_coeff, T& d_slab, 
								   vector<T>& inc_angle, vector<T>& gamma_angle, T& na, T& nf, T& Vcell, vector<string>& dipole_files, string& epsilon_file)
{ 
	string line;
	istringstream iss, iss2;
	ifstream infile(file_name.c_str());
	if (!infile.is_open()) nrerror("read_spectrum_input_data: File not found!");
	
	
	getline(infile, line, ':');
	infile >> d_slab;
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<T>( iss ), istream_iterator<T>(), back_inserter( inc_angle ) );
	iss.clear();
	getline(infile, line, ':');
	getline(infile, line);
	iss.str( line );
	copy( istream_iterator<T>( iss ), istream_iterator<T>(), back_inserter( gamma_angle ) );
	getline(infile, line, ':');
	infile >> na;
	getline(infile, line, ':');
	infile >> nf;
	get_next_from_file(infile, hw_min);
	get_next_from_file(infile, hw_max);
	get_next_from_file(infile, hw_step);
	get_next_from_file(infile, epsilon_inf);
	get_next_from_file(infile, gamma_coeff);
	get_next_from_file(infile, Vcell);
	
	int i;
	int ndipolefiles=-1;
	get_next_from_file(infile, ndipolefiles);
	dipole_files.resize(0);
	
	for (i=0; i<ndipolefiles; i++)
	{
		getline(infile, line, '"');
		getline(infile, line, '"');	
		dipole_files.push_back(line);
	}
	if (ndipolefiles==0)
	{
		getline(infile, line, '"');
		getline(infile, line, '"');	
		epsilon_file=line;
	}
	infile.close();
	return;
};

template <typename T>
void write_absorption_spectrum(T hw_min, T hw_max,  T hw_step, T d_slab,
							   vector<T> inc_angle,
							   vector<T> gamma_angle,
							   T na,
							   T nf,
							   string output_filename, string epsilon_filename)
{
	// Computes spectrum
	
	cout << "Computing absorption spectrum ..." << endl;
	
	ostringstream spectrum_filename;
	for (int i=0; i< inc_angle.size(); i++)
		for (int j=0; j<gamma_angle.size(); j++)
		{
			spectrum_filename.str("");
			spectrum_filename << output_filename
			<< ".absorption.i" << inc_angle[i] << ".g" << gamma_angle[j] << ".txt";
			spectrum_vs_energy(epsilon_filename, spectrum_filename.str(), na, nf, d_slab,
							   inc_angle[i], gamma_angle[j]); 
		}
	return;
}

template<typename T>
void read_dipoles( vector<string> input_files, vector<ublas::vector<complex<T> > > &all_dipoles,
				  vector<T> &all_energies )
{
	complex<T> I(0.0,1.0);
	T hw, re_temp, im_temp;
	int i,dc;
	ublas::vector<complex<T> > cv_temp(3);
	all_dipoles.resize(0);
	all_energies.resize(0);
	
	for (dc=0; dc<input_files.size(); dc++)
	{
		ifstream infile(input_files[dc].c_str());
		if (!infile.is_open()) 
			cout << "write_epsilon: Input file '" << input_files[dc] << "' not found!" << endl;
		else
		{
			cout << "Reading from file: " << input_files[dc].c_str() << endl;
			infile >> hw;
			while (!infile.eof())
			{
				
				all_energies.push_back(hw);
				for (i=0; i<3; i++)
				{
					infile >> re_temp >> im_temp;
					cv_temp[i]=re_temp+I*im_temp;
				}
				all_dipoles.push_back(cv_temp);
				infile >> hw;
			}
			infile.close();
		}
	}
	cout << "N. of dipoles read from file: " << all_dipoles.size() << endl;
	cout << "N. of energies read from file: " << all_energies.size() << endl;
	return;
}

template<typename T>
bool write_epsilon(	vector<string> input_files, string output_file, vector<ublas::vector<complex<T> > > &all_dipoles,
				   vector<T> &all_energies, T hw_step, T epsilon_inf, T cell_volume, T gamma_coeff)
{
	cout << "Writing epsilon ..." << endl;
	complex<T> I(0.0,1.0);
	read_dipoles(input_files, all_dipoles, all_energies);
	if (all_energies.size()==0)
		return(false);
	T hw;
	int i,j,n;
	ublas::matrix<complex<T>,	ublas::column_major> Crystal_epsilon(3,3);
	
	ofstream outfile(output_file.c_str());
	outfile << "En(eV)	Re[e11]	Im[e11]	Re[e12]	Im[e12]	Re[e13]	Im[e13]	Re[e22]	Im[e22]	Re[e23]	Im[e23]	Re[e33]	Im[e33]";
	
	T hw_min=all_energies[0]-0.2;
	T hw_max=all_energies[all_energies.size()-1]+0.2;
	for (hw=hw_min; hw<=hw_max; hw+=hw_step)
	{
		outfile << endl;
		for (i=0; i<3; i++)	for (j=0; j<3; j++)
		{	
			// Writes the backgorund isotropic epsilon
			if (i!=j) 
				Crystal_epsilon(i,j)=0.0;
			else
				Crystal_epsilon(i,i)=epsilon_inf;
			// Adds the excitonic epsilon	
			for (n=0; n<all_energies.size(); n++)
			{
				//  The extra factor 2 needed for spin degeneracy (which is correct for singlet states) must be included 
				//  in the dipole moment value. 
				//				gamma_coeff = 0.03*hw; // for high energy
				Crystal_epsilon(i,j) +=	T(8.0*Pi/1.602)*all_dipoles[n][i]*conj(all_dipoles[n][j])*all_energies[n]/cell_volume
				/(T(pow(all_energies[n],2))-hw*hw-I*gamma_coeff*hw);
			}
		}
		outfile << hw ;
		for (i=0; i<3; i++)	for (j=i; j<3; j++)
			outfile << " " << Crystal_epsilon(i,j).real() << " " << Crystal_epsilon(i,j).imag();
	}
	outfile.close(); 
	return(true);
}

template<typename T>
bool write_epsilon(	vector<string> input_files, string output_file, vector<ublas::vector<complex<T> > > &all_dipoles,
				   vector<T> &all_energies, T hw_min, T hw_max, T hw_step, T epsilon_inf, T cell_volume, T gamma_coeff)
{
	cout << "Writing epsilon ..." << endl;
	complex<T> I(0.0,1.0);
	read_dipoles(input_files, all_dipoles, all_energies);
	if (all_energies.size()==0)
		return(false);
	T hw;
	int i,j,n;
	ublas::matrix<complex<T>,	ublas::column_major> Crystal_epsilon(3,3);
	
	ofstream outfile(output_file.c_str());
	outfile << "En(eV)	Re[e11]	Im[e11]	Re[e12]	Im[e12]	Re[e13]	Im[e13]	Re[e22]	Im[e22]	Re[e23]	Im[e23]	Re[e33]	Im[e33]";
	
	//	T hw_min=all_energies[0]-0.2;
	//	T hw_max=all_energies[all_energies.size()-1]+0.2;
	for (hw=hw_min; hw<=hw_max; hw+=hw_step)
	{
		outfile << endl;
		for (i=0; i<3; i++)	for (j=0; j<3; j++)
		{	
			// Writes the backgorund isotropic epsilon
			if (i!=j) 
				Crystal_epsilon(i,j)=0.0;
			else
				Crystal_epsilon(i,i)=epsilon_inf;
			// Adds the excitonic epsilon	
			for (n=0; n<all_energies.size(); n++)
			{
				//  The extra factor 2 needed for spin degeneracy (which is correct for singlet states) must be included 
				//  in the dipole moment value. 
				//				gamma_coeff = 0.03*hw; // for high energy
				Crystal_epsilon(i,j) +=	T(8.0*Pi/1.602)*all_dipoles[n][i]*conj(all_dipoles[n][j])*all_energies[n]/cell_volume
				/(T(pow(all_energies[n],2))-hw*hw-I*gamma_coeff*hw);
			}
		}
		outfile << hw ;
		for (i=0; i<3; i++)	for (j=i; j<3; j++)
			outfile << " " << Crystal_epsilon(i,j).real() << " " << Crystal_epsilon(i,j).imag();
	}
	outfile.close(); 
	return(true);
}

template void read_spectrum_input_data(string file_name, double& hw_min, double& hw_max, double& hw_step, double& epsilon_inf, double& gamma_coeff, double& d_slab, 
									   vector<double>& inc_angle, vector<double>& gamma_angle, double& na, double& nf);
template void read_spectrum_input_data_full(string file_name, double& hw_min, double& hw_max, double& hw_step, double& epsilon_inf, double& gamma_coeff, double& d_slab, 
											vector<double>& inc_angle, vector<double>& gamma_angle, double& na, double& nf, double& cell_volume, vector<string>& dipole_files, string& epsilon_file);							  


template void read_dipoles( vector<string> input_files, vector<ublas::vector<complex<double> > > &all_dipoles,
						   vector<double> &all_energies );

template bool write_epsilon(	vector<string> input_files, string output_file, vector<ublas::vector<complex<double> > > &all_dipoles,
							vector<double> &all_energies, double hw_step,
							double epsilon_inf, double cell_volume, double gamma_coeff);
template bool write_epsilon(	vector<string> input_files, string output_file, vector<ublas::vector<complex<double> > > &all_dipoles,
							vector<double> &all_energies, double hw_min, double hw_max, double hw_step,
							double epsilon_inf, double cell_volume, double gamma_coeff);				   

template void write_absorption_spectrum(double hw_min, double hw_max,  double hw_step, double d_slab,
										vector<double> inc_angle,
										vector<double> gamma_angle,
										double na,
										double nf,
										string output_filename, string epsilon_filename);

template void read_spectrum_input_data(string file_name, float& hw_min, float& hw_max, float& hw_step,  float& epsilon_inf, float& gamma_coeff,  float& d_slab, 
									   vector<float>& inc_angle, vector<float>& gamma_angle, float& na, float& nf);
template void read_spectrum_input_data_full(string file_name, float& hw_min, float& hw_max, float& hw_step,  float& epsilon_inf, float& gamma_coeff,  float& d_slab, 
											vector<float>& inc_angle, vector<float>& gamma_angle, float& na, float& nf, float& cell_volume, vector<string>& dipole_files, string& epsilon_file);									   

template void read_dipoles( vector<string> input_files, vector<ublas::vector<complex<float> > > &all_dipoles,
						   vector<float> &all_energies );

template bool write_epsilon(	vector<string> input_files, string output_file, vector<ublas::vector<complex<float> > > &all_dipoles,
							vector<float> &all_energies, float hw_step,
							float epsilon_inf, float cell_volume, float gamma_coeff);
template bool write_epsilon(	vector<string> input_files, string output_file, vector<ublas::vector<complex<float> > > &all_dipoles,
							vector<float> &all_energies, float hw_min, float hw_max, float hw_step,
							float epsilon_inf, float cell_volume, float gamma_coeff);							


template void write_absorption_spectrum(float hw_min, float hw_max,  float hw_step, float d_slab,
										vector<float> inc_angle,
										vector<float> gamma_angle,
										float na,
										float nf,
										string output_filename, string epsilon_filename);


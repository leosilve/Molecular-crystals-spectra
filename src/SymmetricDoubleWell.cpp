/* Written by Leonardo Silvestri 2007-2013 */

#include <DEIntegrator.h>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>  
#include <SymmetricDoubleWell.h>
//#include <math.h>
#include <OPutils.h>
#include <N_Functions.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
namespace math = boost::math;

#ifndef ZERO_TOL
#define ZERO_TOL 1e-10
#endif

double find_genFCF(unsigned nug, unsigned nue, double beta, double precision)
{
	FunctionObject f(nug, nue, beta);
	double a;
	for (a=10.0; a>1; a-=0.1)
		if (abs(f(a))>precision*0.01)
			break;
//	cout << a << endl;
//	std::cout << integral << "\n";
//	std::cout << integral << ", " << errorEstimate << ", " << evaluations << "\n";
	return(DEIntegrator<FunctionObject>::Integrate(f, -a, a, precision));
}


// Computes the symmetric double well hamiltonian matrix element between two harmonic oscillator eigenfunctions 
// with quantum numbers m,n. The perturbation to the harmonic oscillator of energy hw0 is 
// a Gaussian with peak height A (in units of energy) and width inversely proportional to alpha (in units of energy)
template <typename FLOAT>
double Hmn(FLOAT hw0, FLOAT A, FLOAT alpha, int m, int n)
{
	if ((m+n)%2==1)
		return(0.0);
	FLOAT sum=0.0;
	for (int l=0; l<=min(m,n); l++)
		if ((m-l)%2==0 && (n-l)%2==0)
			sum+=1.0 / math::factorial<FLOAT>(l) * pow(2.0*hw0/(alpha+hw0),l)*pow(-alpha/(alpha+hw0),(m+n-2*l)/2)
					/math::factorial<FLOAT>((m-l)/2)/math::factorial<FLOAT>((n-l)/2);
	sum*=A*sqrt(math::factorial<FLOAT>(m)*math::factorial<FLOAT>(n)/pow(2.0,m+n)*hw0/(alpha+hw0));
	if (m==n)	sum+=(n+0.5)*hw0;
	return(sum);
}

template <typename FLOAT>
void find_SDW_eigenvectors(FLOAT hw0, FLOAT A, FLOAT alpha, int n, 
						   vector<FLOAT> &eigenvalues, vector<FLOAT>  &eigenvectors)
{
	
	cout << "Finding DW eigenstates ... " << endl;

	int i,j;
	
	if (eigenvectors.size()!=n*n)
		nrerror("Find_SW_eigenvectors: error in the dimension of eigvec");
		
	for (i=0; i<n; i++)
		for (j=i; j<n; j++)
			eigenvectors[j*n+i]=Hmn(hw0, A, alpha, i, j);

	syev_cpp('V', 'U', n, eigenvectors, eigenvalues, 'O');

	cout << "Check eigenvalues from syev_cpp!!!!" << endl;
// Eigenvectors are already normalized
/*	double norm_L2;
	for (i=0; i<n; i++)
	{
		norm_L2=0.0;
		for (j=0; j<n; j++)
			norm_L2+=pow(fabs(eigenvectors[i][j]),2);
		for (j=0; j<n; j++)
			eigenvectors[i][j]/=sqrt(norm_L2);
	}*/
	return;
}

void initialize_SDW_parameters(string DW_datafile, double& DW_hw0, double& DW_alpha, double& DW_A, double& DW_Vd,
							   double& DW_hwe, double& DW_lambdae, unsigned int& DW_basis_dim, unsigned int& DW_nug_max,
							   unsigned int& DW_nue_max, vector<vector<double> >& DW_FCF_table,
								vector<double>&	DW_levels)
	{
		cout << "Reading DW data ..." << endl;
		int i,j,n,m;
		double temp;
		string /*genFCF_file,*/ line;
		ifstream infile(DW_datafile.c_str());
		if (!infile.is_open())
			nrerror("Phonon Cloud constructor: DW_datafile not found!");
		getline(infile, line, ':');
		infile >> DW_hw0;
		getline(infile, line, ':');
		infile >> DW_alpha;
		getline(infile, line, ':');
		infile >> DW_A;
		getline(infile, line, ':');
		infile >> DW_Vd;
		getline(infile, line, ':');
		infile >> DW_hwe;
		getline(infile, line, ':');
		infile >> DW_lambdae;
		getline(infile, line, ':');
		infile >> DW_basis_dim;
		getline(infile, line, ':');
		infile >> DW_nug_max;
		getline(infile, line, ':');
		infile >> DW_nue_max;
		infile.close();

// Finds the generalized Franck Condon Factors needed for SDW Hamiltonian
		cout << "Finding generalized Franck Condon Factors ..." << endl;
		vector<vector<double> > FCFsw(DW_nug_max);
		double precision=1e-5;
		for (i=0; i<DW_nug_max; i++)	
		{
			FCFsw[i].resize(DW_nue_max);
			for (j=0; j<DW_nue_max; j++)
				FCFsw[i][j]=find_genFCF(i,j,DW_hw0/DW_hwe,precision);
		}
		ofstream outfile("genFCF.txt"); 
		for (n=0; n<DW_nug_max; n++)	for (m=0; m<DW_nue_max; m++)
			outfile <<  n << "\t" <<  m << "\t" << FCFsw[n][m] << endl ;
		outfile.close();

// Finds SDW eigenstates
		vector<double> eigvec(DW_basis_dim*DW_basis_dim);
		double sum;
		find_SDW_eigenvectors(DW_hw0, DW_A, DW_alpha, DW_basis_dim,  DW_levels, eigvec);
		ofstream outfile2("Double Well Energy Levels.txt"); 
		for (n=0; n<DW_levels.size(); n++)
		{
			DW_levels[n]-=DW_Vd;
			outfile2 << n<< "\t" << DW_levels[n] << endl;
		}
		outfile2.close();

// Finds final Franck Condon Factors between the SDW and a parabola
		cout << "Finding final Franck Condon Factors ..." << endl;
		DW_FCF_table.resize(DW_nug_max);
		for (n=0; n<DW_nug_max; n++)	
		{
			DW_FCF_table[n].resize(DW_nue_max);
			for (m=0; m<DW_nue_max; m++)
			{
				sum=0.0;
				for (j=0; j<DW_basis_dim; j++)
//					sum+=eigvec[n][j]*FCF(j, hw0, m, hwe, lambdae);
					sum+=eigvec[n*DW_basis_dim+j]*FCFsw[j][m];
				if (fabs(sum)>ZERO_TOL) 
					DW_FCF_table[n][m]=sum;
				else
					DW_FCF_table[n][m]=0.0;
			}
		}
		outfile.open("FCF.txt");
		for (n=0; n<DW_nug_max; n++)	for (m=0; m<DW_nue_max; m++)
			outfile <<  n << "\t" <<  m << "\t" << DW_FCF_table[n][m] << endl ;
		outfile.close();
		
		cout << "Check FCF^2 sum for double well states" << endl;
		for (n=0; n<DW_nug_max; n++)	
		{
			sum=0.0;
			for (m=0; m<DW_nue_max; m++)	sum+=pow(DW_FCF_table[n][m],2);
			cout << n << "\t" << sum << endl;
		}
		cout << "Check FCF^2 sum for excited harmonic well states" << endl;
		for (m=0; m<DW_nue_max; m++)	
		{
			sum=0.0;
			for (n=0; n<DW_nug_max; n++)	 sum+=pow(DW_FCF_table[n][m],2);
			cout << m << "\t" << sum << endl;
		}
	return;
}


template double Hmn(float hw0, float A, float alpha, int m, int n);
template void find_SDW_eigenvectors(float hw0, float A, float alpha, int n, vector<float> &eigenvalues, vector<float>  &eigenvectors); 

template double Hmn(double hw0, double A, double alpha, int m, int n);
template void find_SDW_eigenvectors(double hw0, double A, double alpha, int n, vector<double> &eigenvalues, vector<double>  &eigenvectors); 
/* Written by Leonardo Silvestri 2007-2013 */

#ifndef _ME_H
#define _ME_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <string>
#include <vector>
//#include <OPutils.h>
#include <lattice.h>
 
namespace ublas = boost::numeric::ublas;

///////////////////////////////////////////////////////////////////////////
/////////          Molecular Exciton with phonon cloud class
///////////////////////////////////////////////////////////////////////////

template <class FLOAT>
class ME {
	
	typedef ublas::vector<FLOAT>	Vector;
	
private:
	int		nmodes;				// number of vibrational modes
	int		max_cloud;			// Max number of molecules on which the cloud is spreaded
	Vector	k;					// Exciton wave vector
	int		mol;				// inequivalent molecule on which electronic excitation resides 
	std::vector<std::vector<int> > cloud;		// phonon clouds for each mode
	int		np;					// total number of particles on which excitations effectively reside
	
public:
	// Short exciton constructor
	ME() {
		k.resize(3);
		mol=0; 
		np=0;
		nmodes=0;
		max_cloud=0;
		cloud.resize(0);
	};
	// Long exciton constructor with multiple modes
	ME(FLOAT _kx, FLOAT _ky, FLOAT _kz, int _mol, std::vector<std::vector<int> >& _cloud);
	
	// Computes the number of particles for the exciton and modifies np accordingly
	void set_np();
	void set_k_3d(Vector param); 
	void set_k_3d(FLOAT _kx, FLOAT _ky, FLOAT _kz);
	void set_mol(int _mol) {mol=_mol;};
	void set_cloud(int _mode, int _npos, int _nvib);
	void set_nmodes(int _nmodes);
	void set_max_cloud(int _max_cloud);
	
	Vector	get_k()			const {return(k);};
	int		get_mol()		const {return(mol);};
	int		get_np()		const {return(np);};
	int		get_nmodes()	const {return(nmodes);};
	int		get_max_cloud() const {return(max_cloud);};
	int		get_cloud_size()		const { return(cloud.size()); };
	int		get_cloud_size(int i)	const { return(cloud[i].size()); };
	std::string get_string()				const ;
	int			get_cloud(int mode, int i)	const;
	
	bool operator==(const ME<FLOAT>& exc) const;
	
	// Returns the number of vibrational quanta at relative position r
	int get_vib(int mode, lattice<FLOAT>*	latticePtr, Vector r);

	
};

// Screen output of a molecular exciton 
template <typename FLOAT>
std::ostream &operator<<(std::ostream& stream, const ME<FLOAT>& exc);   

#endif

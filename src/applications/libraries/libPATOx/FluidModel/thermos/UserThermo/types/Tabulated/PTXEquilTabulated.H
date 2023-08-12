/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::PTXEquilTabulated

Description
    Tabulated thermo in function of pressure, temperature and elemental mole fractions

SourceFiles
    PTXEquilTabulated.H
    PTXEquilTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef PTXEquilTabulated_H
#define PTXEquilTabulated_H

#include "basicTabulated.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class PTXEquilTabulated
:
public basicTabulated
{

 protected:
  //- Current Argon mass fractions [-]
  scalar currentYAr_;

  // table pointers
  scalarList * pressureTable;
  scalarList * temperatureTable;
  List<scalarList> * elementsMassFractionTable;
  scalarList * argonMassFractionTable;
  List<List<scalarList>> * eTable;
  List<List<scalarList>> * psiTable;
  List<List<scalarList>> * alphaTable;
  List<List<scalarList>> * cpTable;
  List<List<scalarList>> * cvTable;
  List<List<scalarList>> * sigmaTable;
  List<List<scalarList>> * muTable;
  List<List<List<scalarList>>> * equilDiffFluxPTable;
  List<List<List<scalarList>>> * equilDiffFluxTTable;
  List<List<List<scalarList>>> * equilDiffFluxXTable;


  //- Argon index in the mixture
  int indexAr_;

  //- Update all thermodynamic and transport properties for current conditions in function of pressure, temperature and elemental mass fractions
  void updateCurrentPTY(scalar& p, scalar& T, scalarList& Y) ;

  //- Update all thermodynamic and transport properties for current conditions in function of density, sensible energy and elemental mass fractions
  void updateCurrentRhoEY(scalar& rho, scalar& Es, scalarList& Y);

  //- Update all thermodynamic and transport properties on the boundaries
  virtual void updateBoundaries();

  //- Trilinear interpolation
  // @param x     		- x table
  // @param y     		- y table
  // @param z     		- z table
  // @param f_xyz 		- f(x,y,z) table
  // @param x0    		- x0 value
  // @param y0    		- y0 value
  // @param z0    		- z0 value
  // @return      		- f(x0,y0,z0) linearly interpolated value
  scalar trilinearInterpolation(scalarList x, scalarList y, scalarList z, List<List<scalarList> > f_xyz, scalar x0, scalar y0, scalar z0);

  //- Generate and write the thermo table
  //  @param fileName	- thermo table file name
  virtual void generateTable(const word tableName);

  //- Read the thermo table
  //  @param fileName	- thermo table file name
  virtual void readTable(const word tableName);

  //- Find the temperature for a given sensible energy
  //  @param rho			- Density [J/kg]
  //  @param Es				- Sensible energy [J/kg]
  //  @param p				- Initial pressure [Pa]
  //									- on return, computed pressure
  //  @param T				- Initial temperature [K]
  //									- on return, computed temperature
  //  @param Y				- Elemental mass fractions [-]
  //	@param tol			- tolerance for iterative Newton method
  //	@param maxIter	- maximum iterations for iterative Newton method
  void getPTYfromRhoEY(scalar rho, scalar Es, scalar p, scalar T, scalarList Y, scalar tol = 1e-4, scalar maxIter = 1e2);

  //- Compute the diffusive fluxes for element mass and energy
  //	@param fluxP		- on return, the List of fluxes with respect
  //										to pressure (number of elements + 1)
  //	@param fluxT		- on return, the List of fluxes with respect
  //										to temperature (number of elements + 1)
  //	@param fluxX		- on return, the List of fluxes with respect
  //										to element mole fractions
  //										((number of elements + 1) * number of elements)
  void getElementFluxes(scalarList& fluxP, scalarList& fluxT, scalarList& fluxX);

  //- Compute the species diffusion velocities
  //- @param fluxP		- local elements mass fluxes per pressure gradient (list of scalars)
  //- @param fluxT		- local elements mass fluxes per temperature gradient (list of scalars)
  //- @param fluxX		- local elements mass fluxes per elements mole fractions gradient (list of scalars)
  //-	@param gradP		- local gradient of pressure (vector)
  //-	@param gradT		- local gradient of temperature (vector)
  //-	@param gradXe		- local gradient of elements mole fractions (list of vectors)
  //-	@param diffMF		- on return, the List of element difusion velocities
  void getDiffMassFluxes(const scalarList& fluxP, const scalarList& fluxT, const scalarList& fluxX, const vector& gradP, const vector& gradT, const List<vector>& gradXe, List<vector>& diffMF);

  //- Compute the diffusive heat flux in finite-rate chemistry
  //- @param eFluxP		- local species mass fluxes per pressure gradient (list of scalars)
  //- @param eFluxT		- local species mass fluxes per temperature gradient (list of scalars)
  //- @param eFluxX		- local species mass fluxes per elements mole fractions gradient (list of scalars)
  //-	@param gradP		- local gradient of pressure (vector)
  //-	@param gradT		- local gradient of temperature (vector)
  //-	@param gradXe		- local gradient of elements mole fractions (list of vectors)
  //-	@param qDiff		- on return, the vector of diffusive heat flux (vector)
  void getDiffHeatFlux(const scalar& eFluxP, const scalar& eFluxT, const scalarList& eFluxX, vector& gradP, vector& gradT, List<vector>& gradXe, vector& qDiff);





 public:

  //- Runtime type information
  TypeName("PTXEquilTabulated");


  // Constructors

  //- Construct from mesh and phase name
  PTXEquilTabulated(const fvMesh& mesh, const word& dictName);

  //- Destructor
  virtual ~PTXEquilTabulated();


  // Member Functions

  //- Update thermodynamic and transport properties for the whole field
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "PTXEquilTabulated.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

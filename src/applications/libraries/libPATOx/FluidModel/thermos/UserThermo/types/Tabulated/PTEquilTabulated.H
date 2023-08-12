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
    Foam::PTEquilTabulated

Description
    Tabulated thermo in function of pressure and temperature

SourceFiles
    PTEquilTabulated.H
    PTEquilTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef PTEquilTabulated_H
#define PTEquilTabulated_H

#include "basicTabulated.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class PTEquilTabulated
:
public basicTabulated
{

 protected:

  // table pointers
  scalarList * pressureTable;
  scalarList * temperatureTable;
  List<scalarList> * eTable;
  List<scalarList> * psiTable;
  List<scalarList> * alphaTable;
  List<scalarList> * cpTable;
  List<scalarList> * cvTable;
  List<scalarList> * sigmaTable;
  List<scalarList> * muTable;

  //- Update all thermodynamic and transport properties for current conditions in function of pressure and temperature
  void updateCurrentPT(scalar& p, scalar& T) ;

  //- Update all thermodynamic and transport properties for current conditions in function of density and sensible energy
  void updateCurrentRhoE(scalar& rho, scalar& Es);

  //- Update all thermodynamic and transport properties on the boundaries
  virtual void updateBoundaries();

  //- Bilinear interpolation
  // @param x   			- x table
  // @param y   			- y table
  // @param f_xy			- f(x,y) table
  // @param x0				- x0 value
  // @param y0				- y0 value
  // @return					- f(x0,y0) linearly interpolated value
  scalar bilinearInterpolation(scalarList x, scalarList y, List<scalarList> f_xy, scalar x0, scalar y0);

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
  //	@param tol			- tolerance for iterative Newton method
  //	@param maxIter	- maximum iterations for iterative Newton method
  void getPTfromRhoE(scalar& rho, scalar& Es, scalar& p, scalar& T, scalar tol = 1e-4, scalar maxIter = 1e2);


 public:

  //- Runtime type information
  TypeName("PTEquilTabulated");


  // Constructors

  //- Construct from mesh and phase name
  PTEquilTabulated(const fvMesh& mesh, const word& dictName);

  //- Destructor
  virtual ~PTEquilTabulated();


  //- Member Functions

  //- Update thermodynamic and transport properties for the whole field
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "PTEquilTabulated.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

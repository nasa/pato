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
    Foam::basicTabulated

Description
    Mother class for tabulated thermo

SourceFiles
    basicTabulated.H
    basicTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef basicTabulated_H
#define basicTabulated_H

#include "basicUserThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class basicTabulated Declaration
\*---------------------------------------------------------------------------*/

class basicTabulated
:
public basicUserThermo
{

 protected:
  //- Current specif heat at constant volume [J/kg/K]
  scalar currentCv_;

  //- Current specif heat at constant pressure [J/kg/K]
  scalar currentCp_;

  //- Thermo IO dictionary
  IOdictionary thermoDict_;

  //- Tabulated data file name
  const word tableName_;

  //- Tabulated data file path
  const fileName tableFileName_;


  //- Update all thermodynamic and transport properties on the boundaries
  //   ---   PURE VIRTUAL
  virtual void updateBoundaries() = 0;

  //- Compute the element composition (mass fractions, mole fractions, densities)
  //	@param Xe			- on return, the List of element mole fractions
  //	@param Ye			- on return, the List of element mass fractions
  //	@param rhoe		- on return, the List of element densities
  void getElementsComposition(scalarList& Xe, scalarList& Ye, scalarList& rhoe);

  //- Compute the species composition (mass fractions, mole fractions, densities)
  //	@param Xs			- on return, the List of species mole fractions
  //	@param Ys			- on return, the List of species mass fractions
  //	@param rhos		- on return, the List of species densities
  void getSpeciesComposition(scalarList& Xs, scalarList& Ys, scalarList& rhos);

  //- Compute the mixture gas constant
  //	@param mixR		- on return, the mixture gas constant
  void getMixtureR(scalar& mixR);

  //- Compute the temperature
  //	@param mixMM	- on return, the mixture molar mass
  void getMixtureMolarMass(scalar& mixMM);

//	//- Compute the temperature
//	//	@param T			- on return, the temperature
//	void getT(scalar& T);

//	//- Compute the mixture energy
//	//	@param e			- on return, the mixture energy
//	void getE(scalar& e);

//	//- Compute the mixture ethalpy
//	//	@param h			- on return, the mixture enthalpy
//	void getH(scalar& h);

//	//- Compute the pressure
//	//	@param p			- on return, the pressure
//	void getP(scalar& p);

//	//- Compute the density
//	//	@param rho		- on return, the density
//	void getRho(scalar& rho);

//	//- Compute the compressibility
//	//	@param psi		- on return, the compressibility
//	void getPsi(scalar& psi);

//	//- Compute the dynamic viscosity
//	//	@param mu			- on return, the dynamic viscosity
//	void getMu(scalar& mu);

//	//- Compute the kinematic viscosity
//	//	@param nu			- on return, the kinematic viscosity
//	void getNu(scalar& nu);

//	//- Compute the electrical conductivity
//	//	@param sigma	- on return, the electrical conductivity
//    virtual void getSigma(scalar& sigma) = 0;

//	//- Compute the specific heats ratio (isentropic exponent)
//	//	@param gamma	- on return, the isentropic exponent
//	void getGamma(scalar& gamma);

//	//- Compute the thermal conductivity
//	//	@param kappa	- on return, the thermal conductivity
//	void getKappa(scalar& kappa);

//    //- Compute the specific heat at constant volume [J/kg/K]
//    //	@param cv	- on return, the specific heat at constant volume
//    virtual void getCv(scalar& cv);

//	//- Compute the thermal diffusivity
//	//	@param alpha	- on return, the thermal diffusivity
//    virtual void getAlpha(scalar& alpha) = 0;

  //- Generate and write the thermo table
  //  @param fileName  - thermo table file name
  virtual void generateTable(const word tableName) = 0;

  //- Read the thermo table
  //  @param fileName  - thermo table file name
  virtual void readTable(const word tableName) = 0;

  //- Change the environment variable in a file name
  //  @param fileNameOrigin  - origin file name
  //  @return                - return the full path file name
  fileName changeEnviVar(const fileName fileNameOrigin_);

 public:

  //- Runtime type information
  TypeName("basicTabulated");


  // Constructors

  //- Construct from mesh and phase name
  basicTabulated(const fvMesh& mesh, const word& dictName);

  //- Destructor
  virtual ~basicTabulated();


  //- Member Functions
  virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "basicTabulated.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

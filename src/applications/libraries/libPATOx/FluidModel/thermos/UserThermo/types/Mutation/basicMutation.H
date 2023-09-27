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
    Foam::basicMutation

Description
    Mother class for Mutation++ thermo

SourceFiles
    basicMutation.H
    basicMutation.C

\*---------------------------------------------------------------------------*/

#ifndef basicMutation_H
#define basicMutation_H

#include "basicUserThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class basicMutation Declaration
\*---------------------------------------------------------------------------*/

class basicMutation
:
public basicUserThermo
{

 protected:

  //- Update all thermodynamic and transport properties for current conditions
  //   ---   PURE VIRTUAL
  virtual void updateCurrent() = 0;

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

  //- Compute the temperature
  //	@param T			- on return, the temperature
  void getT(scalar& T);

  //- Compute the mixture energy
  //	@param e			- on return, the mixture energy
  void getE(scalar& e);

  //- Compute the mixture ethalpy
  //	@param h			- on return, the mixture enthalpy
  void getH(scalar& h);

  //- Compute the pressure
  //	@param p			- on return, the pressure
  void getP(scalar& p);

  //- Compute the density
  //	@param rho		- on return, the density
  void getRho(scalar& rho);

  //- Compute the compressibility
  //	@param psi		- on return, the compressibility
  void getPsi(scalar& psi);

  //- Compute the dynamic viscosity
  //	@param mu			- on return, the dynamic viscosity
  void getMu(scalar& mu);

  //- Compute the kinematic viscosity
  //	@param nu			- on return, the kinematic viscosity
  void getNu(scalar& nu);

  //- Compute the electrical conductivity
  //	@param sigma	- on return, the electrical conductivity
  void getSigma(scalar& sigma);

  //- Compute the species thermal diffusion ratios
  //	@param tDiffR	- on return, the species net mass production rates (list of scalars)
  void getThermalDiffRatios(scalarList& tDiffR);

 public:

  //- Runtime type information
  TypeName("basicMutation");


  // Constructors

  //- Construct from mesh and phase name
  basicMutation(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~basicMutation();


  //- Member Functions

  //- Update thermodynamic and transport properties in the field
  virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "basicMutation.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

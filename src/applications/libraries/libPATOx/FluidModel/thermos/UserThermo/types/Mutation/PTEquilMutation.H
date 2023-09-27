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
    Foam::PTEquilMutation

Description
    Mutation thermo in function of pressure and temperature

SourceFiles
    PTEquilMutation.H
    PTEquilMutation.C

\*---------------------------------------------------------------------------*/

#ifndef PTEquilMutation_H
#define PTEquilMutation_H

#include "basicMutation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class PTEquilMutation
:
public basicMutation
{

 protected:

  //- Update all thermodynamic and transport properties for current conditions
  virtual void updateCurrent();

  //- Update all thermodynamic and transport properties on the boundaries
  virtual void updateBoundaries();

  //- Set the state for the current conditions
  //	@param var1 	- the pressure or the density
  //	@param var2 	- the temperature or the mixture specific energy (mass)
  //	@param vars 	- integer that specifies the set of inputs
  //									vars = 0 :		var1 = rho	--		var2 = e
  //									vars = 1 :		var1 = p		--		var2 = T
  virtual void setState(const scalar& var1, const scalar& var2, const int vars);

  //- Compute the specific heats ratio (isentropic exponent)
  //	@param gamma	- on return, the isentropic exponent
  void getGamma(scalar& gamma);

  //- Compute the thermal conductivity
  //	@param kappa	- on return, the thermal conductivity
  void getKappa(scalar& kappa);

  //- Compute the thermal diffusivity
  //	@param alpha	- on return, the thermal diffusivity
  void getAlpha(scalar& alpha);

 public:

  //- Runtime type information
  TypeName("PTEquilMutation");


  // Constructors

  //- Construct from mesh and phase name
  PTEquilMutation(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PTEquilMutation();


  //- Member Functions

  //- Update thermodynamic and transport properties for the whole field
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "PTEquilMutation.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

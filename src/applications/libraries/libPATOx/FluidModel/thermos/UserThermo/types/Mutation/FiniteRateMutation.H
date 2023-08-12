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
    Foam::FiniteRateMutation

Description
    Mutation thermo in function of pressure, temperature and species mole fractions

SourceFiles
    FiniteRateMutation.H
    FiniteRateMutation.C

\*---------------------------------------------------------------------------*/

#ifndef FiniteRateMutation_H
#define FiniteRateMutation_H

#include "basicMutation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class FiniteRateMutation
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
  //	@param Ys   	- the List of element mass fractions
  //	@param vars 	- integer that specifies the set of inputs
  //									vars = 0 :		var1 = rho	--		var2 = e
  //									vars = 1 :		var1 = p		--		var2 = T
  virtual void setState(const scalar& var1, const scalar& var2, const scalarList& Ys, const int vars);

  //- Compute the specific heats ratio (isentropic exponent)
  //	@param gamma	- on return, the isentropic exponent
  void getGamma(scalar& gamma);

  //- Compute the thermal conductivity
  //	@param kappa	- on return, the thermal conductivity
  void getKappa(scalar& kappa);

  //- Compute the thermal diffusivity
  //	@param alpha	- on return, the thermal diffusivity
  void getAlpha(scalar& alpha);

  //- Compute the species net mass production rates
  //	@param omegas	- on return, the species net mass production rates (list of scalars)
  void getOmegasMass(scalarList& omegas);

  //- Compute the species diffusion driving forces
  //	@param gradP	- local gradient of pressure (vector)
  //	@param gradT	- local gradient of temperature (vector)
  //	@param gradXs	- local gradient of species mole fractions (list of vectors)
  //	@param dFs		- on return, the List of species diffusion driving forces
  void getDiffDrivingForces(vector& gradP, vector& gradT, List<vector>& gradXs, List<vector>& dFs);

  //- Compute the species diffusion velocities
  //	@param gradP	- local gradient of pressure (vector)
  //	@param gradT	- local gradient of temperature (vector)
  //	@param gradXs	- local gradient of species mole fractions (list of vectors)
  //	@param diffVs	- on return, the List of species diffusion velocities
  void getDiffVelocities(vector& gradP, vector& gradT, List<vector>& gradXs, List<vector>& diffVs);

  //- Compute the diffusive heat flux in finite-rate chemistry
  //	@param rhos		- local species densities (list of scalars)
  //	@param diffVs - local diffusion velocities (list of vectors)
  //	@param qDiff	- on return, the vector of diffusive heat flux (vector)
  void getDiffHeatFlux(const scalarList& rhos, const List<vector>& diffVs, vector& qDiff);

 public:

  //- Runtime type information
  TypeName("FiniteRateMutation");


  // Constructors

  //- Construct from mesh and phase name
  FiniteRateMutation(const fvMesh& mesh, const word& dictName);

  //- Destructor
  virtual ~FiniteRateMutation();


  //- Member Functions

  //- Update thermodynamic and transport properties for the whole field
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "FiniteRateMutation.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

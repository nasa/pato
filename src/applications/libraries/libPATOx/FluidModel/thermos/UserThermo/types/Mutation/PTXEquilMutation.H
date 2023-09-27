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
    Foam::PTXEquilMutation

Description
    Mutation thermo in function of pressure, temperature and elemental mole fractions

SourceFiles
    PTXEquilMutation.H
    PTXEquilMutation.C

\*---------------------------------------------------------------------------*/

#ifndef PTXEquilMutation_H
#define PTXEquilMutation_H

#include "basicMutation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class PTXEquilMutation
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
  //	@param Ye   	- the List of element mass fractions
  //	@param vars 	- integer that specifies the set of inputs
  //									vars = 0 :		var1 = rho	--		var2 = e
  //									vars = 1 :		var1 = p		--		var2 = T
  virtual void setState(const scalar& var1, const scalar& var2, const scalarList& Ye, const int vars);

  //- Compute the specific heats ratio (isentropic exponent)
  //	@param gamma	- on return, the isentropic exponent
  void getGamma(scalar& gamma);

  //- Compute the thermal conductivity
  //	@param kappa	- on return, the thermal conductivity
  void getKappa(scalar& kappa);

  //- Compute the thermal diffusivity
  //	@param alpha	- on return, the thermal diffusivity
  void getAlpha(scalar& alpha);

  //- Compute the diffusive fluxes for element mass and energy
  //	@param fluxP	- on return, the List of fluxes with respect
  //									to pressure (number of elements + 1)
  //	@param fluxT	- on return, the List of fluxes with respect
  //									to temperature (number of elements + 1)
  //	@param fluxX	- on return, the List of fluxes with respect
  //									to element mole fractions
  //									((number of elements + 1) * number of elements)
  void getElementFluxes(scalarList& fluxP, scalarList& fluxT, scalarList& fluxX);

  //- Compute the species diffusion velocities
  //- @param fluxP	- local elements mass fluxes per pressure gradient (list of scalars)
  //- @param fluxT	- local elements mass fluxes per temperature gradient (list of scalars)
  //- @param fluxX	- local elements mass fluxes per elements mole fractions gradient (list of scalars)
  //-	@param gradP	- local gradient of pressure (vector)
  //-	@param gradT	- local gradient of temperature (vector)
  //-	@param gradXe	- local gradient of elements mole fractions (list of vectors)
  //-	@param diffMF	- on return, the List of element diffusion mass fluxes
  void getDiffMassFluxes(const scalarList& fluxP, const scalarList& fluxT, const scalarList& fluxX, const vector& gradP, const vector& gradT, const List<vector>& gradXe, List<vector>& diffMF);

  //- Compute the diffusive heat flux in finite-rate chemistry
  //- @param eFluxP	- local species mass fluxes per pressure gradient (list of scalars)
  //- @param eFluxT	- local species mass fluxes per temperature gradient (list of scalars)
  //- @param eFluxX	- local species mass fluxes per elements mole fractions gradient (list of scalars)
  //-	@param gradP	- local gradient of pressure (vector)
  //-	@param gradT	- local gradient of temperature (vector)
  //-	@param gradXe	- local gradient of elements mole fractions (list of vectors)
  //-	@param qDiff	- on return, the vector of diffusive heat flux (vector)
  void getDiffHeatFlux(const scalar& eFluxP, const scalar& eFluxT, const scalarList& eFluxX, vector& gradP, vector& gradT, List<vector>& gradXe, vector& qDiff);

 public:

  //- Runtime type information
  TypeName("PTXEquilMutation");


  // Constructors

  //- Construct from mesh and phase name
  PTXEquilMutation(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PTXEquilMutation();


  //- Member Functions

  //- Update thermodynamic and transport properties for the whole field
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "PTXEquilMutation.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

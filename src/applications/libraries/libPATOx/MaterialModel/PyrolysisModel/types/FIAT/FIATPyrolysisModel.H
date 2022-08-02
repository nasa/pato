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
    along with OpenFOAM.  If FIATt, see <http://www.gnu.org/licenses/>.

Class
    Foam::FIATPyrolysisModel

Description

SourceFiles
    FIATPyrolysisModel.H
    FIATPyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef FIATPyrolysisModel_H
#define FIATPyrolysisModel_H

#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FIATPyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class FIATPyrolysisModel
:
public simplePyrolysisModel
{

 protected:

  // Mesh reference
  const fvMesh& mesh_;

  // resin volume fraction [-]
  scalar GAMMA;

  // porosity of virgin material [-]
  scalar PHI;

  // initial density of component i [kg/m3]
  List<scalar> RVI;
  PtrList<dimensionedScalar> RVI_;

  // residual density of component i [kg/m3]
  List<scalar> RCI;
  PtrList<dimensionedScalar> RCI_;

  // pre-exponential factor for component i [1/s]
  List<scalar> AI;
  PtrList<dimensionedScalar> AI_;

  // density factor exponent for component i [-]
  List<scalar> PSII;
  PtrList<dimensionedScalar> PSII_;

  // activation temperature [K]
  List<scalar> EI;
  PtrList<dimensionedScalar> EI_;

  // minimum temperature for reaction i [K]
  List<scalar> TRACI;
  PtrList<dimensionedScalar> TRACI_;

  // composite mixture densities [kg/m3]
  PtrList<volScalarField> rhoI;

  // Temperature [K]
  volScalarField& Ta_;

  // Solid density [K]
  volScalarField& rho_s_;

  // Total pyrolysis rates [kg/m3/s]
  volScalarField& piTotal_;

  //- Virgin solid density [kg/m^3]
  volScalarField& rho_v_;

  //- Charred solid density [kg/m^3]
  volScalarField& rho_c_;

 public:

  //- Runtime type information
  TypeName("FIAT");

  // Constructors

  //- Construct from mesh and phase name
  FIATPyrolysisModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~FIATPyrolysisModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

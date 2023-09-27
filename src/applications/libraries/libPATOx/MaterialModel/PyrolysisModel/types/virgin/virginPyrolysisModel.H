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
    along with OpenFOAM.  If LinearArrheniust, see <http://www.gnu.org/licenses/>.

Class
    Foam::virginPyrolysisModel

Description

SourceFiles
    virginPyrolysisModel.H
    virginPyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef virginPyrolysisModel_H
#define virginPyrolysisModel_H

#include "simplePyrolysisModel.H"
#include "simpleMaterialChemistryModel.H"
#include "simpleEnergyModel.H"
#include "PList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class virginPyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class virginPyrolysisModel
:
public simplePyrolysisModel
{

 protected:

  //- Number of solid phases
  scalar& nSolidPhases_;

  //- Number of pyrolysis reactions per phase
  PList<scalar> nPyroReac;

  //- Total pyrolysis gas production [kg/m3/s]
  volScalarField& piTotal_;

  //- Advancement of the pyrolysis reactions: 1 to 0
  volScalarField& tau;

  //- Solid virgin mass density [kg/m3]
  volScalarField& rho_v_;

  //- Solid charred mass density [kg/m3]
  volScalarField& rho_c_;

  //- Unity volScalarField
  volScalarField& oneVolScalarField;

  //- Number total of the pyrolysis reactions
  int sizeXsi;

  //- Hash table: column 0 = phase number, column 1 = phase reaction number
  RectangularMatrix<label> hashXsi;

  //- Factor F of the pyrolysis reactions
  PList<dimensionedScalar> Fp;

  //- Solid initial volume fraction [-]
  PList<volScalarField> solidEpsI_;

  //- Solid initial mass density [kg/m3]
  PList<volScalarField> solidRhoI_;

  //- Solid volume fraction [-]
  PList<volScalarField> solidEps_;

  //- Solid mass density [kg/m3]
  PList<volScalarField> solidRho_;

 public:

  //- Runtime type information
  TypeName("virgin");

  // Constructors

  //- Construct from mesh and phase name
  virginPyrolysisModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~virginPyrolysisModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Initialize
  virtual void initialize();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

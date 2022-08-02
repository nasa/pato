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
    Foam::FixedPressureGammaMassModel

Description

SourceFiles
    FixedPressureGammaMassModel.H
    FixedPressureGammaMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef FixedPressureGammaMassModel_H
#define FixedPressureGammaMassModel_H

#include "simpleMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FixedPressureGammaMassModel Declaration
\*---------------------------------------------------------------------------*/

class FixedPressureGammaMassModel
:
public simpleMassModel
{

 protected:

  //- Internal pressure [Pa]
  scalar intP_;

  //- Boundary pressure [Pa]
  List<Tuple2<word,scalar>> bcP_;

  //- Pressure [Pa]
  volScalarField& P_;

  //- Internal convective part [s]
  scalar intGamma_;

  //- Boundary convective part [s]
  List<Tuple2<word,scalar>> bcGamma_;

  //- Convective part [s]
  volTensorField Gamma;

  //- Symmetrical of convective part [s]
  volSymmTensorField Gamma_symm;

 public:

  //- Runtime type information
  TypeName("FixedPressureGamma");

  // Constructors

  //- Construct from mesh and phase name
  FixedPressureGammaMassModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~FixedPressureGammaMassModel();

  //- Member Functions

  //- Update thermodynamic and transport properties in the field
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

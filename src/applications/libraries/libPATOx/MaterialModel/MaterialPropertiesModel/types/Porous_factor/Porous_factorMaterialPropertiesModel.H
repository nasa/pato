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
    along with OpenFOAM.  If Porous_factor_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::Porous_polynomial_k_UQMaterialPropertiesModel

Description

SourceFiles
    Porous_factorMaterialPropertiesModel.H
    Porous_factorMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef Porous_factorMaterialPropertiesModel_H
#define Porous_factorMaterialPropertiesModel_H

#include "PorousMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Porous_polynomial_k_UQMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class Porous_factorMaterialPropertiesModel
:
public PorousMaterialPropertiesModel
{

 protected:

  //- Start model initialization after parent class
  const Switch startModelInit_;

  //- Field factors
  List<Tuple2<word,scalar>> fieldFactors_;
  List<scalar> factorList_virgin_;
  List<scalar> factorList_char_;

 public:

  //- Runtime type information
  TypeName("Porous_factor");

  // Constructors

  //- Construct from mesh and phase name
  Porous_factorMaterialPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~Porous_factorMaterialPropertiesModel();

  //- Member Functions

  //- Initialize the model
  virtual void initialize();

  //- Update the model
  virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

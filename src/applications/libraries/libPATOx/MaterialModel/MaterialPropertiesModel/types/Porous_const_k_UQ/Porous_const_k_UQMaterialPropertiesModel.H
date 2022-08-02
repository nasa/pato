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
    along with OpenFOAM.  If Porous_const_k_UQt, see <http://www.gnu.org/licenses/>.

Class
    Foam::Porous_const_k_UQMaterialPropertiesModel

Description

SourceFiles
    Porous_const_k_UQMaterialPropertiesModel.H
    Porous_const_k_UQMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef Porous_const_k_UQMaterialPropertiesModel_H
#define Porous_const_k_UQMaterialPropertiesModel_H

#include "PorousMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Porous_const_k_UQMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class Porous_const_k_UQMaterialPropertiesModel
:
public PorousMaterialPropertiesModel
{

 protected:

  //- Constant virgin thermal conductivity (UQ) [W/m/K]
  scalar k_const_v_;

  //- Constant charred thermal conductivity (UQ) [W/m/K]
  scalar k_const_c_;

 public:

  //- Runtime type information
  TypeName("Porous_const_k_UQ");

  // Constructors

  //- Construct from mesh and phase name
  Porous_const_k_UQMaterialPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~Porous_const_k_UQMaterialPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- update the thermal conductivity
  void updatek();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

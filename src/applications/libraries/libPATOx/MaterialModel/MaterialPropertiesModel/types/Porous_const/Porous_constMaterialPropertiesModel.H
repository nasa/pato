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
    Foam::Porous_constMaterialPropertiesModel

Description

SourceFiles
    Porous_constMaterialPropertiesModel.H
    Porous_constMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef Porous_constMaterialPropertiesModel_H
#define Porous_constMaterialPropertiesModel_H

#include "simpleMaterialPropertiesModel.H"
#include "CommonMaterialPropertiesModel.H"

#define bindFun_Porous_const(fun,name) std::bind(&Foam::Porous_constMaterialPropertiesModel::fun, this, name)
#define newMapField_Porous_const(name,type,init,update) new MapField(name,type,bindFun_Porous_const(init,name),bindFun_Porous_const(update,name))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Porous_constMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class Porous_constMaterialPropertiesModel
:
public simpleMaterialPropertiesModel,
public CommonMaterialPropertiesModel
{
 protected:


 public:

  //- Runtime type information
  TypeName("Porous_const");

  // Constructors

  //- Construct from mesh and phase name
  Porous_constMaterialPropertiesModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~Porous_constMaterialPropertiesModel();

  //- Member Functions

  // Initialize functions

  //- Initialize the porous material properties
  virtual void initialize();

  //- Initialize scalar
  void init_scalar(word name);

  //- Initialize tensor
  void init_tensor(word name);

  //- Initialize tensor
  void init_k(word name);

  // Update functions

  //- Update the porous material properties
  virtual void update();

  //- Empty function
  void update_empty(word name);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

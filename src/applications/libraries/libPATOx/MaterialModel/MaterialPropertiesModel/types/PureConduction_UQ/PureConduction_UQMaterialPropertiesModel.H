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
    along with OpenFOAM.  If PureConduction_UQ_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::PureConduction_UQMaterialPropertiesModel

Description
    PureConduction material properties (k, cp ,rho_s) with a constant value

SourceFiles
    PureConduction_UQMaterialPropertiesModel.H
    PureConduction_UQMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef PureConduction_UQMaterialPropertiesModel_H
#define PureConduction_UQMaterialPropertiesModel_H

#include "simpleMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PureConduction_UQMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class PureConduction_UQMaterialPropertiesModel
:
public simpleMaterialPropertiesModel
{

 protected:

  //- Material Properties directory
  const fileName materialPropertiesDirectory;

  //- Dictionary constantProperties
  IOdictionary constantPropertiesDictionary;

  //- Solid thermal conductivity [W/m/K]
  volTensorField& k_;

  //- Solid heat capacity [J/kg/K]
  volScalarField& cp_;

  //- Total solid density [kg/m3]
  volScalarField& rho_s_;

  //- Constant value of the thermal conductivity [W/m/K]
  scalar kConst_;

  //- Constant value of the heat capacity [J/kg/K]
  scalar cpConst_;

  //- Constant value of the solid mass density [kg/m3]
  scalar rho_sConst_;

 public:

  //- Runtime type information
  TypeName("PureConduction_UQ");

  // Constructors

  //- Construct from mesh and phase name
  PureConduction_UQMaterialPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PureConduction_UQMaterialPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

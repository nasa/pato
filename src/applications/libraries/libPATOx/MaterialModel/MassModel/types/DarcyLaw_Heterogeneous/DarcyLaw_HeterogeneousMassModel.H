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
    Foam::DarcyLaw_HeterogeneousMassModel

Description

SourceFiles
    DarcyLaw_HeterogeneousMassModel.H
    DarcyLaw_HeterogeneousMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef DarcyLaw_HeterogeneousMassModel_H
#define DarcyLaw_HeterogeneousMassModel_H

#include "DarcyLawMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DarcyLaw_HeterogeneousMassModel Declaration
\*---------------------------------------------------------------------------*/

class DarcyLaw_HeterogeneousMassModel
:
public DarcyLawMassModel
{

 protected:

  //- Heterogeneous rate [kg/m3/s] - updated in MaterialChemistryModel
  const volScalarField& omegaHeterogeneousRate_;

 public:

  //- Runtime type information
  TypeName("DarcyLaw_Heterogeneous");

  // Constructors

  //- Construct from mesh and phase name
  DarcyLaw_HeterogeneousMassModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~DarcyLaw_HeterogeneousMassModel();

  //- Member Functions

  //- Update the mass conservation fields
  virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

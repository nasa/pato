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
    Foam::constantLambdaBlowingModel

Description

SourceFiles
    iterativeBprimeBlowingCorrectionModel.H
    iterativeBprimeBlowingCorrectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef iterativeBprimeBlowingCorrectionModel_H
#define iterativeBprimeBlowingCorrectionModel_H

#include "simpleBlowingCorrectionModel.H"
#include "mathFunctions.H"
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class iterativeBprimeBlowingCorrectionModel Declaration
\*---------------------------------------------------------------------------*/

class iterativeBprimeBlowingCorrectionModel
:
public simpleBlowingCorrectionModel
{

 protected:

  //- Energy Model
  simpleEnergyModel& energyModel;

  //- Folder with the previous simulation results of Bprime
  const fileName folderPreviousBprime_;

  //- BC values of the previous simulation results of Bprime
  List<scalarList> valuesPreviousBprime_;

  //- Times of the previous simulation results of Bprime
  scalarList timesPreviousBprime_;

 public:

  //- Runtime type information
  TypeName("iterativeBprime");

  // Constructors

  //- Construct from mesh and phase name
  iterativeBprimeBlowingCorrectionModel(
      const fvMesh& mesh,
      const word& phaseName,
      const word& regionName,
      const label& currentPatchID,
      const dictionary dict
  );

  //- Destructor
  virtual ~iterativeBprimeBlowingCorrectionModel();

  //- Member Functions

  //- Initialize
  virtual void initialize();

  //- Update
  virtual void update();

  //- write in os
  virtual void write(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

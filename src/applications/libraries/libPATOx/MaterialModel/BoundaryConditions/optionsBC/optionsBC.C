/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "optionsBC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optionsBC::optionsBC(const fvMesh& mesh, const dictionary& dict, const label& patchID)
{
  int index = 0;
  wordList list_available =  simpleOptionModel::fvMeshConstructorTablePtr_->sortedToc();
  forAll(list_available, i) {
    Switch option_flag = dict.lookupOrDefault<Switch>(list_available[i],"no");
    if (option_flag) {
      autoPtr<simpleOptionModel> option(simpleOptionModel::New(mesh, dict, patchID, list_available[i]));
      options_list_autoPtr.append(option);
      options_list_ptr.append(&options_list_autoPtr[index]());
      index++;
    }
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optionsBC::update()
{
  forAll(options_list_ptr, i) {
    options_list_ptr[i]->update();
  }
}

void Foam::optionsBC::update(label faceI)
{
  forAll(options_list_ptr, i) {
    options_list_ptr[i]->update(faceI);
  }
}

// ************************************************************************* //

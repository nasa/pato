/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "basicUserThermo.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


inline Foam::volScalarField& Foam::basicUserThermo::meshLookupOrConstructScalar
(
    const fvMesh& mesh,
    const char* name
)
{
  if (!mesh.objectRegistry::foundObject<volScalarField>(name)) {
    volScalarField* fPtr
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volScalarField&>
         (
             mesh.objectRegistry::lookupObject<volScalarField>(name)
         );
}

inline Foam::volScalarField& Foam::basicUserThermo::meshLookupOrConstructScalar
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedScalar& dimScalar_
)
{
  if (!mesh.objectRegistry::foundObject<volScalarField>(name)) {
    volScalarField* fPtr
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimScalar_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volScalarField&>
         (
             mesh.objectRegistry::lookupObject<volScalarField>(name)
         );
}

inline Foam::volScalarField& Foam::basicUserThermo::meshLookupOrConstructScalar
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedScalar& dimScalar_,
    const wordList& BC_
)
{
  if (!mesh.objectRegistry::foundObject<volScalarField>(name)) {
    volScalarField* fPtr
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimScalar_,
            BC_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volScalarField&>
         (
             mesh.objectRegistry::lookupObject<volScalarField>(name)
         );
}

inline Foam::volVectorField& Foam::basicUserThermo::meshLookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name
)
{
  if (!mesh.objectRegistry::foundObject<volVectorField>(name)) {
    volVectorField* fPtr
    (
        new volVectorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volVectorField&>
         (
             mesh.objectRegistry::lookupObject<volVectorField>(name)
         );
}

inline Foam::volVectorField& Foam::basicUserThermo::meshLookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedVector& dimVec_
)
{
  if (!mesh.objectRegistry::foundObject<volVectorField>(name)) {
    volVectorField* fPtr
    (
        new volVectorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimVec_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volVectorField&>
         (
             mesh.objectRegistry::lookupObject<volVectorField>(name)
         );
}

inline Foam::volVectorField& Foam::basicUserThermo::meshLookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedVector& dimVec_,
    const wordList& BC_
)
{
  if (!mesh.objectRegistry::foundObject<volVectorField>(name)) {
    volVectorField* fPtr
    (
        new volVectorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimVec_,
            BC_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volVectorField&>
         (
             mesh.objectRegistry::lookupObject<volVectorField>(name)
         );
}

inline Foam::volTensorField& Foam::basicUserThermo::meshLookupOrConstructTensor
(
    const fvMesh& mesh,
    const char* name
)
{
  if (!mesh.objectRegistry::foundObject<volTensorField>(name)) {
    volTensorField* fPtr
    (
        new volTensorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volTensorField&>
         (
             mesh.objectRegistry::lookupObject<volTensorField>(name)
         );
}

inline Foam::volTensorField& Foam::basicUserThermo::meshLookupOrConstructTensor
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedTensor& dimTens_
)
{
  if (!mesh.objectRegistry::foundObject<volTensorField>(name)) {
    volTensorField* fPtr
    (
        new volTensorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimTens_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volTensorField&>
         (
             mesh.objectRegistry::lookupObject<volTensorField>(name)
         );
}

inline Foam::volTensorField& Foam::basicUserThermo::meshLookupOrConstructTensor
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedTensor& dimTens_,
    const wordList& BC_
)
{
  if (!mesh.objectRegistry::foundObject<volTensorField>(name)) {
    volTensorField* fPtr
    (
        new volTensorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimTens_,
            BC_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volTensorField&>
         (
             mesh.objectRegistry::lookupObject<volTensorField>(name)
         );
}

Foam::wordList Foam::basicUserThermo::eBoundaryTypes()
{
  const volScalarField::Boundary& tbf =
      this->T_.boundaryField();

  wordList hbt = tbf.types();

  forAll(tbf, patchi) {
    if (isA<fixedValueFvPatchScalarField>(tbf[patchi])) {
      hbt[patchi] = fixedValueFvPatchScalarField::typeName;
    } else if
    (
        isA<zeroGradientFvPatchScalarField>(tbf[patchi])
    ) {
      hbt[patchi] = fixedGradientFvPatchScalarField::typeName;
    } else if (isA<fixedGradientFvPatchScalarField>(tbf[patchi])) {
      hbt[patchi] = fixedGradientFvPatchScalarField::typeName;
    } else if (hbt[patchi]=="empty") {
      hbt[patchi] = "empty";
    } else {
      FatalErrorInFunction << tbf[patchi] << "temperature boundary type not implemented yet." << exit(FatalError);
    }
  }

  return hbt;
}

// ************************************************************************* //


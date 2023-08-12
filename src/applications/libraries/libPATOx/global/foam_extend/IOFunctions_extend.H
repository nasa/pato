/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam_extend_

Description
    IO functions

Details

\*---------------------------------------------------------------------------*/

#if defined(FOAM_EXTEND)
#ifndef IOFunctionsFoamExtend_H
#define IOFunctionsFoamExtend_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define namespaceFoam
#undef FALSE
#undef TRUE
#define namespaceFoam_extend_
#include "declare_macros_foam_extend.H" // Redeclare Foam Extend 4.1 macros
#include "Foam_extend_fvCFD.H"

namespace Foam_extend_
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline volScalarField& meshLookupScalar
(
    const fvMesh& mesh,
    const char* name,
    const word modelName_
)
{
  if (!mesh.objectRegistry::foundObject<volScalarField>(name)) {
    FatalErrorInFunction << name << " not found in the mesh database. Please choose another " << modelName_ << " solver." << exit(FatalError);
  }

  return const_cast<volScalarField&>
         (
             mesh.objectRegistry::lookupObject<volScalarField>(name)
         );
}

inline volVectorField& meshLookupVector
(
    const fvMesh& mesh,
    const char* name,
    const word modelName_
)
{
  if (!mesh.objectRegistry::foundObject<volVectorField>(name)) {
    FatalErrorInFunction << name << " not found in the mesh database. Please choose another " << modelName_ << " solver." << exit(FatalError);
  }

  return const_cast<volVectorField&>
         (
             mesh.objectRegistry::lookupObject<volVectorField>(name)
         );
}

inline volTensorField& meshLookupTensor
(
    const fvMesh& mesh,
    const char* name,
    const word modelName_
)
{
  if (!mesh.objectRegistry::foundObject<volTensorField>(name)) {
    FatalErrorInFunction << name << " not found in the mesh database. Please choose another " << modelName_ << " solver." << exit(FatalError);
  }

  return const_cast<volTensorField&>
         (
             mesh.objectRegistry::lookupObject<volTensorField>(name)
         );
}


template<typename Type>
inline Type& meshLookupOrConstruct
(
    const fvMesh& mesh,
    const char* name
)
{
  if (!mesh.objectRegistry::foundObject<Type>(name)) {
    Type* fPtr
    (
        new Type
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

  return const_cast<Type&>
         (
             mesh.objectRegistry::lookupObject<Type>(name)
         );
}

template<typename Type1, typename Type2>
inline Type1& meshLookupOrConstruct
(
    const fvMesh& mesh,
    const char* name,
    const dimensioned<Type2>& dim_,
    const IOobject::readOption& readOpt_ = IOobject::READ_IF_PRESENT
)
{
  if (!mesh.objectRegistry::foundObject<Type1>(name)) {
    Type1* fPtr
    (
        new Type1
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                readOpt_,
                IOobject::NO_WRITE
            ),
            mesh,
            dim_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<Type1&>
         (
             mesh.objectRegistry::lookupObject<Type1>(name)
         );
}

template<typename Type1, typename Type2>
inline Type1& meshLookupOrConstruct
(
    const fvMesh& mesh,
    const char* name,
    const dimensioned<Type2>& dim_,
    const word& BC_,
    const IOobject::readOption& readOpt_ = IOobject::READ_IF_PRESENT
)
{
  if (!mesh.objectRegistry::foundObject<Type1>(name)) {
    Type1* fPtr
    (
        new Type1
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                readOpt_,
                IOobject::NO_WRITE
            ),
            mesh,
            dim_,
            BC_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<Type1&>
         (
             mesh.objectRegistry::lookupObject<Type1>(name)
         );
}

inline volScalarField& meshLookupOrConstructScalar
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

inline volScalarField& meshLookupOrConstructScalar
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

inline volScalarField& meshLookupOrConstructScalar
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedScalar& dimScalar_,
    const word& BC_
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



inline volVectorField& meshLookupOrConstructVector
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

inline volVectorField& meshLookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name,
    const dimensionedVector dimVec_
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
inline volTensorField& meshLookupOrConstructTensor
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

inline volTensorField& meshLookupOrConstructTensor
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

inline volSymmTensorField& meshLookupOrConstructSymmTensor
(
    const fvMesh& mesh,
    const char* name,
    const dimensioned<symmTensor>& dimSymmTens_
)
{
  if (!mesh.objectRegistry::foundObject<volSymmTensorField>(name)) {
    volSymmTensorField* fPtr
    (
        new volSymmTensorField
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
            dimSymmTens_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<volSymmTensorField&>
         (
             mesh.objectRegistry::lookupObject<volSymmTensorField>(name)
         );
}

inline surfaceVectorField& meshLookupOrConstructSurfaceVector
(
    const fvMesh& mesh,
    const char* name
)
{
  if (!mesh.objectRegistry::foundObject<surfaceVectorField>(name)) {
    surfaceVectorField* fPtr
    (
        new surfaceVectorField
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

  return const_cast<surfaceVectorField&>
         (
             mesh.objectRegistry::lookupObject<surfaceVectorField>(name)
         );
}

//template<class Type>
inline surfaceVectorField& meshLookupOrConstructSurfaceVector
(
    const fvMesh& mesh,
    const char* name,
    const surfaceVectorField& svf
    //tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> linearInterpolate_
//    const GeometricField<Type, fvsPatchField, surfaceMesh>& linearInterpolate_
)
{
  if (!mesh.objectRegistry::foundObject<surfaceVectorField>(name)) {
    surfaceVectorField* fPtr
    (
        new surfaceVectorField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            svf
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<surfaceVectorField&>
         (
             mesh.objectRegistry::lookupObject<surfaceVectorField>(name)
         );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam_extend_

#include "declare_macros_openfoam.H" // Redeclare OpenFOAM 7 macros

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
#endif

// ************************************************************************* //

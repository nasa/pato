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

#include "GradPTimeControlModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GradPTimeControlModel::GradPTimeControlModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleTimeControlModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
massModel_(meshLookupOrConstructModel<simpleMassModel>(mesh_,dictName_,"Mass")),
T_(meshLookupOrConstructScalar(mesh,"Ta")),
vG_(massModel_.vG()),
vS_
(
    IOobject
    (
        "vS",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(vG_)
),
adjustTimeStep
(
    const_cast<Time&>(mesh_.time()).controlDict().lookup("adjustTimeStep")
),
stitchingSolidPatchName
(
    const_cast<Time&>(mesh_.time()).controlDict().lookupOrDefault<word>("stitchingSolidPatch", "undefined")
),
maxCo
(
    readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("maxCo"))
),
adjustStartTime
(
    readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("adjustStartTime"))
),
minDeltaTw_(0) // NAmrofel 26/11/2019
{

  maxDeltaT = GREAT;

  if (const_cast<Time&>(mesh_.time()).controlDict().found("maxDeltaT")) {
    maxDeltaT = readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("maxDeltaT"));
  }

  REVlength = 1e-3;

  if (const_cast<Time&>(mesh_.time()).controlDict().found("REVlength")) {
    REVlength = readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("REVlength"));
  }

  minDeltaT = SMALL;

  if (const_cast<Time&>(mesh_.time()).controlDict().found("minDeltaT")) {
    minDeltaT = readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("minDeltaT"));
  }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GradPTimeControlModel::~GradPTimeControlModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::GradPTimeControlModel::update()
{
  /**** Set "CFL" number ****/
  if(this->debug_) {
    Info << "start -- Foam::GradPTimeControlModel::update()" << endl;
  }
  scalar CoNum = 0.0;
  vS_ = linearInterpolate(vG_);

  if (mesh_.nInternalFaces()) {
    tmp<surfaceScalarField> SfUfbyDelta_tmp =
        mesh_.surfaceInterpolation::deltaCoeffs() * mag(vS_) / (REVlength);
    surfaceScalarField& SfUfbyDelta = const_cast<surfaceScalarField&>(SfUfbyDelta_tmp());
    CoNum = max(SfUfbyDelta).value() * const_cast<Time&>(mesh_.time()).deltaT().value();
  }

  // ***** Adjust time step ***** //

  if (adjustTimeStep && (const_cast<Time&>(mesh_.time()).value() > adjustStartTime)) {
    scalar maxDeltaTFact = maxCo / (CoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.2 * maxDeltaTFact), 1.05);
    scalar deltaTphys =
        max
        (
            min
            (
                deltaTFact * const_cast<Time&>(mesh_.time()).deltaT().value(),
                maxDeltaT
            ),
            minDeltaT
        );
    const_cast<Time&>(mesh_.time()).setDeltaT(deltaTphys);
  }
  if(this->debug_) {
    Info << "end -- Foam::GradPTimeControlModel::update()" << endl;
  }
}

scalar Foam::GradPTimeControlModel::updateMinDeltaTw() // NAmrofel 26/11/2019
{
  // ***** Calculate 2T criteria ***** //
 if (const_cast<Time&>(mesh_.time()).controlDict().found("stitchingSolidPatch"))
{
  Tws_n = T_.boundaryField()[mesh_.boundaryMesh().findPatchID(stitchingSolidPatchName)];
  Tws_nless1 = T_.oldTime().boundaryField()[mesh_.boundaryMesh().findPatchID(stitchingSolidPatchName)];

  deltaTws = mag( (Tws_n - Tws_nless1) / (mesh_.time()).deltaT().value() );
  minDeltaTw_ = min(deltaTws); 

  Info << "minDeltaTw" << minDeltaTw_ << "  K/s";
}
  return minDeltaTw_;
}

// ************************************************************************* //

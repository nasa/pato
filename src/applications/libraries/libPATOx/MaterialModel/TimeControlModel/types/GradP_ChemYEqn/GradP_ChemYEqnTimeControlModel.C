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

#include "GradP_ChemYEqnTimeControlModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GradP_ChemYEqnTimeControlModel::GradP_ChemYEqnTimeControlModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleTimeControlModel(mesh, regionName),
massModel_(refModel<simpleMassModel>()),
vG_(massModel_.refVolField<vector>("vG")),
vS_(createSurfaceField<vector>("vS",linearInterpolate(vG_))),
MaterialChemistryModel_(refModel<simpleMaterialChemistryModel>()),
Y_(MaterialChemistryModel_.massFractions()),
Yold_(MaterialChemistryModel_.oldMassFractions()),
speciesNames_(MaterialChemistryModel_.speciesNames()),
dtChem_(MaterialChemistryModel_.dtChem()),
adjustTimeStep
(
    const_cast<Time&>(mesh_.time()).controlDict().lookup("adjustTimeStep")
),
maxCo
(
    readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("maxCo"))
),
adjustStartTime
(
    readScalar(const_cast<Time&>(mesh_.time()).controlDict().lookup("adjustStartTime"))
),
chemTransEulerStepLimiter(simpleTimeControlModel::materialDict_.subDict("TimeControl").template lookupOrDefault<Switch>("chemTransEulerStepLimiter","no"))
{
  // Extract indices of C(gr) for in-depth oxidation of the solid
  iCs_ = -1;
  for (int i = 0; i < Y_.size(); i++) {
    if (speciesNames_[i] == "C(gr)") {
      iCs_ = i;
    }
  }

  if (iCs_<0) {
    FatalErrorInFunction << "C(gr) not found in the mixture." << exit(FatalError);
  }

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

  // set tolerance on mass fraction variation rate
  Ythreshold = 1e-3;
  if (mesh_.time().controlDict().found("Ythreshold")) {
    Ythreshold = readScalar(mesh_.time().controlDict().lookup("Ythreshold"));
  }

  dYtolMin = 0.01;
  if (mesh_.time().controlDict().found("dYtolMin")) {
    dYtolMin = readScalar(mesh_.time().controlDict().lookup("dYtolMin"));
  }

  dYtolMax = 0.05;
  if (mesh_.time().controlDict().found("dYtolMax")) {
    dYtolMax = readScalar(mesh_.time().controlDict().lookup("dYtolMax"));
  }

  minChemDeltaT = SMALL;
  if (mesh_.time().controlDict().found("minChemDeltaT")) {
    minChemDeltaT = readScalar(mesh_.time().controlDict().lookup("minChemDeltaT"));
  }

  dtChemYEqn =0;
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GradP_ChemYEqnTimeControlModel::~GradP_ChemYEqnTimeControlModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GradP_ChemYEqnTimeControlModel::update()
{
  /*** Compute chem deltaT  ***/
  tmp<volScalarField> Yini_tmp = 0.0 * Y_[0];
  volScalarField& Yini = const_cast<volScalarField&>(Yini_tmp());
  scalar dYmax = 0;
  scalar speciesMax = 0;
  scalar cellMax = 0;

  for (label i = 0; i < Y_.size(); i++) {
    if (i != iCs_) {
      Yini = Yold_[i];
      forAll(Yini, celli) {
        if (Yini[celli] > Ythreshold) {
          if (fabs((Yini[celli] - Y_[i][celli]) / Yini[celli])  > dYmax) {
            dYmax = fabs((Yini[celli] - Y_[i][celli]) / Yini[celli]);
            speciesMax = i;
            cellMax = celli;
          }
        }
      }
    }
  }

  Info << "dYmax(species: " << Y_[speciesMax].name()
       << "; in cell: " << cellMax
       << ") = " << dYmax
       << " | Y = " << Y_[speciesMax][cellMax]
       << nl;

  Info << "dYmax " << dYmax << nl;
  if (dYmax > dYtolMax) {
    dtChemYEqn = mesh_.time().deltaT().value() / 1.2;
  }

  if (dYmax < dYtolMin) {
    dtChemYEqn = mesh_.time().deltaT().value() * 1.2;
  }

  // retain the min chemical time step bewteen chemistry (dtChem) and transport (dtChemYEqn)
  if (chemTransEulerStepLimiter) {
    dtChemYEqn = min(dtChemYEqn, dtChem_);
  }

  /**** Set "CFL" number ****/

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
    const_cast<Time&>(mesh_.time()).setDeltaT(max(min(dtChemYEqn, deltaTphys), minChemDeltaT));

  }
}

scalar Foam::GradP_ChemYEqnTimeControlModel::updateMinDeltaTw()
{
  return 0;
}

// ************************************************************************* //

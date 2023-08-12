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
    along with OpenFOAM.  If FIATt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FIATPyrolysisModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FIATPyrolysisModel::FIATPyrolysisModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simplePyrolysisModel(mesh, dictName),
mesh_(mesh),
GAMMA(readScalar(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("GAMMA"))),
PHI(readScalar(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("PHI"))),
RVI(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("RVI")),
RCI(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("RCI")),
AI(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("AI")),
PSII(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("PSII")),
EI(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("EI")),
TRACI(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").lookup("TRACI")),
piTotal_(createVolField<scalar>("piTotal",dimensionedScalar("0",dimMass/dimVolume/dimTime,0))),
rho_v_(createVolField<scalar>("rho_v",dimensionedScalar("0",dimMass/dimVolume,0))),
rho_c_(createVolField<scalar>("rho_c",dimensionedScalar("0",dimMass/dimVolume,0))),
energyModel(refModel<simpleEnergyModel>()),
Ta_(energyModel.refVolField<scalar>("Ta")),
rho_s_(energyModel.refVolField<scalar>("rho_s"))
{
  int numberComponents = 3;
  if (RVI.size()!=numberComponents || RCI.size()!=numberComponents
      || AI.size()!=numberComponents  || PSII.size()!=numberComponents
      || EI.size()!=numberComponents  || TRACI.size()!=numberComponents) {
    FatalErrorInFunction
        <<  "Size of RVI, RCI, AI, PSII, EI or TRACI is different to " << numberComponents << " in \"" << simplePyrolysisModel::name() << "\" file." << nl
        << "Only old model of FIAT pyrolysis is implemented."
        << exit(FatalError);
  }

  RVI_.resize(numberComponents);
  RCI_.resize(numberComponents);
  AI_.resize(numberComponents);
  PSII_.resize(numberComponents);
  EI_.resize(numberComponents);
  TRACI_.resize(numberComponents);

  forAll(RVI, i) {
    RVI_.set
    (
        i,
        new dimensionedScalar
        (
            "RVI"+std::to_string(i), dimensionSet(1,-3,0,0,0), RVI[i]
        )
    );

    RCI_.set
    (
        i,
        new dimensionedScalar
        (
            "RCI"+std::to_string(i), dimensionSet(1,-3,0,0,0), RCI[i]
        )
    );

    AI_.set
    (
        i,
        new dimensionedScalar
        (
            "AI"+std::to_string(i), dimensionSet(0,0,-1,0,0), AI[i]
        )
    );


    PSII_.set
    (
        i,
        new dimensionedScalar
        (
            "PSII"+std::to_string(i), dimless, PSII[i]
        )
    );

    EI_.set
    (
        i,
        new dimensionedScalar
        (
            "EI"+std::to_string(i), dimensionSet(0,0,0,1,0), EI[i]
        )
    );
  }

  wordList BCgradZero(mesh.boundaryMesh().size());
  forAll(BCgradZero, BCI) {
    BCgradZero[BCI]= "zeroGradient";
  }

  rhoI.resize(numberComponents);
  forAll(rhoI, rhoII) {

    rhoI.set
    (
        rhoII,
        new volScalarField
        (
            IOobject
            (
                "rhoI[" + std::to_string(rhoII) + "]",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            RVI_[rhoII],
            BCgradZero
        )
    );
  }

  dimensionedScalar rho_unity("rho_unity", dimMass/dimVolume, 1.0);
  rho_v_ = (1-PHI)*(GAMMA*(RVI[0]+RVI[1]) + (1-GAMMA)*RVI[2])*rho_unity;
  rho_c_ = (1-PHI)*(GAMMA*(RCI[0]+RCI[1]) + (1-GAMMA)*RCI[2])*rho_unity;

  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FIATPyrolysisModel::~FIATPyrolysisModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FIATPyrolysisModel::update()
{

  forAll(rhoI, rhoII) {
    volScalarField rhoI_old = rhoI[rhoII];

    if(!this->dynamicMesh_) {
      solve
      (
          fvm::ddt(rhoI[rhoII])
          + AI_[rhoII]* RVI_[rhoII] * pow((rhoI[rhoII] - RCI_[rhoII])/RVI_[rhoII] , PSII_[rhoII]) * exp(-EI_[rhoII]/Ta_),
          "Xsii"
      );
    } else {
      solve
      (
          fvm::ddt(rhoI[rhoII]) - fvm::div(mesh_.phi(), rhoI[rhoII])
          + AI_[rhoII]* RVI_[rhoII] * pow((rhoI[rhoII] - RCI_[rhoII])/RVI_[rhoII] , PSII_[rhoII]) * exp(-EI_[rhoII]/Ta_),
          "Xsii"
      );
    }
    forAll(Ta_, cellI) {
      if(Ta_[cellI]<TRACI[rhoII]) {
        rhoI[rhoII][cellI]=rhoI_old[cellI];
      }
    }
  }

  piTotal_ = (1-PHI)*GAMMA*(- fvc::ddt(rhoI[0]) - fvc::ddt(rhoI[1]));
  rho_s_ = (1-PHI)*(GAMMA*(rhoI[0]+rhoI[1]) + (1-GAMMA)*rhoI[2]);
}

void Foam::FIATPyrolysisModel::initialize()
{}

// ************************************************************************* //

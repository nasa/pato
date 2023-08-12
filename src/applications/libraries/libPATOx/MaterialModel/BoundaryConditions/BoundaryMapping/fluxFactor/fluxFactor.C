/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fluxFactor.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFactor::fluxFactor
(
    const dictionary& dict,
    const simpleModel& simpModel,
    const word& fluxMapName,
    const word& pressureMapName
)
  :
dict_(dict),
simpleModel_(const_cast<simpleModel&>(simpModel)),
mesh_(simpleModel_.mesh()),
indexDistance_i_(0),
indexDistance_iPlus1_(0),
f_D_(0.0),
value_(0.0),
currentDistanceIndex_(0.0),
fluxMap_(simpleModel_.createVolFieldIfNotFound<scalar>(simpleModel_,fluxMapName.c_str(), dimensionedScalar("0", dimless, scalar(0.0)))),
pressureMap_(simpleModel_.createVolFieldIfNotFound<scalar>(simpleModel_,pressureMapName.c_str(), dimensionedScalar("0", dimless, scalar(0.0)))),
fluxFactorMapFileName_(dict_.lookup("fluxFactorMapFileName")),
fluxFactorProjection_(dict_.lookupOrDefault<Switch>("fluxFactorProjection","no")),
fluxFactorCenter_(dict_.lookup("fluxFactorCenter")),
fluxFactorNormal_(dict_.lookup("fluxFactorNormal")),
fluxFactorThreshold_(dict_.lookup("fluxFactorThreshold")),
pointMotionDirection_(dict_.lookup("pointMotionDirection")),
moreThanThresholdFlag_(dict_.lookupOrDefault<Switch>("moreThanThresholdFlag",true))
{
  // read and store the boundary layer properties into the RAM for faster access
  Info << simpleModel::getTabLevel() << "Reading flux Map" << nl;

  // OFstream fluxFactorMapOutputFile ("fluxFactorMapOutputFile");  // opens an output file
  IFstream fluxFactorMap(changeEnviVar(fluxFactorMapFileName_));  // opens an input file
  IFstream fluxFactorMapTemp(changeEnviVar(fluxFactorMapFileName_));  // opens a temporary input file

  if (fluxFactorMap.good() == false) { // checks that the input file is opened
    FatalErrorInFunction << fluxFactorMapFileName_ << " not found" << exit(FatalError);
  }

  Info << simpleModel::getTabLevel() << "The fluxFactorMap must have 3 columns: "
       "mapping variable (e.g. Y), flux factor, pressure factor"
       << nl;

  int columnTableFF = 3;
  int rawTableFF = 0;
  scalar tempFF, rawTableFracFF, rawTableIntFF;
  int i_rawFF = 0;

  while (true) {
    fluxFactorMapTemp >> tempFF;
    if (fluxFactorMapTemp.eof() == 1) {
      break;
    }
    i_rawFF++;
  }

  rawTableFracFF = modf(static_cast<scalar>(i_rawFF) / static_cast<scalar>(columnTableFF), &rawTableIntFF);
  if (rawTableFracFF != 0) {
    Info << fluxFactorMapFileName_ << " does not seem to have 3 columns." << nl;
    FatalErrorInFunction << rawTableIntFF << " lines and "
                         << rawTableFracFF * columnTableFF
                         << " 'double precision' have been read."
                         <<  exit(FatalError);
  } else {
    rawTableFF = i_rawFF / columnTableFF;
    Info << simpleModel::getTabLevel() << "The fluxFactorMap has been read (" << rawTableFF << " lines)." << nl;
  }

//
//    RectangularMatrix<scalar> fluxFactorTable(rawTableFF, columnTableFF);
  fluxFactorTable_.setSize(rawTableFF, columnTableFF);
  for (int x = 0; x < rawTableFF; x++) {
    for (int i = 0; i < columnTableFF; i++) {
      fluxFactorMap >> fluxFactorTable_[x][i];
    }
  }

  // the call format is fluxFactorAxi.value(distance);
  //  Info << "factor1 = " << fluxFactorAxi.value(0.03) << nl;
  //  Info << "The flux mapping factors have been read." << nl;
  //  Info << "Process killed after reading the flux factor map." << nl; return (0);  // exits without error
  // END OF THE IO TEST

  // Compute the flux factor (or transfer coefficient factor) for the faces of the ablating suface for 2D and 3D axi-symmetrical heat fluxes
  forAll(mesh_.boundaryMesh(), patchI) {
    tmp<vectorField> c_tmp = mesh_.Cf().boundaryField()[patchI] - fluxFactorCenter_;
    vectorField& c = const_cast<vectorField&>(c_tmp());

    tmp<scalarField> Distance_tmp = (c & fluxFactorNormal_);
    scalarField& Distance = const_cast<scalarField&>(Distance_tmp());

    vector vectorX(1, 0, 0);
    vector vectorY(0, 1, 0);
    vector vectorZ(0, 0, 1);
    if (fluxFactorProjection_) {
      // distance to center, on the plane perpendicular to fluxFactorDirection
      // (assuming that fluxFactorNormal is either vectorX, vectorY, or vectorZ,
      // as required in the input file).
      Distance =
          sqrt
          (
              (c & vectorX) * (c & vectorX)
              + (c & vectorY) * (c & vectorY)
              + (c & vectorZ) * (c & vectorZ)
              - Distance * Distance
          );
    }

    forAll(fluxMap_.boundaryField()[patchI], faceI) {
      fluxMap_.boundaryFieldRef()[patchI][faceI] =
          fluxValue(Distance[faceI]);
      pressureMap_.boundaryFieldRef()[patchI][faceI] =
          pressureValue(Distance[faceI]);
    }
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFactor::~fluxFactor()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// indexDi: Returns the Distance interpolation factor and updates indexDi_ and indexDiPlus1_
void Foam::fluxFactor::update(scalar Distance_)
{
  if ( (Distance_ < fluxFactorTable_[0][0]) | (Distance_ > fluxFactorTable_[fluxFactorTable_.m() - 1][0])) {
    Info << "Distance out of fluxFactor table range" << nl;
    if (Distance_ < fluxFactorTable_[0][0]) {
      Info << "Distance = "
           << Distance_
           << " < "
           << fluxFactorTable_[0][0]
           << nl;
    }

    if (Distance_ > fluxFactorTable_[fluxFactorTable_.m() - 1][0]) {
      Info << "Distance = "
           << Distance_
           << " > "
           << fluxFactorTable_[fluxFactorTable_.m() - 1][0]
           << nl;
    }

    ::exit(1);
  }

  int i_Distance(0);
  while (fluxFactorTable_[i_Distance][0] < Distance_) {
    i_Distance++;
  }

  indexDistance_i_ = i_Distance;

  if (Distance_ == fluxFactorTable_[indexDistance_i_][0]) {
    indexDistance_iPlus1_ = indexDistance_i_;
    f_D_ = 1;
  } else {
    indexDistance_iPlus1_ = indexDistance_i_ - 1;
    f_D_ =
        (
            fluxFactorTable_[indexDistance_iPlus1_][0] -
            Distance_
        ) /
        (
            fluxFactorTable_[indexDistance_iPlus1_][0] -
            fluxFactorTable_[indexDistance_i_][0]
        );
  }

  // return f_D_;
  currentDistanceIndex_ = Distance_;
}

Foam::scalar Foam::fluxFactor::fluxValue(scalar Distance_)
{
  if (currentDistanceIndex_ != Distance_) {
    update(Distance_);
  }

  return ffactor();
}

Foam::scalar Foam::fluxFactor::pressureValue(scalar Distance_)
{
  if (currentDistanceIndex_ != Distance_) {
    update(Distance_);
  }

  return pfactor();
}

const Foam::fileName& Foam::fluxFactor::fluxFactorMapFileName() const
{
  return fluxFactorMapFileName_;
}

const Foam::dimensionedScalar& Foam::fluxFactor::fluxFactorThreshold() const
{
  return fluxFactorThreshold_;
}

const Foam::vector& Foam::fluxFactor::pointMotionDirection() const
{
  return pointMotionDirection_;
}

Foam::volScalarField& Foam::fluxFactor::fluxMap()
{
  return fluxMap_;
}

Foam::volScalarField& Foam::fluxFactor::pressureMap()
{
  return pressureMap_;
}

const Foam::Switch& Foam::fluxFactor::fluxFactorProjection() const
{
  return fluxFactorProjection_;
}

const Foam::vector& Foam::fluxFactor::fluxFactorCenter() const
{
  return fluxFactorCenter_;
}

const Foam::vector& Foam::fluxFactor::fluxFactorNormal() const
{
  return fluxFactorNormal_;
}

const Foam::Switch& Foam::fluxFactor::moreThanThresholdFlag() const
{
  return moreThanThresholdFlag_;
}

// * * * * * * * * * * * * * * * END Member Functions  * * * * * * * * * * * * * //

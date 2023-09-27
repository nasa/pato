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

Class
    Foam::fluxFactor

Description
    Flux factor map reader class
    Extracts/interpolates the flux factor as a function of the distance to a reference point (based on a 2D matrix: distance, factor)

\*---------------------------------------------------------------------------*/

#ifndef fluxFactor_H
#define fluxFactor_H

#include "scalar.H"
#include "RectangularMatrix.H"
#include "fvCFD.H"
#include "IOFunctions.H"
#include "simpleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            fluxFactor class Declaration
\*---------------------------------------------------------------------------*/

class fluxFactor
{

// private
  const dictionary& dict_;
  simpleModel& simpleModel_;
  const fvMesh& mesh_;
  int indexDistance_i_;                       // Index of distance inf: distance_i >= t
  int indexDistance_iPlus1_;                  // Index of distance sup: distance_i <= t
  scalar f_D_;                                // Interpolation coefficient
  scalar value_;
  scalar currentDistanceIndex_;
  RectangularMatrix<scalar> fluxFactorTable_; // 2 columns: distance, factor
  volScalarField& fluxMap_;
  volScalarField& pressureMap_;
  fileName fluxFactorMapFileName_;
  Switch fluxFactorProjection_;
  vector fluxFactorCenter_;
  vector fluxFactorNormal_;
  dimensionedScalar fluxFactorThreshold_;
  vector pointMotionDirection_;
  Switch moreThanThresholdFlag_; // if moreThanThreshold == true => in B' bc: if (fluxMap > fluxFactorThreshold):  normal_ = (1 / (nf & pointMotionDirection_)) * pointMotionDirection_;

 public:

  // Constructors
  fluxFactor();

  //- Construct from components
  fluxFactor(const dictionary& dict, const simpleModel& simpModel, const word& fluxMapName="fluxMap", const word& pressureMapName="pressureMap");


  // Destructor

  virtual ~fluxFactor();

  // Access Functions
  int indexDistance_i() const {
    return indexDistance_i_;
  }
  int indexDistance_iPlus1() const {
    return indexDistance_iPlus1_;
  }
  double f_D() const {
    return f_D_;
  }
  scalar fvalue_i() const {
    return fluxFactorTable_[indexDistance_i_][1];
  }
  scalar fvalue_iPlus1() const {
    return fluxFactorTable_[indexDistance_iPlus1_][1];
  }
  scalar ffactor() const {
    return fvalue_i() * f_D() + fvalue_iPlus1() * (1 - f_D());
  }

  scalar pvalue_i() const {
    return fluxFactorTable_[indexDistance_i_][2];
  }
  scalar pvalue_iPlus1() const {
    return fluxFactorTable_[indexDistance_iPlus1_][2];
  }
  scalar pfactor() const {
    return pvalue_i() * f_D() + pvalue_iPlus1() * (1 - f_D());
  }

  // Member Functions
  // returns the Distance interpolation factor and updates the indexes of Di and Di+1
  void update(scalar Distance_);
  scalar fluxValue(scalar Distance_);
  scalar pressureValue(scalar Distance_);
  const fileName& fluxFactorMapFileName() const;
  const dimensionedScalar& fluxFactorThreshold() const;
  volScalarField& fluxMap();
  volScalarField& pressureMap();
  const vector& pointMotionDirection() const;
  const Switch& fluxFactorProjection() const;
  const vector& fluxFactorCenter() const;
  const vector& fluxFactorNormal() const;
  const Switch& moreThanThresholdFlag() const;

}; // End class MotionInput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

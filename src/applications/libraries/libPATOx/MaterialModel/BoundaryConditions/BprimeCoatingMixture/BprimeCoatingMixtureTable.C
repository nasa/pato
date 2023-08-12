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

#include "BprimeCoatingMixtureTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BprimeCoatingMixtureTable::BprimeCoatingMixtureTable
(
    RectangularMatrix<scalar>& BprimeCoatingMixtureTable
)
  :
BprimeCoatingMixtureTable_(BprimeCoatingMixtureTable),
indexPiMin_(0),
indexPiMax_(0),
indexPiPlus1Min_(0),
indexPiPlus1Max_(0),
f_p_(0.0),
indexBpGiMin_(0),
indexBpGiMax_(0),
indexBpGiPlus1Min_(0),
indexBpGiPlus1Max_(0),
f_BpG_(0.0),
f_BpGi_(0.0),
f_BpGiPlus1_(0.0),
indexTi_(0.0),
indexTiPlus1_(0.0),
f_T_(0.0),
BpC_(0.0),
hw_(0.0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BprimeCoatingMixtureTable::~BprimeCoatingMixtureTable()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// indexTi: Returns the Temperature interpolation factor and updates indexTi_ and indexTiPlus1_
Foam::scalar Foam::BprimeCoatingMixtureTable::indexTi(scalar T_, int i_min_, int i_max_)
{
  if (BprimeCoatingMixtureTable_[i_min_][3] >= BprimeCoatingMixtureTable_[i_max_][3]) { // decreasing temperature table
    if
    (
        (T_ > BprimeCoatingMixtureTable_[i_min_][3]) ||
        (T_ < BprimeCoatingMixtureTable_[i_max_][3])
    ) {
//      Info << "Temperature out of B' table range" << nl;

      if (T_ > BprimeCoatingMixtureTable_[i_min_][3]) {
//        Info << "T = " << T_ << " > " << BprimeCoatingMixtureTable_[i_min_][3] << nl;
        T_ = BprimeCoatingMixtureTable_[i_min_][3];

      }
      if (T_ < BprimeCoatingMixtureTable_[i_max_][3]) {
//        Info << "T = " << T_ << " < " << BprimeCoatingMixtureTable_[i_max_][3] << nl;
        T_ = BprimeCoatingMixtureTable_[i_max_][3];

      }
//            ::exit(1);
    }

    int i_Ti(i_min_);
    while (BprimeCoatingMixtureTable_[i_Ti][3] > T_) {
      i_Ti++;
    }
    indexTi_ = i_Ti;
    if (T_ == BprimeCoatingMixtureTable_[indexTi_][3]) {
      indexTiPlus1_ = indexTi_;
      f_T_ = 1;
      return f_T_;
    } else {
      indexTiPlus1_ = indexTi_ - 1;
      f_T_ = (T_ - BprimeCoatingMixtureTable_[indexTi_][3]) / (BprimeCoatingMixtureTable_[indexTiPlus1_][3] - BprimeCoatingMixtureTable_[indexTi_][3]);
    }
  } // end decreasing temperature table

  else { // increasing temperature table
    if
    (
        (T_ > BprimeCoatingMixtureTable_[i_max_][3]) ||
        (T_ < BprimeCoatingMixtureTable_[i_min_][3])
    ) {
//      Info << "Temperature out of B' table range" << nl;

      if (T_ > BprimeCoatingMixtureTable_[i_max_][3]) {
//        Info << "T = " << T_ << " > " << BprimeCoatingMixtureTable_[i_max_][3] << nl;
        T_ = BprimeCoatingMixtureTable_[i_max_][3];

      }

      if (T_ < BprimeCoatingMixtureTable_[i_min_][3]) {
//        Info << "T = " << T_ << " < " << BprimeCoatingMixtureTable_[i_min_][3] << nl;
        T_ = BprimeCoatingMixtureTable_[i_min_][3];

      }
//            ::exit(1);
    }

    int i_Ti(i_min_);
    while (BprimeCoatingMixtureTable_[i_Ti][3] < T_) {
      i_Ti++;
    }

    indexTiPlus1_ = i_Ti;
    if (T_ == BprimeCoatingMixtureTable_[indexTiPlus1_][3]) {
      indexTi_ = indexTiPlus1_;
      f_T_ = 1;
      return f_T_;
    } else {
      indexTi_ = indexTiPlus1_ - 1;
      f_T_ = (T_ - BprimeCoatingMixtureTable_[indexTi_][3]) / (BprimeCoatingMixtureTable_[indexTiPlus1_][3] - BprimeCoatingMixtureTable_[indexTi_][3]);
    }
  } // end increasing temperature table

// Info << "f_T = " << f_T_ << nl;
  return f_T_;
}

// indexPi: Returns the interpolation factor f and updates IndexPiMin_, IndexPiMax_, IndexPiPlus1Min_, IndexPiPlus1Max_
Foam::scalar Foam::BprimeCoatingMixtureTable::indexPi(scalar p_)
{
  int i_Pi(0);

  if
  (
      (p_ > BprimeCoatingMixtureTable_[0][0]) ||
      (p_ < BprimeCoatingMixtureTable_[BprimeCoatingMixtureTable_.m() - 1][0])
  ) {
//    Info << "Pressure out of B' table range" << nl;

    if (p_ > BprimeCoatingMixtureTable_[0][0]) {
//      Info<< "p = " << p_
//          << " > " << BprimeCoatingMixtureTable_[0][0]
//          << nl;
      p_ =  BprimeCoatingMixtureTable_[0][0];
    }


    if (p_ < BprimeCoatingMixtureTable_[BprimeCoatingMixtureTable_.m() - 1][0]) {
//      Info<< "p = " << p_
//          << " < " << BprimeCoatingMixtureTable_[BprimeCoatingMixtureTable_.m() - 1][0]
//          << nl;
      p_ = BprimeCoatingMixtureTable_[BprimeCoatingMixtureTable_.m() - 1][0];
    }


//        Info<< BprimeCoatingMixtureTable_ << endl;

//        ::exit(1);
  }

  while (BprimeCoatingMixtureTable_[i_Pi][0] > p_) {
    i_Pi++;
  }

  indexPiMin_ = i_Pi;

  while (BprimeCoatingMixtureTable_[i_Pi][0] == BprimeCoatingMixtureTable_[indexPiMin_][0]) {
    i_Pi++;
    if (i_Pi == BprimeCoatingMixtureTable_.m())
      break;
  }

  indexPiMax_ = i_Pi - 1;

  if (p_ == BprimeCoatingMixtureTable_[indexPiMin_][0]) {
    indexPiPlus1Min_ = indexPiMin_;
    indexPiPlus1Max_ = indexPiMax_;
    f_p_ = 1;
  } else {
    indexPiPlus1Max_ = indexPiMin_ - 1;
    i_Pi = indexPiMin_ - 1;

    while (BprimeCoatingMixtureTable_[i_Pi][0] == BprimeCoatingMixtureTable_[indexPiPlus1Max_][0]) {
      i_Pi--;
      if (i_Pi == -1)
        break;
    }

    indexPiPlus1Min_ = i_Pi + 1;

    f_p_ =
        (::log(p_) - ::log(BprimeCoatingMixtureTable_[indexPiMin_][0])) /
        (
            ::log(BprimeCoatingMixtureTable_[indexPiPlus1Min_][0]) -
            ::log(BprimeCoatingMixtureTable_[indexPiMin_][0])
        );
  }
  // Info << "f_BprimeCoatingMixtureP = " << f_p_ << nl;
  return f_p_;
}

// indexPi: Returns the interpolation factor f_BpG_ and updates IndexBpGMin_, IndexBpGMax_, IndexBpGPlus1Min_, IndexBpGPlus1Max_
Foam::scalar Foam::BprimeCoatingMixtureTable::indexBprimeCoatingMixtureG(scalar BprimeCoatingMixtureG_, int i_min_, int i_max_)
{
  int i_BpGi(i_min_);

  /* To be removed later ... avoids problems for cases with no blowing or 'negative blowing', or excessive blowing */
  BprimeCoatingMixtureG_ = max(BprimeCoatingMixtureTable_[i_max_][1], BprimeCoatingMixtureG_); // CheckPoint
  /*   *   *   *   *   *   *   *   *   *   */

  if ((BprimeCoatingMixtureG_ > BprimeCoatingMixtureTable_[i_min_][1])) {
//    Info<< "B'g out of B' table range : " << BprimeCoatingMixtureG_
//        <<" > " << BprimeCoatingMixtureTable_[i_min_][1]
//        << nl;
//    Info<< "Set to max value available to force resolution "
//        "-- PREDICTIONS WILL BE INCORRECT -- "
//        "Update the B' table or use Mutation++"
//        << nl;

    BprimeCoatingMixtureG_ = min(BprimeCoatingMixtureTable_[i_min_][1], BprimeCoatingMixtureG_); // CheckPoint
  }

  while (BprimeCoatingMixtureTable_[i_BpGi][1] > BprimeCoatingMixtureG_) {
    i_BpGi++;
  }

  indexBpGiMin_ = i_BpGi;

  while (BprimeCoatingMixtureTable_[i_BpGi][1] == BprimeCoatingMixtureTable_[indexBpGiMin_][1]) {
    i_BpGi++;

    if (i_BpGi == i_max_ + 1)
      break;
  }

  indexBpGiMax_ = i_BpGi - 1;

  if (BprimeCoatingMixtureG_ == BprimeCoatingMixtureTable_[indexBpGiMin_][1]) {
    indexBpGiPlus1Min_ = indexBpGiMin_;
    indexBpGiPlus1Max_ = indexBpGiMax_;
    f_BpG_ = 1;
  } else {
    indexBpGiPlus1Max_ = indexBpGiMin_ - 1;
    i_BpGi = indexBpGiMin_ - 1;

    while (BprimeCoatingMixtureTable_[i_BpGi][1] == BprimeCoatingMixtureTable_[indexBpGiPlus1Max_][1]) {
      i_BpGi--;
      if (i_BpGi == i_min_ - 1)
        break;
    }

    indexBpGiPlus1Min_ = i_BpGi + 1;

    // linear interpolation
    f_BpG_=
        (BprimeCoatingMixtureG_ - BprimeCoatingMixtureTable_[indexBpGiMin_][1]) /
        (
            BprimeCoatingMixtureTable_[indexBpGiPlus1Min_][1]-
            BprimeCoatingMixtureTable_[indexBpGiMin_][1]
        );
    // log interpolation - only works when there are no zero BprimeCoatingMixtureG values in the B' table
    // f_BpG_= (::log(BprimeCoatingMixtureG_) - ::log(BprimeCoatingMixtureTable_[indexBpGiMin_][1])) / (::log(BprimeCoatingMixtureTable_[indexBpGiPlus1Min_][1])- ::log(BprimeCoatingMixtureTable_[indexBpGiMin_][1]));
  }

  return f_BpG_;
}

// BprimeCoatingMixtureC: Returns B'c (and updates hw) knowing p, B'g, T (interpolates B' tables)
Foam::scalar Foam::BprimeCoatingMixtureTable::BprimeCoatingMixtureC(scalar p_, scalar BprimeCoatingMixtureG_, scalar T_)
{
  // local variables
  scalar
  BpC_j_i_, BpC_jPlus1_i_, BpC_i_, BpC_iPlus1_,
            f_Ti_, hw_j_i_, hw_jPlus1_i_, hw_i_, hw_iPlus1_;

  // pressure indexes
  indexPi(p_);

  // BprimeCoatingMixtureG indexes for pressure i
  indexBprimeCoatingMixtureG(BprimeCoatingMixtureG_, indexPiMin_, indexPiMax_);

  // temperature indexes for sub-table i,j
  f_Ti_ = indexTi(T_, indexBpGiMin_, indexBpGiMax_);

  // B'c and hw for sub-table i,j
  BpC_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // temperature indexes for sub-table i,j+1
  f_Ti_ = indexTi(T_, indexBpGiPlus1Min_, indexBpGiPlus1Max_);

  // B'c for sub-table i,j+1
  BpC_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];
  // calculate B'c of table i
  BpC_i_ = f_BpG_ * BpC_jPlus1_i_ + (1 - f_BpG_) * BpC_j_i_;
  hw_i_  = f_BpG_ * hw_jPlus1_i_ + (1 - f_BpG_) * hw_j_i_;

  // BprimeCoatingMixtureG indexes for pressure i+1
  indexBprimeCoatingMixtureG(BprimeCoatingMixtureG_, indexPiPlus1Min_, indexPiPlus1Max_);

  // temperature indexes for sub-table i+1,j
  f_Ti_ = indexTi(T_, indexBpGiMin_, indexBpGiMax_);

  // B'c for sub-table i,j
  BpC_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // temperature indexes for sub-table i+1,j+1
  f_Ti_ = indexTi(T_, indexBpGiPlus1Min_, indexBpGiPlus1Max_);

  // B'c for sub-table i+1,j+1
  BpC_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // calculate B'c of table i+1
  BpC_iPlus1_ = f_BpG_ * BpC_jPlus1_i_ + (1 - f_BpG_) * BpC_j_i_;
  hw_iPlus1_  = f_BpG_ * hw_jPlus1_i_ + (1 - f_BpG_) * hw_j_i_;

  // B'c
  hw_  = f_p_ * hw_iPlus1_ + (1. - f_p_) * hw_i_;
  BpC_ = f_p_ * BpC_iPlus1_ + (1. - f_p_) * BpC_i_;

  return BpC_;
}

// BprimeCoatingMixtureC: Returns B'c (and updates hw) knowing p, B'g, T (interpolates B' tables using log)
Foam::scalar Foam::BprimeCoatingMixtureTable::BprimeCoatingMixtureClog(scalar p_, scalar BprimeCoatingMixtureG_, scalar T_)
{
  // local variables
  scalar BpC_j_i_, BpC_jPlus1_i_, BpC_i_, BpC_iPlus1_, f_Ti_, hw_j_i_, hw_jPlus1_i_, hw_i_, hw_iPlus1_;
  scalar petit = 1e-8;

  // pressure indexes
  indexPi(p_);

  // BprimeCoatingMixtureG indexes for pressure i
  indexBprimeCoatingMixtureG(BprimeCoatingMixtureG_, indexPiMin_, indexPiMax_);

  // temperature indexes for sub-table i,j
  f_Ti_ = indexTi(T_, indexBpGiMin_, indexBpGiMax_);

  // B'c and hw for sub-table i,j
  BpC_j_i_ = ::exp(f_Ti_ * ::log(petit + BprimeCoatingMixtureTable_[indexTiPlus1_][2]) + (1 - f_Ti_) * ::log(petit + BprimeCoatingMixtureTable_[indexTi_][2]));
  hw_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // temperature indexes for sub-table i,j+1
  f_Ti_ = indexTi(T_, indexBpGiPlus1Min_, indexBpGiPlus1Max_);

  // B'c for sub-table i,j+1
  BpC_jPlus1_i_ = ::exp(f_Ti_ * ::log(petit + BprimeCoatingMixtureTable_[indexTiPlus1_][2]) + (1 - f_Ti_) * ::log(petit + BprimeCoatingMixtureTable_[indexTi_][2]));
  hw_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // calculate B'c of table i
  BpC_i_ = ::exp(f_BpG_ * ::log(petit + BpC_jPlus1_i_) + (1 - f_BpG_) * ::log(petit + BpC_j_i_));
  hw_i_  = f_BpG_ * hw_jPlus1_i_ + (1 - f_BpG_) * hw_j_i_;

  // BprimeCoatingMixtureG indexes for pressure i+1
  indexBprimeCoatingMixtureG(BprimeCoatingMixtureG_, indexPiPlus1Min_, indexPiPlus1Max_);

  // temperature indexes for sub-table i+1,j
  f_Ti_ = indexTi(T_, indexBpGiMin_, indexBpGiMax_);

  // B'c for sub-table i,j
  BpC_j_i_ = ::exp(f_Ti_ * ::log(petit + BprimeCoatingMixtureTable_[indexTiPlus1_][2]) + (1 - f_Ti_) * ::log(petit + BprimeCoatingMixtureTable_[indexTi_][2]));
  hw_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // temperature indexes for sub-table i+1,j+1
  f_Ti_ = indexTi(T_, indexBpGiPlus1Min_, indexBpGiPlus1Max_);

  // B'c for sub-table i+1,j+1
  BpC_jPlus1_i_ = ::exp(f_Ti_ * ::log(petit + BprimeCoatingMixtureTable_[indexTiPlus1_][2]) + (1 - f_Ti_) * ::log(petit + BprimeCoatingMixtureTable_[indexTi_][2]));
  hw_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // calculate B'c of table i+1
  BpC_iPlus1_ = ::exp(f_BpG_ * ::log(petit + BpC_jPlus1_i_) + (1 - f_BpG_) * ::log(petit + BpC_j_i_));
  hw_iPlus1_  = f_BpG_ * hw_jPlus1_i_ + (1 - f_BpG_) * hw_j_i_;

  // B'c
  BpC_ = ::exp(f_p_ * ::log(petit + BpC_iPlus1_) + (1 - f_p_) * ::log(petit + BpC_i_));
  hw_ =   f_p_ * hw_iPlus1_ + (1 - f_p_) * hw_i_;

  return BpC_;
}

// hw: Returns hw (and updates B'c) knowing p, B'g, T (interpolates B' tables)
Foam::scalar Foam::BprimeCoatingMixtureTable::hw(scalar p_, scalar BprimeCoatingMixtureG_, scalar T_)
{
// local variables
  scalar BpC_j_i_, BpC_jPlus1_i_, BpC_i_, BpC_iPlus1_, f_Ti_, hw_j_i_, hw_jPlus1_i_, hw_i_, hw_iPlus1_;

  // pressure indexes
  indexPi(p_);

  // BprimeCoatingMixtureG indexes for pressure i
  indexBprimeCoatingMixtureG(BprimeCoatingMixtureG_, indexPiMin_, indexPiMax_);

  // temperature indexes for sub-table i,j
  f_Ti_ = indexTi(T_, indexBpGiMin_, indexBpGiMax_);

  // B'c and hw for sub-table i,j
  BpC_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // temperature indexes for sub-table i,j+1
  f_Ti_ = indexTi(T_, indexBpGiPlus1Min_, indexBpGiPlus1Max_);

  // B'c for sub-table i,j+1
  BpC_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // calculate B'c of table i
  BpC_i_ = f_BpG_ * BpC_jPlus1_i_ + (1 - f_BpG_) * BpC_j_i_;
  hw_i_  = f_BpG_ * hw_jPlus1_i_ + (1 - f_BpG_) * hw_j_i_;

  // BprimeCoatingMixtureG indexes for pressure i+1
  indexBprimeCoatingMixtureG(BprimeCoatingMixtureG_, indexPiPlus1Min_, indexPiPlus1Max_);

  // temperature indexes for sub-table i+1,j
  f_Ti_ = indexTi(T_, indexBpGiMin_, indexBpGiMax_);

  // B'c for sub-table i,j
  BpC_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_j_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // temperature indexes for sub-table i+1,j+1
  f_Ti_ = indexTi(T_, indexBpGiPlus1Min_, indexBpGiPlus1Max_);

  // B'c for sub-table i+1,j+1
  BpC_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][2] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][2];
  hw_jPlus1_i_ = f_Ti_ * BprimeCoatingMixtureTable_[indexTiPlus1_][4] + (1 - f_Ti_) * BprimeCoatingMixtureTable_[indexTi_][4];

  // calculate B'c of table i+1
  BpC_iPlus1_ = f_BpG_ * BpC_jPlus1_i_ + (1 - f_BpG_) * BpC_j_i_;
  hw_iPlus1_  = f_BpG_ * hw_jPlus1_i_ + (1 - f_BpG_) * hw_j_i_;

  // B'c
  BpC_ = f_p_ * BpC_iPlus1_ + (1 - f_p_) * (BpC_i_);
  hw_  = f_p_ * hw_iPlus1_ + (1 - f_p_) * hw_i_;

  return hw_;
}


// ************************************************************************* //

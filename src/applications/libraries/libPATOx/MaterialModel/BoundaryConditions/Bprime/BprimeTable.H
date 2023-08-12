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
    Foam::Bprime

Description
    B prime table class
    Extracts B'c

\*---------------------------------------------------------------------------*/

#ifndef BprimeTable_H
#define BprimeTable_H

#include "fvCFD.H"
#include "scalar.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Bprime class Declaration
\*---------------------------------------------------------------------------*/

class BprimeTable
{

// private

  RectangularMatrix<scalar> BprimeTable_; // BprimeTable itself
  int indexPiMin_;                        // First index of Pi : Pi <= P
  int indexPiMax_;                        // Last index of Pi : Pi <= P
  int indexPiPlus1Min_;                   // First index of P_{i+1} : P_{i+1} >= P
  int indexPiPlus1Max_;                   // Last index of P_{i+1} : P_{i+1} >= P
  scalar f_p_;                            // Interpolation coefficient for the pression
  int indexBpGiMin_;                      // First index of BpGi : BpGi <= BpG
  int indexBpGiMax_;                      // Last index of BpGi : BpGi <= BpG
  int indexBpGiPlus1Min_;                 // First index of BpG_{i+1} : BpG_{i+1} >= BpG
  int indexBpGiPlus1Max_;                 // Last index of BpG_{i+1} : BpG_{i+1} >= BpG
  scalar f_BpG_;                          // Interpolation coefficient for BpG
  scalar f_BpGi_;                         // Interpolation coefficient for BpGi
  scalar f_BpGiPlus1_;                    // Interpolation coefficient for BpGiPlus1
  int indexTi_;                           // Index of Ti : Ti <= T
  int indexTiPlus1_;                      // Index of T_{i+1} : T_{i+1} >= T
  scalar f_T_;                            // Interpolation coefficient for T
  scalar BpC_;                            // B'c
  scalar hw_;                             // h_w
  scalar small_;                          // small number


 public:

  // Constructors
  BprimeTable();

  //- Construct from components
  BprimeTable(RectangularMatrix<scalar>& BprimeTable);


  // Destructor
  virtual ~BprimeTable();

  // Access Functions
  int indexPiMin() const {
    return indexPiMin_;
  }
  int indexPiMax() const {
    return indexPiMax_;
  }
  int indexPiPlus1Min() const {
    return indexPiPlus1Min_;
  }
  int indexPiPlus1Max() const {
    return indexPiPlus1Max_;
  }
  int indexBpGiMin() const {
    return indexBpGiMin_;
  }
  int indexBpGiMax() const {
    return indexBpGiMax_;
  }
  int indexBpGiPlus1Min() const {
    return indexBpGiPlus1Min_;
  }
  int indexBpGiPlus1Max() const {
    return indexBpGiPlus1Max_;
  }
  int indexTi() const {
    return indexTi_;
  }
  int indexTiPlus1() const {
    return indexTiPlus1_;
  }
  scalar Ti() const {
    return BprimeTable_[indexTi_][3];
  }
  scalar TiPlus1() const {
    return BprimeTable_[indexTiPlus1_][3];
  }
  scalar BpC() const {
    return BpC_;
  }
  scalar hw() const {
    return hw_;
  }

  // Member Functions
  // returns the Temperature interpolation factor and updates the indexes of Ti and Ti+1
  scalar indexTi(scalar T_, int i_min_, int i_max_);
  // returns the pressure interpolation factor and updates the indexes of the pressure subtables
  scalar indexPi(scalar p_);
  // returns the BprimeG interpolation factor and updates the indexes of the BprimeG subtables contained in a given pressure subtable
  scalar indexBprimeG(scalar BprimeG_, int i_min_, int i_max_);
  // returns BprimeC with 'semi-log' interpolation
  scalar BprimeC(scalar p_, scalar BprimeG_, scalar T_);
  // returns BprimeC with log interpolation
  scalar BprimeClog(scalar p_, scalar BprimeG_, scalar T_);
  // returns hw
  scalar hw(scalar p_, scalar BprimeG_, scalar T_);

}; // End class Bprime
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

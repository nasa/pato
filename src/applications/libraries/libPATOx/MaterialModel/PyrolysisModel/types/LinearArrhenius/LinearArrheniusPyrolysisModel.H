/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    along with OpenFOAM.  If LinearArrheniust, see <http://www.gnu.org/licenses/>.

Class
    Foam::LinearArrheniusPyrolysisModel

Description

SourceFiles
    LinearArrheniusPyrolysisModel.H
    LinearArrheniusPyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef LinearArrheniusPyrolysisModel_H
#define LinearArrheniusPyrolysisModel_H

#include "simplePyrolysisModel.H"
#include "simpleMaterialChemistryModel.H"
#include "simpleEnergyModel.H"
#include "PList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class LinearArrheniusPyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class LinearArrheniusPyrolysisModel
:
public simplePyrolysisModel
{

 protected:

  //- Number of solid phases
  scalar& nSolidPhases_;

  //- Number of pyrolysis reactions per phase
  PList<scalar> nPyroReac;

  //- Advancement of the pyrolysis reactions: 1 to 0
  volScalarField& tau;

  //- Total pyrolysis gas production [kg/m3/s]
  volScalarField& piTotal_;

  //- Solid virgin mass density [kg/m3]
  volScalarField& rho_v_;

  //- Solid charred mass density [kg/m3]
  volScalarField& rho_c_;

  //- Sum of the norm of factor F of the pyrolysis reactions
  volScalarField& normFp_sum;

  //- Species/elements pyrolysis gas production [kg/m3/s]
  PList<volScalarField> pi_;

  //- Pyrolysis gas production per reaction [kg/m3/s]
  PList<volScalarField> piPyroReac_;

  //- Material Chemistry Model
  simpleMaterialChemistryModel& materialChemistryModel;

  //- Species names
  const wordList& speciesNames_;

  //- Elements names
  const wordList& elementNames_;

  //- Energy Model
  simpleEnergyModel& energyModel;

  //- Temperature [K]
  const volScalarField& Ta_;

  //- Solid mass density [kg/m3]
  volScalarField& rho_s_;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar& R_;

  //- Number total of the pyrolysis reactions
  int sizeXsi;

  //- Hash table: column 0 = phase number, column 1 = phase reaction number
  RectangularMatrix<label> hashXsi;

  //- List of all the pyrolysis reactions for every phases
  PList<volScalarField> Xsi;

  //- Factor F of the pyrolysis reactions
  PList<dimensionedScalar> Fp ;

  //- Norm of factor F of the pyrolysis reactions
  PList<volScalarField> normFp ;

  //- Factor A of the pyrolysis reactions
  PList<dimensionedScalar> Ap;

  //- Factor E of the pyrolysis reactions
  PList<dimensionedScalar> Ep;

  //- Factor m of the pyrolysis reactions
  PList<dimensionedScalar> mp;

  //- Factor n of the pyrolysis reactions
  PList<dimensionedScalar> np;

  //- Factor T of the pyrolysis reactions
  PList<dimensionedScalar> Tp;

  //- Factor h of the pyrolysis reactions
  PList<dimensionedScalar>& hp;

  //- VolScalarField of factor A of the pyrolysis reactions
  PList<volScalarField> AField;

  //- Zeta (element) of the pyrolysis reactions
  PList<PList<scalar>> Zp;

  //- Gamma (species) of the pyrolysis reactions
  PList<PList<scalar>> Yp;

  //- Material chemistry type
  word MaterialChemistryType_;

  //- Solid initial volume fraction [-]
  PList<volScalarField> solidEpsI_;

  //- Solid initial mass density [kg/m3]
  PList<volScalarField> solidRhoI_;

  //- Solid volume fraction [-]
  PList<volScalarField> solidEps_;

  //- Solid mass density [kg/m3]
  PList<volScalarField> solidRho_;

  //- Unity volScalarField
  volScalarField& oneVolScalarField;

  //- Minimum value of the pyrolysis gas production
  scalar minPiPyro_;

 public:

  //- Runtime type information
  TypeName("LinearArrhenius");

  // Constructors

  //- Construct from mesh and phase name
  LinearArrheniusPyrolysisModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~LinearArrheniusPyrolysisModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Initialize
  virtual void initialize();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

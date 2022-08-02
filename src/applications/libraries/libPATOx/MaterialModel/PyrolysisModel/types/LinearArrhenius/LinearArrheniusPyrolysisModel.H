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
#include "simpleMaterialPropertiesModel.H"
#include "simpleMaterialChemistryModel.H"

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

  //- Mesh reference
  const fvMesh& mesh_;

  //- Dictionary name
  const word& dictName_;

  //- Material Properties Model
  simpleMaterialPropertiesModel& materialPropertiesModel;

  //- Material Chemistry Model
  simpleMaterialChemistryModel& materialChemistryModel;

  //- Number of solid phases
  scalar& nSolidPhases_;

  //- Number total of the pyrolysis reactions
  scalar sizeXsi;

  //- Number of pyrolysis reactions per phase
  PtrList<scalar> nPyroReac;

  //- Hash table: column 0 = phase number, column 1 = phase reaction number
  RectangularMatrix<label> hashXsi;

  //- List of all the pyrolysis reactions for every phases
  PtrList<volScalarField> Xsi;

  //- Factor F of the pyrolysis reactions
  PtrList<dimensionedScalar> Fp ;

  //- Norm of factor F of the pyrolysis reactions
  PtrList<volScalarField> normFp ;

  //- Factor A of the pyrolysis reactions
  PtrList<dimensionedScalar> Ap;

  //- Factor E of the pyrolysis reactions
  PtrList<dimensionedScalar> Ep;

  //- Factor m of the pyrolysis reactions
  PtrList<dimensionedScalar> mp;

  //- Factor n of the pyrolysis reactions
  PtrList<dimensionedScalar> np;

  //- Factor T of the pyrolysis reactions
  PtrList<dimensionedScalar> Tp;

  //- Factor h of the pyrolysis reactions
  PtrList<dimensionedScalar>& hp;

  //- VolScalarField of factor A of the pyrolysis reactions
  PtrList<volScalarField> AField;

  //- Zeta (element) of the pyrolysis reactions
  RectangularMatrix<scalar> Zp;

  //- Gamma (species) of the pyrolysis reactions
  RectangularMatrix<scalar> Yp;

  //- Advancement of the pyrolysis reactions: 1 to 0
  volScalarField& tau;

  //- Solid initial volume fraction [-]
  PtrList<volScalarField>& solidEpsI;

  //- Solid initial mass density [kg/m3]
  PtrList<volScalarField>& solidRhoI;

  //- Temperature [K]
  volScalarField& Ta_;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar& R_;

  //- Total pyrolysis gas production [kg/m3/s]
  volScalarField& piTotal_;

  //- Species/elements pyrolysis gas production [kg/m3/s]
  PtrList<volScalarField>& pi_;

  //- Pyrolysis gas production per reaction [kg/m3/s]
  PtrList<volScalarField>& piPyroReac_;

  //- Material properties directory
  const fileName materialPropertiesDirectory;

  //- Material properties dictionary
  IOdictionary constantPropertiesDictionary;

  //- Material chemistry type
  word MaterialChemistryType_;

  //- Solid mass density [kg/m3]
  PtrList<volScalarField>& solidRho_;

  //- Solid volume fraction [-]
  PtrList<volScalarField>& solidEps_;

  //- Solid mass density [kg/m3]
  volScalarField& rho_s_;

  //- Solid virgin mass density [kg/m3]
  volScalarField& rho_v_;

  //- Solid charred mass density [kg/m3]
  volScalarField& rho_c_;

  //- Unity volScalarField
  volScalarField oneVolScalarField;

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

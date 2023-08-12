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

#include "basicUserThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(basicUserThermo, 0);
defineRunTimeSelectionTable(basicUserThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicUserThermo::basicUserThermo
(
    const fvMesh& mesh,
    const word& dictName
):
mesh_(mesh),
currentP_(1e4),
currentE_(1e5),
currentT_(1e3),
currentRho_(1.),
currentPsi_(1.),
currentMu_(1.),
currentNu_(1.),
currentSigma_(1.),
currentGamma_(1.),
currentKappa_(1.),
currentAlpha_(1.),
p_(meshLookupOrConstructScalar(mesh, "p")),
T_(meshLookupOrConstructScalar(mesh, "T")),
e_
(
    IOobject
    (
        "e",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", pow(dimLength, 2)/pow(dimTime, 2), scalar(0.0)),
    eBoundaryTypes()
),
rho_(meshLookupOrConstructScalar(mesh, "rho", dimensionedScalar("0", dimMass/pow3(dimLength), scalar(0.0)))),
psi_(meshLookupOrConstructScalar(mesh, "psi", dimensionedScalar("0", pow(dimTime, 2)/pow(dimLength, 2), scalar(0.0)))),
sigma_(meshLookupOrConstructScalar(mesh, "sigma", dimensionedScalar("0", dimensionSet(1, -3, 3, 0, 0, 2, 0), scalar(0.0)))),
gamma_(meshLookupOrConstructScalar(mesh, "gamma", dimensionedScalar("0", dimless, scalar(0.0)))),
alpha_(meshLookupOrConstructScalar(mesh, "alpha", dimensionedScalar("0", dimensionSet(1, -1, -1, 0, 0), scalar(0.0)))),
mu_(meshLookupOrConstructScalar(mesh, "mu", dimensionedScalar("0", dimMass/dimLength/dimTime, scalar(0.0)))),
nu_(meshLookupOrConstructScalar(mesh, "nu", dimensionedScalar("0", pow(dimMass,2)/dimTime, scalar(0.0)))),
kappa_(meshLookupOrConstructScalar(mesh, "kappa", dimensionedScalar("0", dimMass*dimLength/(pow3(dimTime)*dimTemperature), scalar(0.0)))),
kappaEff_(meshLookupOrConstructScalar(mesh, "kappaEff", dimensionedScalar("0", dimMass*dimLength/(pow3(dimTime)*dimTemperature), scalar(0.0)))),
alphaEff_(meshLookupOrConstructScalar(mesh, "alphaEff", dimensionedScalar("0", dimMass/dimLength/dimTime, scalar(0.0)))),
alphahe_(meshLookupOrConstructScalar(mesh, "alphahe", dimensionedScalar("0", dimMass/dimLength/dimTime, scalar(0.0)))),
isInConstructor_(1)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::basicUserThermo> Foam::basicUserThermo::New
(
    const fvMesh& mesh,
    const word& dictName
)
{
  IOdictionary thermoDict_
  (
      IOobject
      (
          dictName,
          mesh.time().constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
      )
  );

  word typeName(thermoDict_.lookup("type"));
  Info << " --- Mutation++:  Thermo model selecting type \"" << typeName << "\"" << endl << endl;

  typename basicUserThermo::fvMeshConstructorTable::iterator cstrIter =
      basicUserThermo::fvMeshConstructorTablePtr_->find(typeName);

  if (cstrIter == basicUserThermo::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  basicUserThermo::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << basicUserThermo::typeName << " type "
        << typeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<basicUserThermo>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicUserThermo::~basicUserThermo()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

//- void to configure the mixture pointer for Mutation++
void Foam::basicUserThermo::configureMixture()
{
  Info << " --- Mutation++:  ===  Configuring mixture --> Started " << endl;

  Mutation::MixtureOptions opts(mixtureName_);
  opts.setStateModel(stateModelName_);

  this->pMutationMix_.reset(new Mutation::Mixture(opts));

  this->nSpecies_ = pMutationMix_->nSpecies();
  this->nElements_ = pMutationMix_->nElements();

  this->speciesList_.resize(nSpecies_);
  this->speciesMolarMassList_.resize(nSpecies_);
  this->elementsList_.resize(nElements_);

  this->currentYs_.resize(nSpecies_);
  this->currentYe_.resize(nElements_);
  this->currentXs_.resize(nSpecies_);
  this->currentXe_.resize(nElements_);
  this->currentRhos_.resize(nSpecies_);
  this->currentRhoe_.resize(nElements_);

  this->Ys_.resize(nSpecies_);
  this->Ye_.resize(nElements_);
  this->Xs_.resize(nSpecies_);
  this->Xe_.resize(nElements_);
  this->rhos_.resize(nSpecies_);
  this->rhoe_.resize(nElements_);

  forAll(speciesList_,is) {
    this->speciesList_[is] = this->pMutationMix_->speciesName(is);
    this->speciesMolarMassList_[is] = this->pMutationMix_->speciesMw(is);
  }

  forAll(elementsList_,ie) {
    this->elementsList_[ie] = this->pMutationMix_->elementName(ie);
  }

  this->nMassEqns_ = this->pMutationMix_->nMassEqns();
  this->nEnergyEqns_ = this->pMutationMix_->nEnergyEqns();

  Info << " --- Mutation++:  Mixture name :               " << this->mixtureName_ << endl;
  Info << " --- Mutation++:  State model name :           " << this->stateModelName_ << endl;
  Info << " --- Mutation++:  Number of elements :         " << this->nElements_ << endl;
  Info << " --- Mutation++:  Number of species  :         " << this->nSpecies_ << endl;
  Info << " --- Mutation++:  List of elements :           ";
  forAll(elementsList_, ie) {
    Info << this->elementsList_[ie];
    if (ie == this->nElements_-1) {
      Info << endl;
    } else {
      Info << " ,  ";
    };
  };
  Info << " --- Mutation++:  List of species  :           ";
  forAll(speciesList_,  is) {
    Info << this->speciesList_[is];
    if (is == this->nSpecies_ -1) {
      Info << endl;
    } else {
      Info << " ,  ";
    };
  };

  Info << " --- Mutation++:  ===  Configuring mixture --> Finished " << endl << endl;
}

// ====================================================== //

void Foam::basicUserThermo::createCompositionFields(const fvMesh& mesh)
{
  for(int ie = 0; ie < this->nElements(); ie++) {
    word nameXe = "Xe[" + this->getElementName(ie) +"]";
    word nameYe = "Ye[" + this->getElementName(ie) +"]";
    word nameRhoe = "rhoe[" + this->getElementName(ie) +"]";
    dimensionedScalar dimScalar_("0",dimless,scalar(0.0));

    Xe_.set(ie, new volScalarField
            (
                IOobject
                (
                    nameXe,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimScalar_
            )
           );

    Ye_.set(ie, new volScalarField
            (
                IOobject
                (
                    nameYe,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimScalar_
            )
           );

    rhoe_.set(ie, new volScalarField
              (
                  IOobject
                  (
                      nameRhoe,
                      mesh.time().timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE
                  ),
                  mesh,
                  dimScalar_
              )
             );
  }

  for(int is = 0; is < this->nSpecies(); is++) {
    word nameXs = "Xs[" + this->getSpeciesName(is) +"]";
    word nameYs = "Ys[" + this->getSpeciesName(is) +"]";
    word nameRhos = "rhos[" + this->getSpeciesName(is) +"]";
    dimensionedScalar dimScalar_("0",dimless,scalar(0.0));

    Xs_.set(is, new volScalarField
            (
                IOobject
                (
                    nameXs,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimScalar_
            )
           );

    Ys_.set(is, new volScalarField
            (
                IOobject
                (
                    nameYs,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimScalar_
            )
           );

    rhos_.set(is, new volScalarField
              (
                  IOobject
                  (
                      nameRhos,
                      mesh.time().timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE
                  ),
                  mesh,
                  dimScalar_
              )
             );
  }
}

// ====================================================== //

void Foam::basicUserThermo::createElementFluxFields(const fvMesh& mesh)
{
  dimensionedVector dimVector1_("vector1", dimMass/dimTime/pow(dimLength,2), Foam::vector(0.,0.,0.));

  diffMassFluxes_.resize(this->nElements());
  elementsMassDiffFluxesGradP_.resize(this->nElements());
  elementsMassDiffFluxesGradT_.resize(this->nElements());
  elementsMassDiffFluxesGradX_.resize(this->nElements()*this->nElements());
  elementsEnergyDiffFluxesGradX_.resize(this->nElements());

  for(int ie = 0; ie < this->nElements(); ie++) {
    word name = "diffMassFlux[" + this->getElementName(ie) +"]";
    diffMassFluxes_.set(ie, new volVectorField
                        (
                            IOobject
                            (
                                name,
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh,
                            dimVector1_
                        )
                       );
  }

  dimensionedScalar dimScalar1_("0", dimTime, scalar(0.0));
  dimensionedScalar dimScalar2_("0", dimMass/(dimTime*dimLength*dimTemperature), scalar(0.0));
  dimensionedScalar dimScalar3_("0", dimMass/(dimMass*dimTime), scalar(0.0));
  dimensionedScalar dimScalar4_("0", pow(dimLength,2)/dimTime, scalar(0.0));
  dimensionedScalar dimScalar5_("0", dimMass*dimLength/(pow3(dimTime)*dimTemperature), scalar(0.0));
  dimensionedScalar dimScalar6_("0", dimMass*dimLength/pow3(dimTime), scalar(0.0));

  for(int ie = 0; ie < this->nElements(); ie++) {
    word nameFP = "massFluxP["   + this->getElementName(ie) + "]";
    word nameFT = "massFluxT["   + this->getElementName(ie) + "]";
    word nameEX = "energyFluxX[" + this->getElementName(ie) + "]";

    elementsMassDiffFluxesGradP_.set(ie, new volScalarField
                                     (
                                         IOobject
                                         (
                                             nameFP,
                                             mesh.time().timeName(),
                                             mesh,
                                             IOobject::NO_READ,
                                             IOobject::NO_WRITE
                                         ),
                                         mesh,
                                         dimScalar1_
                                     )
                                    );

    elementsMassDiffFluxesGradT_.set(ie, new volScalarField
                                     (
                                         IOobject
                                         (
                                             nameFT,
                                             mesh.time().timeName(),
                                             mesh,
                                             IOobject::NO_READ,
                                             IOobject::NO_WRITE
                                         ),
                                         mesh,
                                         dimScalar2_
                                     )
                                    );

    elementsEnergyDiffFluxesGradX_.set(ie, new volScalarField
                                       (
                                           IOobject
                                           (
                                               nameEX,
                                               mesh.time().timeName(),
                                               mesh,
                                               IOobject::NO_READ,
                                               IOobject::NO_WRITE
                                           ),
                                           mesh,
                                           dimScalar6_
                                       )
                                      );

    for(int ie2 = 0; ie2 < this->nElements(); ie2++) {
      word nameFX = "massFluxX[" + this->getElementName(ie) + "][" + this->getElementName(ie2) + "]";

      int index = ie * this->nElements() + ie2;

      elementsMassDiffFluxesGradX_.set(index, new volScalarField
                                       (
                                           IOobject
                                           (
                                               nameFX,
                                               mesh.time().timeName(),
                                               mesh,
                                               IOobject::NO_READ,
                                               IOobject::NO_WRITE
                                           ),
                                           mesh,
                                           dimScalar3_
                                       )
                                      );
    }
  }

  elementsEnergyDiffFluxesGradP_.set(new volScalarField
                                     (
                                         IOobject
                                         (
                                             "energyFluxP",
                                             mesh.time().timeName(),
                                             mesh,
                                             IOobject::NO_READ,
                                             IOobject::NO_WRITE
                                         ),
                                         mesh,
                                         dimScalar4_
                                     )
                                    );

  elementsEnergyDiffFluxesGradT_.set(new volScalarField
                                     (
                                         IOobject
                                         (
                                             "energyFluxT",
                                             mesh.time().timeName(),
                                             mesh,
                                             IOobject::NO_READ,
                                             IOobject::NO_WRITE
                                         ),
                                         mesh,
                                         dimScalar5_
                                     )
                                    );

  qDiff_.set(new volVectorField
             (
                 IOobject
                 (
                     "qDiff",
                     mesh.time().timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE
                 ),
                 mesh,
                 dimensionedVector("vector0", dimMass/pow3(dimTime), Foam::vector(0.,0.,0.))
             )
            );
}

// ====================================================== //

void Foam::basicUserThermo::createFiniteRateFields(const fvMesh& mesh)
{
  dimensionedVector dimVector1_("vector1", dimLength/dimTime, Foam::vector(0.,0.,0.));
  dimensionedScalar dimScalar_("0",dimMass/dimTime/pow3(dimLength),scalar(0.));

  for(int is = 0; is < this->nSpecies(); is++) {
    word name = "Vdiff[" + this->getSpeciesName(is) +"]";

    diffVs_.set(is, new volVectorField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimVector1_
                )
               );

    word name1 = "Omega[" + this->getSpeciesName(is) +"]";

    omegas_.set(is, new volScalarField
                (
                    IOobject
                    (
                        name1,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimScalar_
                )
               );
  }

  qDiff_.set(new volVectorField
             (
                 IOobject
                 (
                     "qDiff",
                     mesh.time().timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE
                 ),
                 mesh,
                 dimensionedVector("vector0", dimMass/pow3(dimTime), Foam::vector(0.,0.,0.))
             )
            );
}

// ====================================================== //

const Foam::word& Foam::basicUserThermo::getMixtureName() const
{
  return mixtureName_;
}

// ====================================================== //

const Foam::word& Foam::basicUserThermo::getStateModelName() const
{
  return stateModelName_;
}

// ====================================================== //

const int& Foam::basicUserThermo::nSpecies() const
{
  return nSpecies_;
}

// ====================================================== //

const int& Foam::basicUserThermo::nElements() const
{
  return nElements_;
}

// ====================================================== //

const int& Foam::basicUserThermo::nMassEqns() const
{
  return nMassEqns_;
}

// ====================================================== //

const int& Foam::basicUserThermo::nEnergyEqns() const
{
  return nEnergyEqns_;
}
// ====================================================== //

const Foam::wordList& Foam::basicUserThermo::getElementsList() const
{
  return elementsList_;
}

// ====================================================== //

Foam::word& Foam::basicUserThermo::getElementName(const int ie)
{
  return elementsList_[ie];
}

// ====================================================== //

const Foam::wordList& Foam::basicUserThermo::getSpeciesList() const
{
  return speciesList_;
}

// ====================================================== //

Foam::word& Foam::basicUserThermo::getSpeciesName(const int is)
{
  return speciesList_[is];
}

// ====================================================== //

const Foam::scalarList& Foam::basicUserThermo::speciesMolarMassList() const
{
  return speciesMolarMassList_;
}

// ====================================================== //

Foam::scalar& Foam::basicUserThermo::speciesMolarMass(const int is)
{
  return speciesMolarMassList_[is];
}

// ====================================================== //

Foam::scalar& Foam::basicUserThermo::RU()
{
  return RU_;
}

// ====================================================== //

Foam::volScalarField& Foam::basicUserThermo::p()
{
  return p_;
}

// ====================================================== //

Foam::volScalarField& Foam::basicUserThermo::e()
{
  return e_;
}

// ====================================================== //

Foam::volScalarField& Foam::basicUserThermo::T()
{
  return T_;
}

// ====================================================== //

const Foam::volScalarField& Foam::basicUserThermo::rho() const
{
  return rho_;
}

// ====================================================== //

Foam::volScalarField& Foam::basicUserThermo::rho()
{
  return rho_;
}

// ====================================================== //

const Foam::volScalarField& Foam::basicUserThermo::psi() const
{
  return psi_;
}

// ====================================================== //

Foam::volScalarField& Foam::basicUserThermo::sigma()
{
  return sigma_;
}

// ====================================================== //

const Foam::volScalarField& Foam::basicUserThermo::gamma() const
{
  return gamma_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsMassFractions()
{
  return Ye_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsMoleFractions()
{
  return Xe_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsDensities()
{
  return rhoe_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::speciesMassFractions()
{
  return Ys_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::speciesMoleFractions()
{
  return Xs_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::speciesDensities()
{
  return rhos_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsMassDiffFluxesGradP()
{
  return elementsMassDiffFluxesGradP_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsMassDiffFluxesGradT()
{
  return elementsMassDiffFluxesGradT_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsMassDiffFluxesGradX()
{
  return elementsMassDiffFluxesGradX_;
}

// ====================================================== //

Foam::autoPtr<Foam::volScalarField>& Foam::basicUserThermo::elementsEnergyDiffFluxesGradP()
{
  return elementsEnergyDiffFluxesGradP_;
}

// ====================================================== //

Foam::autoPtr<Foam::volScalarField>& Foam::basicUserThermo::elementsEnergyDiffFluxesGradT()
{
  return elementsEnergyDiffFluxesGradT_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::elementsEnergyDiffFluxesGradX()
{
  return elementsEnergyDiffFluxesGradX_;
}

// ====================================================== //

Foam::PtrList<Foam::volVectorField>& Foam::basicUserThermo::elementsDiffusionMassFluxes()
{
  return diffMassFluxes_;
}

// ====================================================== //

Foam::PtrList<Foam::volVectorField>& Foam::basicUserThermo::speciesDiffusionVelocities()
{
  return diffVs_;
}

// ====================================================== //

Foam::PtrList<Foam::volScalarField>& Foam::basicUserThermo::speciesNetMassProductionRates()
{
  return omegas_;
}

// ====================================================== //

Foam::autoPtr<Foam::volVectorField>& Foam::basicUserThermo::diffHeatFlux()
{
  return qDiff_;
}


// ====================================================== //

const Foam::volScalarField& Foam::basicUserThermo::alpha() const
{
  return alpha_;
}

// ====================================================== //

const Foam::scalarField& Foam::basicUserThermo::alpha
(
    const label patchi
) const
{
  return alpha_.boundaryField()[patchi];
}

// ====================================================== //

Foam::tmp<Foam::volScalarField>  Foam::basicUserThermo::mu() const
{
  return mu_;
}

// ====================================================== //


Foam::tmp<Foam::scalarField>  Foam::basicUserThermo::mu(const label patchi) const
{
  return mu_.boundaryField()[patchi];
}

// ====================================================== //


Foam::tmp<Foam::volScalarField>  Foam::basicUserThermo::nu() const
{
  return nu_;
}

// ====================================================== //


Foam::tmp<Foam::scalarField>  Foam::basicUserThermo::nu(const label patchi) const
{
  return nu_.boundaryField()[patchi];
}

// ====================================================== //


Foam::tmp<Foam::volScalarField>  Foam::basicUserThermo::kappa() const
{
  return kappa_;
}

// ====================================================== //

Foam::tmp<Foam::scalarField>  Foam::basicUserThermo::kappa
(
    const label patchi
) const
{
  return kappa_.boundaryField()[patchi];
}

// ====================================================== //

Foam::tmp<Foam::volScalarField>  Foam::basicUserThermo::alphahe() const
{
  return alphahe_;
}

// ====================================================== //

Foam::tmp<Foam::scalarField>  Foam::basicUserThermo::alphahe
(
    const label patchi
) const
{
  return alphahe_.boundaryField()[patchi];
}

// ====================================================== //

Foam::tmp<Foam::volScalarField>  Foam::basicUserThermo::kappaEff
(
    const volScalarField& kappat
) const
{
  // kappaEff = k + kt
  return kappaEff_;
}

// ====================================================== //

Foam::tmp<Foam::scalarField>  Foam::basicUserThermo::kappaEff
(
    const scalarField& kappat,
    const label patchi
) const
{
  return kappaEff_.boundaryField()[patchi];
}

// ====================================================== //

Foam::tmp<Foam::volScalarField>  Foam::basicUserThermo::alphaEff
(
    const volScalarField& alphat
) const
{
  // alphaEff = alpha + alphat with alpha = k/Cv
  return alphaEff_;
}

// ====================================================== //

Foam::tmp<Foam::scalarField>  Foam::basicUserThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
  return alphaEff_.boundaryField()[patchi];
}

// ************************************************************************* //


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

#include "pureThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pureThermo::pureThermo
(
    Time& runTime
)
  :
basicFluidModel(runTime)
{

#include "createMesh.H" // Create mesh object
  const word dictName = "thermoDict";
  IOdictionary thermoDict
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
  const word thermoType(thermoDict.lookup("type"));

  Info << "- pureThermo: creating pointer to thermo object " << endl;
  autoPtr<basicUserThermo> basicModel_(basicUserThermo::New(mesh,dictName));

  Info << "- pureThermo: creating thermo object " << endl;
  basicUserThermo& thermo = basicModel_();

//- Create the user inputs list (p,T,Y)

  scalarList pressureList(thermoDict.lookup("pressureList"));
  scalarList temperatureList(thermoDict.lookup("temperatureList"));
  scalarList rhoList(thermoDict.lookupOrDefault<scalarList>("rhoList",scalarList(0)));
  scalarList energyList(thermoDict.lookupOrDefault<scalarList>("energyList",scalarList(0)));
  List<scalarList> * speciesMassFractionList(0);
  List<scalarList> * elementsMassFractionList(0);

  if (thermoType=="FiniteRateMutation") {
    speciesMassFractionList = new List<scalarList>(thermoDict.lookup("speciesMassFractionList"));
    List<scalarList>& list = *speciesMassFractionList;
    const wordList speciesName(basicModel_->getSpeciesList());

    forAll(list, i) {
      if (list[i].size()!= basicModel_->nSpecies()) {
        FatalError << "Number of species in speciesMassFractionList(" << list[i].size() << ") not coherent with the number of species in \"" << (word) thermoDict.lookup("mutationMixture") << "\" mutationMixture(" << basicModel_->nSpecies() << ")" << nl;
        FatalError << "\"" << (word) thermoDict.lookup("mutationMixture") << "\" mutationMixture species = " << speciesName << exit(FatalError);
      }

      Info << "- pureThermo: speciesMassFractionList = " << endl;

      forAll(list, i) {
        Info << "\t\t#" << i+1 << ": ";

        forAll(list[i], j) {
          if (list[i][j]!=0) {
            Info << "Y[" << speciesName[j] << "]=" << list[i][j] << " ";
          }
        }

        Info << endl;
      }
    }
  }

  if ((thermoType=="PTXEquilMutation") || (thermoType=="PTXEquilTabulated")) {
    elementsMassFractionList = new List<scalarList>(thermoDict.lookup("elementsMassFractionList"));
    List<scalarList>& list = *elementsMassFractionList;
    const wordList elementsName(basicModel_->getElementsList());

    forAll(list, i) {
      if (list[i].size()!= basicModel_->nElements()) {
        FatalError << "Number of elements in elementsMassFractionList(" << list[i].size() << ") not coherent with the number of elements in \"" << (word) thermoDict.lookup("mutationMixture") << "\" mutationMixture(" << basicModel_->nElements() << ")" << nl;
        FatalError << "\"" << (word) thermoDict.lookup("mutationMixture") << "\" mutationMixture elements = " << elementsName << exit(FatalError);
      }

      Info << "- pureThermo: elementsMassFractionList = " << endl;

      forAll(list, i) {
        Info << "\t\t#" << i+1 << ": ";

        forAll(list[i], j) {
          if (list[i][j]!=0) {
            Info << "Y[" << elementsName[j] << "]=" << list[i][j] << " ";
          }
        }

        Info << endl;
      }
    }
  }



  Info << "- pureThermo: creating fields " << endl;
// Mixture density [kg/m^3]
  volScalarField rho
  (
      IOobject
      (
          "rho",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("0",dimMass/pow3(dimLength),scalar(1.0))
  );

  // Mass averaged bulk velocity of the gas [m/s]
  volVectorField U
  (
      IOobject
      (
          "U",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("0",dimLength/dimTime,vector(0,0,0))
  );

  // rho*U [kg/s/m^2]
  volVectorField rhoU
  (
      IOobject
      (
          "rhoU",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      rho*U
  );

  // Mass flux [kg/s]
  surfaceScalarField phi("phi", fvc::flux(rhoU));


  autoPtr<compressibleUser::turbulenceModel> turbulence
  (
      compressibleUser::turbulenceModel::New
      (
          rho,
          U,
          phi,
          thermo
      )
  );



  Info << "- pureThermo: writing output" << endl;

  if(!isDir("output/"+thermoType)) {
    system("mkdir -p output/"+thermoType);
  }



//- Create the list of output files
  PtrList<PtrList<OFstream> > osList;

  if (thermoType == "FiniteRateMutation") {
    osList.resize(speciesMassFractionList->size());
    List<scalarList>& list = *speciesMassFractionList;

    forAll(osList, sI) {
      osList.set(sI, new PtrList<OFstream>(pressureList.size()));

      forAll(osList[sI], pI) {
        word fileName = "output/" + thermoType + "/" +thermoType + "_mix" + name(sI) + "_P" + name(pressureList[pI]) + ".out";
        osList[sI].set(pI,  new OFstream(fileName));
        const wordList speciesName(basicModel_->getSpeciesList());
        osList[sI][pI] << "//T(K) p(Pa) rho(kg/m3) psi(s2/m2) e(J/kg) alpha(kg/m/s) alphaEff(kg/m/s) sigma(S/m) gamma(-) mu(Pa.s) muEff(Pa.s) ";
        for(int ie = 0; ie < basicModel_->nElements(); ie++) {
          osList[sI][pI] 	<< "Y[" << basicModel_->getElementName(ie) << "](-) ";
        }
        for(int is = 0; is < basicModel_->nSpecies(); is++) {
          osList[sI][pI] 	<< "Y[" << basicModel_->getSpeciesName(is) << "](-) ";
        }
        if (thermoType == "FiniteRateMutation") {
          for(int is = 0; is < basicModel_->nSpecies(); is++) {
            osList[sI][pI] 	<< "omega[" << basicModel_->getSpeciesName(is) << "](kg/m3/s) ";
          }
        }
        osList[sI][pI] << endl;
      }
    }
  } else if ((thermoType=="PTXEquilMutation") || (thermoType=="PTXEquilTabulated")) {
    osList.resize(elementsMassFractionList->size());
    List<scalarList>& list = *elementsMassFractionList;

    forAll(osList, eI) {
      osList.set(eI, new PtrList<OFstream>(pressureList.size()));

      forAll(osList[eI], pI) {
        word fileName = "output/" + thermoType + "/" +thermoType + "_mix" + name(eI) + "_P" + name(pressureList[pI]) + ".out";
        osList[eI].set(pI,  new OFstream(fileName));
        const wordList elementsName(basicModel_->getElementsList());
        osList[eI][pI] << "//T(K) p(Pa) rho(kg/m3) psi(s2/m2) e(J/kg) alpha(kg/m/s) alphaEff(kg/m/s) sigma(S/m) gamma(-) mu(Pa.s) muEff(Pa.s) ";
        for(int ie = 0; ie < basicModel_->nElements(); ie++) {
          osList[eI][pI] 	<< "Y[" << basicModel_->getElementName(ie) << "](-) ";
        }
        for(int is = 0; is < basicModel_->nSpecies(); is++) {
          osList[eI][pI] 	<< "Y[" << basicModel_->getSpeciesName(is) << "](-) ";
        }
        if (thermoType == "FiniteRateMutation") {
          for(int is = 0; is < basicModel_->nSpecies(); is++) {
            osList[eI][pI] 	<< "omega[" << basicModel_->getSpeciesName(is) << "](kg/m3/s) ";
          }
        }
        osList[eI][pI] << endl;
      }
    }
  } else {
    osList.resize(1);
    osList.set(0, new PtrList<OFstream>(pressureList.size()));

    forAll(pressureList, pI) {
      word fileName = "output/" + thermoType + "/" +thermoType + "_P" + name(pressureList[pI]) + ".out";
      osList[0].set(pI, new OFstream(fileName));
      osList[0][pI] << "//T(K) p(Pa) rho(kg/m3) psi(s2/m2) e(J/kg) alpha(kg/m/s) alphaEff(kg/m/s) sigma(S/m) gamma(-) mu(Pa.s) muEff(Pa.s) ";
      for(int ie = 0; ie < basicModel_->nElements(); ie++) {
        osList[0][pI] 	<< "Y[" << basicModel_->getElementName(ie) << "](-) ";
      }
      for(int is = 0; is < basicModel_->nSpecies(); is++) {
        osList[0][pI] 	<< "Y[" << basicModel_->getSpeciesName(is) << "](-) ";
      }
      if (thermoType == "FiniteRateMutation") {
        for(int is = 0; is < basicModel_->nSpecies(); is++) {
          osList[0][pI] 	<< "omega[" << basicModel_->getSpeciesName(is) << "](kg/m3/s) ";
        }
      }
      osList[0][pI]	<< endl;
    }
  }



  Info << "- pureThermo: writing values in output" << endl;
//- Write the thermo values

  forAll(pressureList, pI) {
    forAll(temperatureList, tI) {
      if (thermoType == "FiniteRateMutation") {
        forAll(osList, mixI) {
          basicModel_->isInConstructor_ = true;
          // only first cell
          basicModel_->p()[0] = pressureList[pI];
          basicModel_->T()[0] = temperatureList[tI];
          List<scalarList>& list = *speciesMassFractionList;

          forAll(basicModel_->speciesMassFractions(), specI) {
            basicModel_->speciesMassFractions()[specI][0] = list[mixI][specI];
          }

          basicModel_->update();
          osList[mixI][pI] 	<< basicModel_->T()[0] << " "
                            << basicModel_->p()[0] << " "
                            << basicModel_->rho()[0] << " "
                            << basicModel_->psi()[0] << " "
                            << basicModel_->e()[0] << " "
                            << basicModel_->alpha()[0] << " "
                            << turbulence->alphaEff()()[0] << " "
                            << basicModel_->sigma()[0] << " "
                            << basicModel_->gamma()[0] << " "
                            << basicModel_->mu()()[0] << " "
                            << turbulence->muEff()()[0] << " ";

          for(int ie = 0; ie < basicModel_->nElements(); ie++) {
            osList[mixI][pI] 	<< basicModel_->elementsMassFractions()[ie][0] << " ";
          }
          for(int is = 0; is < basicModel_->nSpecies(); is++) {
            osList[mixI][pI] 	<< basicModel_->speciesMassFractions()[is][0] << " ";
          }
          for(int is = 0; is < basicModel_->nSpecies(); is++) {
            osList[mixI][pI] 	<< basicModel_->speciesNetMassProductionRates()[is][0] << " ";
          }
          osList[mixI][pI]	<< endl;
        }
      } else if ((thermoType=="PTXEquilMutation") || (thermoType=="PTXEquilTabulated")) {
        forAll(osList, mixI) {
          basicModel_->isInConstructor_ = true;
          // only first cell
          basicModel_->p()[0] = pressureList[pI];
          basicModel_->T()[0] = temperatureList[tI];
          List<scalarList>& list = *elementsMassFractionList;

          forAll(basicModel_->elementsMassFractions(), elemI) {
            basicModel_->elementsMassFractions()[elemI][0] = list[mixI][elemI];
          }

          Info << "Compute " << pressureList[pI] << " " << temperatureList[tI] << " " << list[mixI]<< endl;
          basicModel_->update();
          osList[mixI][pI] 	<< basicModel_->T()[0] << " "
                            << basicModel_->p()[0] << " "
                            << basicModel_->rho()[0] << " "
                            << basicModel_->psi()[0] << " "
                            << basicModel_->e()[0] << " "
                            << basicModel_->alpha()[0] << " "
                            << turbulence->alphaEff()()[0] << " "
                            << basicModel_->sigma()[0] << " "
                            << basicModel_->gamma()[0] << " "
                            << basicModel_->mu()()[0] << " "
                            << turbulence->muEff()()[0] << " ";

          for(int ie = 0; ie < basicModel_->nElements(); ie++) {
            osList[mixI][pI] 	<< basicModel_->elementsMassFractions()[ie][0] << " ";
          }
          for(int is = 0; is < basicModel_->nSpecies(); is++) {
            osList[mixI][pI] 	<< basicModel_->speciesMassFractions()[is][0] << " ";
          }
          osList[mixI][pI]	<< endl;
        }
      } else {
        basicModel_->isInConstructor_ = true;
        // only first cell
        basicModel_->p()[0] = pressureList[pI];
        basicModel_->T()[0] = temperatureList[tI];
        basicModel_->update();
        osList[0][pI]		<< basicModel_->T()[0] << " "
                        << basicModel_->p()[0] << " "
                        << basicModel_->rho()[0] << " "
                        << basicModel_->psi()[0] << " "
                        << basicModel_->e()[0] << " "
                        << basicModel_->alpha()[0] << " "
                        << turbulence->alphaEff()()[0] << " "
                        << basicModel_->sigma()[0] << " "
                        << basicModel_->gamma()[0] << " "
                        << basicModel_->mu()()[0] << " "
                        << turbulence->muEff()()[0] << " ";

        for(int ie = 0; ie < basicModel_->nElements(); ie++) {
          osList[0][pI] 	<< basicModel_->elementsMassFractions()[ie][0] << " ";
        }
        for(int is = 0; is < basicModel_->nSpecies(); is++) {
          osList[0][pI] 	<< basicModel_->speciesMassFractions()[is][0] << " ";
        }
        osList[0][pI]	<< endl;
      }
    }
  }

  if (rhoList.size()>0) {
    if (thermoType == "FiniteRateMutation")
    {}
    else if ((thermoType=="PTXEquilMutation") || (thermoType=="PTXEquilTabulated")) {
      if (elementsMassFractionList->size() == 0) {
        FatalErrorInFunction << "elementsMassFractionList is empty." << exit(FatalError);
      }
      osList.resize(elementsMassFractionList->size());
      List<scalarList>& list = *elementsMassFractionList;

      forAll(osList, eI) {
        osList.set(eI, new PtrList<OFstream>(rhoList.size()));

        forAll(osList[eI], pI) {
          word fileName = "output/" + thermoType + "/" +thermoType +  "_mix" + name(eI) + "_R" + name(rhoList[pI]) + ".out";
          osList[eI].set(pI, new OFstream(fileName));
          osList[eI][pI] << "//rho(kg/m3)  e(J/kg) p(Pa) T(K) psi(s2/m2) alpha(kg/m/s) alphaEff(kg/m/s) sigma(S/m) gamma(-) mu(Pa.s) muEff(Pa.s)" << endl;
          forAll(energyList, tI) {
            basicModel_->isInConstructor_ = false;
            // only first cell
            basicModel_->rho()[0] = rhoList[pI];
            basicModel_->e()[0] = energyList[tI];
            forAll(basicModel_->elementsMassFractions(), elemI) {
              basicModel_->elementsMassFractions()[elemI][0] = list[eI][elemI];
            }
            basicModel_->update();
            osList[eI][pI]		                << basicModel_->rho()[0] << " "
                                              << basicModel_->e()[0] << " "
                                              << basicModel_->p()[0] << " "
                                              << basicModel_->T()[0] << " "
                                              << basicModel_->psi()[0] << " "
                                              << basicModel_->alpha()[0] << " "
                                              << turbulence->alphaEff()()[0] << " "
                                              << basicModel_->sigma()[0] << " "
                                              << basicModel_->gamma()[0] << " "
                                              << basicModel_->mu()()[0] << " "
                                              << turbulence->muEff()()[0] << endl;
          }
        }
      }
    } else {
      osList.resize(1);
      osList.set(0, new PtrList<OFstream>(rhoList.size()));
      forAll(rhoList, pI) {
        word fileName = "output/" + thermoType + "/" +thermoType + "_R" + name(rhoList[pI]) + ".out";
        osList[0].set(pI, new OFstream(fileName));
        osList[0][pI] << "//rho(kg/m3)  e(J/kg) p(Pa) T(K) psi(s2/m2) alpha(kg/m/s) alphaEff(kg/m/s) sigma(S/m) gamma(-) mu(Pa.s) muEff(Pa.s)" << endl;

        forAll(energyList, tI) {
          basicModel_->isInConstructor_ = false;
          // only first cell
          basicModel_->rho()[0] = rhoList[pI];
          basicModel_->e()[0] = energyList[tI];
          basicModel_->update();
          osList[0][pI]		                << basicModel_->rho()[0] << " "
                                          << basicModel_->e()[0] << " "
                                          << basicModel_->p()[0] << " "
                                          << basicModel_->T()[0] << " "
                                          << basicModel_->psi()[0] << " "
                                          << basicModel_->alpha()[0] << " "
                                          << turbulence->alphaEff()()[0] << " "
                                          << basicModel_->sigma()[0] << " "
                                          << basicModel_->gamma()[0] << " "
                                          << basicModel_->mu()()[0] << " "
                                          << turbulence->muEff()()[0] << endl;
        }
      }

    }
  }
  Info << "- pureThermo: Finish " << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pureThermo::~pureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pureThermo::updateBefore()
{
}

void Foam::pureThermo::updateAfter()
{
}

// ************************************************************************* //

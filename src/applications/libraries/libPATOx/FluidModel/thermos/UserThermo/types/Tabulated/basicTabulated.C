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

#include "basicTabulated.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(basicTabulated, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicTabulated::basicTabulated
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicUserThermo(mesh, dictName),
thermoDict_
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
),
tableName_(thermoDict_.lookup("tableName")),
tableFileName_(changeEnviVar(thermoDict_.path()+"/"+tableName_))
{
  // Get the elements/species name of the mixture
  basicUserThermo::mixtureName_= thermoDict_.lookupOrDefault<word>("mutationMixture", "air13");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicTabulated::~basicTabulated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicTabulated::getElementsComposition(scalarList& Xe, scalarList& Ye, scalarList& rhoe)
{
}

// ====================================================== //

void Foam::basicTabulated::getSpeciesComposition(scalarList& Xs, scalarList& Ys, scalarList& rhos)
{
}

// ====================================================== //

void Foam::basicTabulated::getMixtureR(scalar& mixR)
{
  scalar mixtureMM;
  this->getMixtureMolarMass(mixtureMM);

  mixR = this->RU() / mixtureMM;
}

// ====================================================== //

void Foam::basicTabulated::getMixtureMolarMass(scalar& mixMM)
{
  //- ATM: interpolate MM TO-DO
}

//// ====================================================== //

//void Foam::basicTabulated::getE(scalar& e)
//{
//	//- ATM: interpolate energy TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getH(scalar& h)
//{
//	//- ATM: interpolate enthalpy TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getT(scalar& T)
//{
//	//- ATM: interpolate T TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getP(scalar& p)
//{
//	//- ATM: interpolate p TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getRho(scalar& rho)
//{
//}

//// ====================================================== //

//void Foam::basicTabulated::getPsi(scalar& psi)
//{
//	scalar rho;
//	this->getRho(rho);

//	scalar p;
//	this->getP(p);

//	psi = rho / p;
//}

//// ====================================================== //

//void Foam::basicTabulated::getMu(scalar& mu)
//{
//	//- ATM: interpolate mu TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getNu(scalar& nu)
//{
//	scalar mu;
//	this->getMu(mu);

//	scalar rho;
//	this->getRho(rho);

//	nu = mu / rho;
//}

//// ====================================================== //

//void Foam::basicTabulated::getKappa(scalar& kappa)
//{
//	//- ATM: interpolate kappa TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getCv(scalar& cv)
//{
//    //- ATM: interpolate cv TO-DO
//}

//// ====================================================== //

//void Foam::basicTabulated::getGamma(scalar& gamma)
//{
//	//- ATM: interpolate gamma TO-DO
//}


// ====================================================== //

Foam::fileName Foam::basicTabulated::changeEnviVar(const fileName filenameOrigin_)
{
  // Change environment variable from filename_
  wordList enviVariables_;
  fileName fileNameNew_ = filenameOrigin_;
  fileName fileName_ = filenameOrigin_;
  int i = 0;
  while(i < 100) {
    i++;
    int first = fileNameNew_.find("$");
    int last = fileNameNew_.find("/");
    if (first < 0 ) {
      break;
    }
    string strNew = fileNameNew_.substr (first+1,last-first-1);
    enviVariables_.append(strNew);
    fileNameNew_=fileNameNew_.replaceAll( fileNameNew_.substr (first,last-first),"");

  }

  forAll(enviVariables_, enviI) {
    fileName envi_dir = getEnv(enviVariables_[enviI]);
    fileName_.replaceAll("$"+enviVariables_[enviI] , envi_dir);
  }
  return fileName_;
}


// ************************************************************************* //

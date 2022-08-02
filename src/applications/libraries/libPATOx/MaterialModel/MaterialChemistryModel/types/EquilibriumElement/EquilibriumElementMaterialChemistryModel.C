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
    along with OpenFOAM.  If EquilibriumElementt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "EquilibriumElementMaterialChemistryModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EquilibriumElementMaterialChemistryModel::EquilibriumElementMaterialChemistryModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialChemistryModel(mesh, dictName),
mixtureMutation(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").lookup("mixture")),
mix_(simpleMaterialChemistryModel::mixture_),
tolZ(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").template lookupOrDefault<scalar>("tolZ",1e-10)),
materialPropertiesDirectory(fileName(simpleMaterialChemistryModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh.time().db().parent(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
elementNames_(simpleMaterialChemistryModel::elementNames_),
speciesNames_(simpleMaterialChemistryModel::speciesNames_),
massFractions_(simpleMaterialChemistryModel::massFractions_),
moleFractions_(simpleMaterialChemistryModel::moleFractions_),
moleFractionGasInMaterial_(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").template lookupOrDefault<List<Tuple2<word,scalar> > >("moleFractionGasInMaterial",List<Tuple2<word,scalar> >(0))),
D0_(constantPropertiesDictionary.lookup("D0")),
eta0_(constantPropertiesDictionary.lookup("eta0")),
T_(meshLookupOrConstructScalar(mesh,"Ta")),
Zdefault_
(
    IOobject
    (
        "Zdefault",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh_
),
init_(init()),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
rho_g_(gasPropertiesModel_.rho_g()),
eps_g_(gasPropertiesModel_.eps_g()),
Dm_(gasPropertiesModel_.Dm()),
massModel_(meshLookupOrConstructModel<simpleMassModel>(mesh,dictName,"Mass")),
mDotG_(massModel_.mDotG()), //meshLookupOrConstructVector(mesh, "mDotG", dimensionedVector("zero", dimMass/(dimLength*dimLength*dimTime), vector(0, 0, 0)))),
mDotGFace_(massModel_.mDotGFace()), //meshLookupOrConstructSurfaceVector(mesh, "mDotGface", linearInterpolate(mDotG_))),
mDotGw_(meshLookupOrConstructScalar(mesh, "mDotGw",dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, scalar(0.0)))),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
pi_(pyrolysisModel_.pi()),
fvOptions_(mesh_)
{
  forAll(T_.boundaryField(), patchI) {
    if (isA<BprimeFvPatchScalarField>(T_.boundaryField()[patchI])) {
      BprimeFvPatchScalarField& T_bf_ = refCast<BprimeFvPatchScalarField>(T_.boundaryFieldRef()[patchI]);
      // Verify that the Bprime patch is using the "on-the-fly" option
      if (T_bf_.BprimeBoundaryConditions_.nameBprimeFile_==word::null) {
        BprimePatches_.append(patchI);
      }
    }
  }

  forAll(BprimePatches_, i) {
    label patchI = BprimePatches_[i];
    BprimeFvPatchScalarField& T_bf_ = refCast<BprimeFvPatchScalarField>(T_.boundaryFieldRef()[patchI]);
    BprimeBC_.append(&T_bf_.BprimeBoundaryConditions_);
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- ZmassFlux_ and Dm_ --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }

  ZmassFlux_.resize(ne_mix);
  Dm_.resize(ne_mix);
  forAll(ZmassFlux_, elementI) {
    ZmassFlux_.set
    (
        elementI,
        new surfaceVectorField
        (
            IOobject
            (
                "ZmassFlux[" + elementNames_[elementI] + "]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mDotGFace_
        )
    );
    Dm_.set
    (
        elementI,
        new volScalarField
        (
            IOobject
            (
                "Dm[" + elementNames_[elementI] + "]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            D0_
        )
    );
  }

  // verify pyrolysis reaction rates
  if (pi_.size() != ne_mix) {
    FatalErrorInFunction << "Elemental pyrolysis reaction rates size(" << pi_.size() <<") != elements mixture size(" << ne_mix << ")" << exit(FatalError);
  }

  // add the fields in the interpolation scheme table
  forAll(massFractions_, elementI) {
    fields.add(massFractions_[elementI]);
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- end --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::EquilibriumElementMaterialChemistryModel::~EquilibriumElementMaterialChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EquilibriumElementMaterialChemistryModel::update()
{

  // update the gradient of the elemental mass fractions
  updateBC();

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- epsphi epsrho_g_ ---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
  }
  // Solve the element conservation equation (equilibrium MaterialChemistry mode)
  tmp<surfaceScalarField> epsphi_tmp = mDotGFace_ & mesh_.Sf(); // Mass flux of pyrolysis gases.
  surfaceScalarField& epsphi = const_cast<surfaceScalarField&>(epsphi_tmp());
  tmp<volScalarField> epsrho_g_tmp = eps_g_ * rho_g_;
  volScalarField& epsrho_g_ = const_cast<volScalarField&>(epsrho_g_tmp());

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- mvConvection ---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
  }
  tmp<fv::convectionScheme<scalar> > mvConvection
  (
      fv::convectionScheme<scalar>::New
      (
          mesh_,
          fields,
          epsphi,
          mesh_.divScheme("div(epsphi,Zi)")
      )
  );


  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- Zt ZtotMassFlux ---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
  }
  tmp<volScalarField> Zt_tmp = 0.0 * massFractions_[0];
  volScalarField& Zt = const_cast<volScalarField&>(Zt_tmp());
  tmp<surfaceVectorField> ZtotMassFlux_tmp = 0.0 * ZmassFlux_[0];
  surfaceVectorField& ZtotMassFlux = const_cast<surfaceVectorField&>(ZtotMassFlux_tmp());

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- massFractions_ ---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
  }
  forAll(massFractions_, elementI) {
    if(simpleMaterialChemistryModel::debug_) {
      Info<<"--- Dmi,Zi[" << elementI << "]---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
    }
    volScalarField& Zi = massFractions_[elementI];
    volScalarField& Dmi = Dm_[elementI]; // constant diffusion coefficient = D0 in constantProperties dictionary
    if(simpleMaterialChemistryModel::debug_) {
      Info<<"--- pi_[" << elementI << "]---  Foam:5	:equilibriumElementConservationMaterialChemistrySolver<MaterialChemistryModel>::update()" << endl;
    }
    const volScalarField& totalZp = pi_[elementI];

    if(simpleMaterialChemistryModel::debug_) {
      Info<<"--- fvScalarMatrix[" << elementI << "]---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
    }
    if (this->dynamicMesh_) {
      fvScalarMatrix ZiEqn
      (
          fvm::ddt(epsrho_g_, Zi)                    // storage
          - fvm::div(fvc::interpolate(epsrho_g_)*mesh_.phi(), Zi)  // mesh motion correction (ALE)
          + mvConvection->fvmDiv(epsphi, Zi)         // convection
          - fvm::laplacian(Dmi / eta0_.value() * epsrho_g_, Zi) // diffusion
          - totalZp                 // elements produced by the pyrolysis reactions
          == fvOptions_(epsrho_g_, Zi)
      );
      ZiEqn.solve("Zi");
    } else {
      fvScalarMatrix ZiEqn
      (
          fvm::ddt(epsrho_g_, Zi)                    // storage
          + mvConvection->fvmDiv(epsphi, Zi)         // convection
          - fvm::laplacian(Dmi / eta0_.value()  * epsrho_g_, Zi) // diffusion
          - totalZp                 // elements produced by the pyrolysis reactions
          == fvOptions_(epsrho_g_, Zi)
      );
      ZiEqn.solve("Zi");
    }

    if(simpleMaterialChemistryModel::debug_) {
      Info<<"--- Zt ZmassFlux_ ZtotMassFlux[" << elementI << "]---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
    }
    Zi.max(0.0);
    forAll(Zi, cellI) {
      if(Zi[cellI]<tolZ) {
        Zi[cellI]=0;
      }
    }
    forAll(Zi.boundaryField(), patchI) {
      forAll(Zi.boundaryField()[patchI],faceI) {
        if(Zi.boundaryField()[patchI][faceI]<tolZ) {
          Zi.boundaryFieldRef()[patchI][faceI]=0;
        }
      }
    }
    Zt += Zi;

    // compute the mass diffusion fluxes of the elements
    // convection - implement as zeroGradient, ie, we need the sub surface cell value at the boundary.
    ZmassFlux_[elementI] =  -  fvc::interpolate(Dmi / eta0_.value() * epsrho_g_) * (mesh_.Sf() / mesh_.magSf()) * fvc::snGrad(Zi); //diffusion
    ZtotMassFlux += ZmassFlux_[elementI]; // average mass diffusion flux
  }
  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- massFractions_ ---  Foam::equilibriumElementConservationMaterialChemistryModel::update()" << endl;
    forAll(Zt, cellI) {
      if (Zt[cellI]==0) {
        Zt.write();
        forAll(massFractions_, elementI) {
          massFractions_[elementI].write();
        }
        FatalErrorInFunction << "Zt["<< cellI << "]=0"<< exit(FatalError);
      }
    }

    forAll(this->mesh_.boundaryMesh(), patchI) {
      forAll(this->mesh_.boundaryMesh()[patchI], faceI) {
        if (Zt.boundaryField()[patchI][faceI]==0) {
          Zt.write();
          forAll(massFractions_, elementI) {
            massFractions_[elementI].write();
          }
          FatalErrorInFunction << "Zt.bf["<< patchI << "][" << faceI << "]=0"<< exit(FatalError);
        }
      }
    }
  }

  forAll(massFractions_, elementI) {
    massFractions_[elementI] /= Zt;
  }

  forAll(T_, cellI) {
    double * Z_Y_to_X_ = new double[ne_mix];

    forAll(massFractions_, elementI) {
      Z_Y_to_X_[elementI] = massFractions_[elementI][cellI];
    }
    mix_().convert<Mutation::Thermodynamics::YE_TO_XE>(Z_Y_to_X_, Z_Y_to_X_);
    forAll(moleFractions_, elementI) {
      moleFractions_[elementI][cellI] = Z_Y_to_X_[elementI];
    }
  }

  forAll(T_.boundaryField(), patchI) {
    forAll(T_.boundaryField()[patchI], faceI) {
      double * Z_Y_to_X_ = new double[ne_mix];

      forAll(massFractions_, elementI) {
        Z_Y_to_X_[elementI] = massFractions_[elementI].boundaryField()[patchI][faceI];
      }
      mix_().convert<Mutation::Thermodynamics::YE_TO_XE>(Z_Y_to_X_, Z_Y_to_X_);
      forAll(moleFractions_, elementI) {
        moleFractions_[elementI].boundaryFieldRef()[patchI][faceI] = Z_Y_to_X_[elementI];
      }
    }
  }
}

Switch Foam::EquilibriumElementMaterialChemistryModel::init()
{
  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- init MaterialChemistry --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }

  if (!isFile(constantPropertiesDictionary.path()+"/constantProperties")) {
    FatalErrorInFunction << "Unknown materialPropertiesDirectory " << constantPropertiesDictionary.path()
                         << exit(FatalError);
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- mixture --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }

  optEq_.reset(new Mutation::MixtureOptions(mixtureMutation));
  optEq_().setStateModel("Equil");
  mix_.reset(new Mutation::Mixture(optEq_()));

  ne_mix=mix_().nElements();
  elementNames_.resize(ne_mix);

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- elementName --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }
  forAll(elementNames_, elemI) {
    elementNames_[elemI]=mix_().elementName(elemI);
  }
  moleFractions_.resize(ne_mix);
  massFractions_.resize(ne_mix);

  // Initialize Mutation++ to compute constantEquilibrium MaterialChemistry, thermodynamics and transport properties
  // Initialization of Mutation++ for constantEquilibrium MaterialChemistry
  p_Z.resize(ne_mix);      // mass fractions of elements
  p_Z0.resize(ne_mix);     // ref mass fractions of elements
//  p_Zsub.resize(ne_mix);   // mass fractions of elements in the subsurface cell
  p_Zx.resize(ne_mix);     // mole fractions of elements

  if (moleFractionGasInMaterial_.size()>0) {
    if(simpleMaterialChemistryModel::debug_) {
      Info<<"--- moleFractionGasInMaterial_ --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
    }
    Info << "Use \"Zx[...]\" from \"moleFractionGasInMaterial\" instead of \"" << (word) constantPropertiesDictionary.path() << "/constantProperties\"." << endl;

    scalar sumZx_ = 0;
    scalarList values(moleFractionGasInMaterial_.size());
    wordList names(moleFractionGasInMaterial_.size());

    forAll(moleFractionGasInMaterial_, gasI) {
      names[gasI] = "Zx["+moleFractionGasInMaterial_[gasI].first()+"]";
      values[gasI] = moleFractionGasInMaterial_[gasI].second();
      sumZx_ += values[gasI];
    }
    if (sumZx_ != 1) {
      forAll(moleFractionGasInMaterial_, gasI) {
        values[gasI]/=sumZx_;
      }
    }
    forAll(moleFractionGasInMaterial_, gasI) {
      Info << names[gasI] << " = " << values[gasI] << endl;
      constantPropertiesDictionary.add(names[gasI] ,values[gasI]);
    }
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- read from constantPropertiesDictionary --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }
  Info << "Reading initial elemental composition of the gas in the sample" << nl;
  scalar sumZx_ = 0;
  for (int i = 0; i < ne_mix; i++) {
    p_Zx[i] =
        constantPropertiesDictionary.lookupOrDefault<scalar>
        (
            "Zx[" + mix_().elementName(i) + "]",
            0.0
        );
    Info << "Zx[" << mix_().elementName(i) << "] = " << p_Zx[i] << nl;
    p_Z[i] = p_Zx[i];
    sumZx_+=p_Zx[i];
  }
  if (sumZx_ != 1) {
    for (int i = 0; i < ne_mix; i++) {
      p_Zx[i]/=sumZx_;
      p_Z[i] = p_Zx[i];
    }
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t copyZ" << endl;
  }
  scalar copy_p_Z[ne_mix];
  forAll(p_Z, elemI) {
    copy_p_Z[elemI]=p_Z[elemI];
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- mix XE_TO_YE --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }
  mix_().convert<Mutation::Thermodynamics::XE_TO_YE>(copy_p_Z, copy_p_Z);  // convert in mass fractions

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- massFractions_ and moleFractions_ --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }
  massFractions_.resize(elementNames_.size());
  moleFractions_.resize(elementNames_.size());
  forAll(massFractions_, elemI) {
    massFractions_.set
    (
        elemI,
        new volScalarField
        (
            IOobject
            (
                "Z["+ elementNames_[elemI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Zdefault_
        )
    );
    moleFractions_.set
    (
        elemI,
        new volScalarField
        (
            IOobject
            (
                "Zx["+ elementNames_[elemI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Zdefault_
        )
    );
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- p_Z --- Foam::equilibriumElementConservationMaterialChemistryModel::equilibriumElementConservationMaterialChemistrySolver" << endl;
  }
  forAll(p_Z, elemI) {
    moleFractions_[elemI] == p_Zx[elemI];
    moleFractions_[elemI].boundaryFieldRef() ==  p_Zx[elemI];
    p_Z[elemI]=copy_p_Z[elemI];
    p_Z0[elemI] = p_Z[elemI];
    massFractions_[elemI] == p_Z[elemI];
    massFractions_[elemI].boundaryFieldRef() ==  p_Z[elemI];
  }

  return true;
}

void Foam::EquilibriumElementMaterialChemistryModel::updateBC()
{
  // Imposing fixedGradient on the movingPatch: Z[i]
  // note: we use the special BC 'positiveGradient'
  // gradient() = (Ziw - Ziw.patchInternalField()) * this->patch().deltaCoeffs() * pos((Ziw - Ziw.patchInternalField()));
  for (int elemI = 0; elemI < ne_mix; elemI++) {
    forAll(mesh_.boundaryMesh(), patchJ) {

      fvPatchScalarField& Zi = massFractions_[elemI].boundaryFieldRef()[patchJ];
      if (isA<fixedGradientFvPatchScalarField>(Zi)) {
        fixedGradientFvPatchScalarField& zMovingPatch =
            refCast<fixedGradientFvPatchScalarField>
            (
                Zi
            );

        // update the gradient at the wall
        scalarField zPatchSize(zMovingPatch.size());

        if (!foundInList(patchJ,BprimePatches_)) {
          // updating BL values with constant elemental mass fraction
          forAll(zPatchSize, faceI) {
            zMovingPatch[faceI] = p_Z0[elemI];
          }
        } else { //  all the patches with Bprime
          if(simpleMaterialChemistryModel::debug_) {
            Info<<"--- Ykg_vol ---  Foam::equilibriumElementConservationMaterialChemistryModel::updateBC()" << endl;
          }
          // Update Ykg_vol
          forAll(BprimeBC_, i) {
            label patchI = BprimePatches_[i];

            tmp<vectorField> nf_tmp =
                - mesh_.Sf().boundaryField()[patchI]
                / mesh_.magSf().boundaryField()[patchI];
            vectorField& nf = const_cast<vectorField&>(nf_tmp());

            int ne_Bprime = BprimeBC_[i]->ne_Bprime;
            forAll(mesh_.boundaryMesh()[patchI], faceI) {
              scalar mDotGw_faceI = mDotGw_.boundaryField()[patchI][faceI];
              // update p_Ykg - from wall mass-fraction of element fluxes
              if (mDotGw_faceI >= 1e-6) { // computing elemental composition when gaz flux is outwards
                scalar p_Ykgt = 0;
                for(int elementI = 0; elementI < ne_Bprime; elementI++) {
                  volScalarField& Ykg_vol = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(BprimeBC_[i]->Ykg_names[elementI]));
                  scalar& p_Ykg = Ykg_vol.boundaryFieldRef()[patchI][faceI];
                  const tmp<scalarField> Zsub_tmp = massFractions_[elementI].boundaryField()[patchI].patchInternalField();
                  const scalarField& Zsub = Zsub_tmp();
                  p_Ykg = (Zsub[faceI]* mDotGw_faceI) + (ZmassFlux_[elementI].boundaryField()[patchI][faceI]& (-nf[i]));
                  if (p_Ykg <= 0) p_Ykg = 0;
                  p_Ykgt += p_Ykg;
                }
                for(int elementI = 0; elementI < ne_Bprime; elementI++) {
                  volScalarField& Ykg_vol = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(BprimeBC_[i]->Ykg_names[elementI]));
                  scalar& p_Ykg = Ykg_vol.boundaryFieldRef()[patchI][faceI];
                  p_Ykg = p_Ykg/ p_Ykgt;
                }
              } else { // just grab constant input values
                for(int elementI = 0; elementI < ne_Bprime; elementI++) {
                  volScalarField& Ykg_vol = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(BprimeBC_[i]->Ykg_names[elementI]));
                  scalar& p_Ykg = Ykg_vol.boundaryFieldRef()[patchI][faceI];
                  p_Ykg = p_Z0[elementI];
                }
              } // flux null or fixed BC value
            }
          }

          if(simpleMaterialChemistryModel::debug_) {
            Info<<"--- T_.correctBoundaryConditions(); ---  Foam::equilibriumElementConservationMaterialChemistryModel::updateBC()" << endl;
          }
          // Correct the temperature (Bprime patches)
          T_.correctBoundaryConditions();

          if(simpleMaterialChemistryModel::debug_) {
            Info<<"--- Z_bc_ ---  Foam::equilibriumElementConservationMaterialChemistryModel::updateBC()" << endl;
          }
          // Update the elemental mass fraction at the patches with Bprime
          forAll(BprimeBC_, i) {
            label patchI = BprimePatches_[i];
            int ns_Bprime = BprimeBC_[i]->ns_Bprime;
            int ne_Bprime = BprimeBC_[i]->ne_Bprime;
            forAll(mesh_.boundaryMesh()[patchI], faceI) {
              scalar mwg = 0.0;
              double * Z_bc_ = new double[ne_Bprime]; // elemental mass fraction at the boundary
              for (int elementI = 0; elementI < ne_Bprime; elementI++) {
                Z_bc_[elementI] = 0;
              }
              for(int speciesI = 0; speciesI < ns_Bprime; speciesI++) {
                if (BprimeBC_[i]->mixBprime_ptr->species(speciesI).phase() == Mutation::Thermodynamics::GAS) {
                  volScalarField& Xw_vol = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(BprimeBC_[i]->Xw_names[speciesI]));
                  const scalar Xw_faceI = Xw_vol.boundaryField()[patchI][faceI];
                  mwg += BprimeBC_[i]->mixBprime_ptr->speciesMw(speciesI) * Xw_faceI;
                  for (int elementI = 0; elementI < ne_Bprime; elementI++) {
                    Z_bc_[elementI] += BprimeBC_[i]->mixBprime_ptr->elementMatrix()(speciesI, elementI) * Xw_faceI;
                  }
                }

              }
              for (int elementI = 0; elementI < ne_Bprime; elementI++) {
                Z_bc_[elementI] *=BprimeBC_[i]->mixBprime_ptr->atomicMass(elementI) / mwg;
              }
              for (int elementI = 0; elementI < ne_Bprime; elementI++) {
                massFractions_[elementI].boundaryFieldRef()[patchI][faceI]= Z_bc_[elementI];
              }
              BprimeBC_[i]->mixBprime_ptr->convert<Mutation::Thermodynamics::YE_TO_XE>(Z_bc_,Z_bc_);
              for (int elementI = 0; elementI < ne_Bprime; elementI++) {
                moleFractions_[elementI].boundaryFieldRef()[patchI][faceI]= Z_bc_[elementI];
              }
            }
          }
        }

      }
    }
  }
}


// ************************************************************************* //

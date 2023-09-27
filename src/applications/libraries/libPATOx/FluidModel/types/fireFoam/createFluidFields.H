// Initialise fluid field pointer lists
pThermo.resize(fluidRegions.size());
rhoFluid.resize(fluidRegions.size());
UFluid.resize(fluidRegions.size());
phiFluid.resize(fluidRegions.size());
gFluid.resize(fluidRegions.size());
hRefFluid.resize(fluidRegions.size());
pRefFluid.resize(fluidRegions.size());
ghFluid.resize(fluidRegions.size());
ghfFluid.resize(fluidRegions.size());
turbulence.resize(fluidRegions.size());
combustion.resize(fluidRegions.size());
radiation.resize(fluidRegions.size());
p_rghFluid.resize(fluidRegions.size());
ph_rghFluid.resize(fluidRegions.size());
KFluid.resize(fluidRegions.size());
dpdtFluid.resize(fluidRegions.size());

initialMassFluid.resize(fluidRegions.size());
MRFfluid.resize(fluidRegions.size());
fluidFvOptions.resize(fluidRegions.size());
fieldsFluid.resize(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
  Info<< "*** Reading fluid mesh thermophysical properties for region "
      << fluidRegions[i].name() << nl << endl;

  Info<< "    Adding to pThermo\n" << endl;

  pThermo.set
  (
      i,
      psiReactionThermo::New(fluidRegions[i]).ptr()
  );

  Info<< "    Adding to rhoFluid\n" << endl;
  rhoFluid.set
  (
      i,
      new volScalarField
      (
          IOobject
          (
              "rho",
              runTime.timeName(),
              fluidRegions[i],
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          pThermo[i].rho()
      )
  );

  Info<< "    Adding to UFluid\n" << endl;
  UFluid.set
  (
      i,
      new volVectorField
      (
          IOobject
          (
              "U",
              runTime.timeName(),
              fluidRegions[i],
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
          ),
          fluidRegions[i]
      )
  );

  Info<< "    Adding to phiFluid\n" << endl;
  phiFluid.set
  (
      i,
      new surfaceScalarField
      (
          IOobject
          (
              "phi",
              runTime.timeName(),
              fluidRegions[i],
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          linearInterpolate(rhoFluid[i]*UFluid[i])
          & fluidRegions[i].Sf()
      )
  );

  Info<< "    Adding MRF\n" << endl;
  MRFfluid.set
  (
      i,
      new IOMRFZoneList(fluidRegions[i])
  );

  Info<< "    Adding to turbulence\n" << endl;
  turbulence.set
  (
      i,
      compressible::turbulenceModel::New
      (
          rhoFluid[i],
          UFluid[i],
          phiFluid[i],
          pThermo[i]
      ).ptr()
  );

  turbulence[i].validate();

  Info<< "    Adding to combustion\n" << endl;
  combustion.set
  (
      i,
      CombustionModel<psiReactionThermo>::New
      (
          pThermo[i],
          turbulence[i]
      ).ptr()
  );

  Info<< "    Adding to gFluid\n" << endl;
  gFluid.set
  (
      i,
      new uniformDimensionedVectorField
      (
          IOobject
          (
              "g",
              runTime.constant(),
              fluidRegions[i],
              IOobject::MUST_READ,
              IOobject::NO_WRITE
          )
      )
  );

  Info<< "    Adding to hRefFluid\n" << endl;
  hRefFluid.set
  (
      i,
      new uniformDimensionedScalarField
      (
          IOobject
          (
              "hRef",
              runTime.constant(),
              fluidRegions[i],
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          dimensionedScalar("hRef", dimLength, 0)
      )
  );

  dimensionedScalar ghRef(- mag(gFluid[i])*hRefFluid[i]);

//  dimensionedScalar ghRef
//  (
//      mag(gFluid[i].value()) > SMALL
//      ? gFluid[i]
//      & (cmptMag(gFluid[i].value())/mag(gFluid[i].value()))*hRefFluid[i]
//      : dimensionedScalar("ghRef", gFluid[i].dimensions()*dimLength, 0)
//  );

  Info<< "    Adding to ghFluid\n" << endl;
  ghFluid.set
  (
      i,
      new volScalarField
      (
          "gh",
          (gFluid[i] & fluidRegions[i].C()) - ghRef
      )
  );

  Info<< "    Adding to ghfFluid\n" << endl;
  ghfFluid.set
  (
      i,
      new surfaceScalarField
      (
          "ghf",
          (gFluid[i] & fluidRegions[i].Cf()) - ghRef
      )
  );

  Info<< "    Adding to pRefFluid\n" << endl;
  pRefFluid.set
  (
      i,
      new uniformDimensionedScalarField
      (
          IOobject
          (
              "pRef",
              runTime.constant(),
              fluidRegions[i],
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          dimensionedScalar("pRef", dimPressure, 0)
      )
  );


  ph_rghFluid.set
  (
      i,
      new volScalarField
      (
          IOobject
          (
              "ph_rgh",
              runTime.timeName(),
              fluidRegions[i],
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
          ),
          fluidRegions[i]
      )
  );


  p_rghFluid.set
  (
      i,
      new volScalarField
      (
          IOobject
          (
              "p_rgh",
              runTime.timeName(),
              fluidRegions[i],
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
          ),
          fluidRegions[i]
      )
  );

  // Force p_rgh to be consistent with p
//  p_rghFluid[i] = pThermo[i].p() - rhoFluid[i]*ghFluid[i];

  fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

  Info<< "    Adding to dpdtFluid\n" << endl;
  dpdtFluid.set
  (
      i,
      new volScalarField
      (
          IOobject
          (
              "dpdt",
              runTime.timeName(),
              fluidRegions[i]
          ),
          fluidRegions[i],
          dimensionedScalar
          (
              "dpdt",
              pThermo[i].p().dimensions()/dimTime,
              0
          )
      )
  );

  initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

  Info<< "    Adding to KFluid\n" << endl;
  KFluid.set
  (
      i,
      new volScalarField
      (
          "K",
          0.5*magSqr(UFluid[i])
      )
  );

  radiation.set
  (
      i,
      radiationModel::New(pThermo[i].T())
  );


  Info<< "    Adding fvOptions\n" << endl;
  fluidFvOptions.set
  (
      i,
      new fv::options(fluidRegions[i])
  );

  fieldsFluid.set
  (
      i,
      new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
  );

  basicMultiComponentMixture& composition = pThermo[i].composition();

  forAll(composition.Y(), j) {
    fieldsFluid[i].add(composition.Y()[i]);
  }

  fieldsFluid[i].add(pThermo[i].he());

  dynamicFvMesh& mesh_ = fluidRegions[i];

  diffY.resize(composition.Y().size());
  forAll(composition.Y(), specI) {
    word diffYName = "diffY["+composition.Y()[specI].name()+"]";
    diffY.set(
        specI,
        new volScalarField
        (
            IOobject
            (
                diffYName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(1,-1,-1,0,0), 0.0)
        )
    );
  }
  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(composition.Y(), specI) {
      word diffYName = "diffY["+composition.Y()[specI].name()+"]";
      forAll(mesh_.boundaryMesh(), patchI) {
        if (isA<coupledMixedFvPatchScalarField>(composition.Y()[specI].boundaryFieldRef()[patchI])) {
          coupledMixedFvPatchScalarField& cm_ =  refCast<coupledMixedFvPatchScalarField>(composition.Y()[specI].boundaryFieldRef()[patchI]);
          word name = composition.Y()[specI].name();
          if (composition.Y()[specI].name()=="C") {
            name = "C(gr)";
          }
          cm_.setTnbrName(name);
          cm_.setKappa(diffYName);

        }
      }
    }
  }
}

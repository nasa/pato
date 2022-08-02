// Initialise fluid field pointer lists
thermoFluid.resize(fluidRegions.size());
reactionFluid.resize(fluidRegions.size());
rhoFluid.resize(fluidRegions.size());
UFluid.resize(fluidRegions.size());
phiFluid.resize(fluidRegions.size());
gFluid.resize(fluidRegions.size());
hRefFluid.resize(fluidRegions.size());
ghFluid.resize(fluidRegions.size());
ghfFluid.resize(fluidRegions.size());
turbulence.resize(fluidRegions.size());
p_rghFluid.resize(fluidRegions.size());
radiation.resize(fluidRegions.size());
KFluid.resize(fluidRegions.size());
dpdtFluid.resize(fluidRegions.size());
//dQFluid.resize(fluidRegions.size());

initialMassFluid.resize(fluidRegions.size());
frozenFlowFluid.resize(fluidRegions.size(), false);

MRFfluid.resize(fluidRegions.size());
fluidFvOptions.resize(fluidRegions.size());
fieldsFluid.resize(fluidRegions.size());


// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
  Info<< "*** Reading fluid mesh thermophysical properties for region "
      << fluidRegions[i].name() << nl << endl;

  Info<< "    Creating combustion model\n" << endl;

  thermoFluid.set
  (
      i,
      rhoReactionThermo::New(fluidRegions[i]).ptr()
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
          thermoFluid[i].rho()
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

  dimensionedScalar ghRef
  (
      mag(gFluid[i].value()) > SMALL
      ? gFluid[i]
      & (cmptMag(gFluid[i].value())/mag(gFluid[i].value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", gFluid[i].dimensions()*dimLength, 0)
  );

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

  Info<< "    Adding to turbulence\n" << endl;
  turbulence.set
  (
      i,
      compressible::turbulenceModel::New
      (
          rhoFluid[i],
          UFluid[i],
          phiFluid[i],
          thermoFluid[i]
      ).ptr()
  );

  reactionFluid.set
  (   i,
      CombustionModel<rhoReactionThermo>::New(thermoFluid[i], turbulence[i]).ptr()
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
  p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

  fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

  radiation.set
  (
      i,
      radiationModel::New(thermoFluid[i].T())
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
              thermoFluid[i].p().dimensions()/dimTime,
              0
          )
      )
  );


//  Info<< "    Adding to dQ\n" << endl;
//  dQFluid.set
//  (
//      i,
//      new volScalarField
//      (
//          IOobject
//          (
//              "dQ",
//              runTime.timeName(),
//              fluidRegions[i],
//              IOobject::NO_READ,
//              IOobject::AUTO_WRITE
//          ),
//          fluidRegions[i],
//          dimensionedScalar
//          (
//              "dQ",
//              dimEnergy/dimTime,
//              0
//          )
//      )
//  );

  const dictionary& pimpleDict =
      fluidRegions[i].solutionDict().subDict("PIMPLE");
  pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

  Info<< "    Adding MRF\n" << endl;
  MRFfluid.set
  (
      i,
      new IOMRFZoneList(fluidRegions[i])
  );

  Info<< "    Adding fvOptions\n" << endl;
  fluidFvOptions.set
  (
      i,
      new fv::options(fluidRegions[i])
  );

  turbulence[i].validate();

  fieldsFluid.set
  (
      i,
      new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
  );

  basicMultiComponentMixture& composition = thermoFluid[i].composition();

  forAll(composition.Y(), j) {
    fieldsFluid[i].add(composition.Y()[i]);
  }

  fieldsFluid[i].add(thermoFluid[i].he());

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

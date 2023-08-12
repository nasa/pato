// Initialise fluid field pointer lists
UFluid.resize(fluidRegions.size());
phiFluid.resize(fluidRegions.size());
pFluid.resize(fluidRegions.size());
turbulence.resize(fluidRegions.size());
transport.resize(fluidRegions.size());
Uf_auto_list.resize(fluidRegions.size());
MRFfluid.resize(fluidRegions.size());
fluidFvOptions.resize(fluidRegions.size());

// Initialise global pressure reference
pRefCell = 0;
pRefValue = 0.0;

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
  Info<< "*** Reading fluid mesh thermophysical properties for region "
      << fluidRegions[i].name() << nl << endl;

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



  Info<< "    Adding to pFluid\n" << endl;
  pFluid.set
  (
      i,
      new volScalarField
      (
          IOobject
          (
              "p",
              runTime.timeName(),
              fluidRegions[i],
              IOobject::NO_READ,
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
          fvc::flux(UFluid[i])
      )
  );

  const dictionary& pimpleDict =
      fluidRegions[i].solutionDict().subDict("PIMPLE");


  setRefCell(pFluid[i], pimpleDict, pRefCell, pRefValue);
  fluidRegions[i].setFluxRequired(pFluid[i].name());

  Info<< "    Adding to transport\n" << endl;
  transport.set
  (
      i,
      new singlePhaseTransportModel
      (
          UFluid[i],
          phiFluid[i]
      )
  );

  Info<< "    Adding to turbulence\n" << endl;
  turbulence.set
  (
      i,
      incompressible::turbulenceModel::New
      (
          UFluid[i],
          phiFluid[i],
          transport[i]
      ).ptr()
  );

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

  Info<< "Constructing face velocity Uf\n" << endl;
  Uf_auto_list[i].reset
  (
      new surfaceVectorField
      (
          IOobject
          (
              "Uf"+name(i),
              runTime.timeName(),
              fluidRegions[i],
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          fvc::interpolate(UFluid[i])
      )
  );
}















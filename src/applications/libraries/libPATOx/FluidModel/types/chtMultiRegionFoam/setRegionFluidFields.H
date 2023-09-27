//  fvMesh& mesh = fluidRegions[i];
dynamicFvMesh& mesh = fluidRegions[i];

rhoThermo& thermo = thermoFluid[i];
thermo.validate("PATOx", "h", "e"); // args.executable() = PATOx

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];

compressible::turbulenceModel& turb = turbulence[i];
volScalarField& K = KFluid[i];
volScalarField& dpdt = dpdtFluid[i];

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

volScalarField& p_rgh = p_rghFluid[i];

const dimensionedVector& g = gFluid[i];
const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

radiationModel& rad = radiation[i];

IOMRFZoneList& MRF = MRFfluid[i];
fv::options& fvOptions = fluidFvOptions[i];

const dimensionedScalar initialMass
(
    "initialMass",
    dimMass,
    initialMassFluid[i]
);

const bool frozenFlow = frozenFlowFluid[i];

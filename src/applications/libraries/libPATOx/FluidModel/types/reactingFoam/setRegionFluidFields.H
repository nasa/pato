//  fvMesh& mesh = fluidRegions[i];
dynamicFvMesh& mesh = fluidRegions[i];

rhoReactionThermo& thermo = thermoFluid[i];
CombustionModel<rhoReactionThermo>& reaction = reactionFluid[i];
basicMultiComponentMixture& composition = thermoFluid[i].composition();
PtrList<volScalarField>& Y = composition.Y();
thermo.validate("PATOx", "h", "e"); //args.executable(), "h", "e");

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
//volScalarField& dQ = dQFluid[i];

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
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields = fieldsFluid[i];

const dimensionedScalar initialMass
(
    "initialMass",
    dimMass,
    initialMassFluid[i]
);

const bool frozenFlow = frozenFlowFluid[i];

const word inertSpecie(thermo.lookup("inertSpecie"));


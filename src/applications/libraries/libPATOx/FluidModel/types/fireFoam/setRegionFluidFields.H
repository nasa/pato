dynamicFvMesh& mesh = fluidRegions[i];

psiReactionThermo& thermo = pThermo[i];
thermo.validate("PATOx", "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
  FatalIOErrorIn("PATOx", thermo)
      << "Inert specie " << inertSpecie << " not found in available species "
      << composition.species()
      << exit(FatalIOError);
}

//const label inertIndex(composition.species()[inertSpecie]);

const dimensionedScalar initialMass
(
    "initialMass",
    dimMass,
    initialMassFluid[i]
);

volScalarField& rho = rhoFluid[i];
volScalarField& p = thermo.p();
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];

IOMRFZoneList& MRF = MRFfluid[i];

compressible::turbulenceModel& turb = turbulence[i];
CombustionModel<psiReactionThermo>& comb = combustion[i];

const dimensionedVector& g = gFluid[i];
const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

const dimensionedScalar& pRef = pRefFluid[i];
volScalarField& p_rgh = p_rghFluid[i];
volScalarField& ph_rgh = ph_rghFluid[i];

#include "phrghEqn.H"

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields = fieldsFluid[i];

volScalarField& dpdt = dpdtFluid[i];
volScalarField& K = KFluid[i];
const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();

radiationModel& rad = radiation[i];

fv::options& fvOptions = fluidFvOptions[i];

//const dictionary& additionalControlsDict = mesh.solutionDict().subDict("additionalControls");
//
//Switch solvePrimaryRegion
//(
//    additionalControlsDict.lookupOrDefault<bool>("solvePrimaryRegion", true)
//);
//
//Switch solvePyrolysisRegion
//(
//    additionalControlsDict.lookupOrDefault<bool>("solvePyrolysisRegion", true)
//);

#include "createClouds.H"
#include "createSurfaceFilmModel.H"

regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();

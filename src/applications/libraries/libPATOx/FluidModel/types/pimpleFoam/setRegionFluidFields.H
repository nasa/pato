dynamicFvMesh& mesh = fluidRegions[i];

singlePhaseTransportModel& transp = transport[i];

volVectorField& U = UFluid[i];
autoPtr<surfaceVectorField>& Uf=Uf_auto_list[i];
surfaceScalarField& phi = phiFluid[i];

incompressible::turbulenceModel& turb = turbulence[i];

volScalarField& p = pFluid[i];

IOMRFZoneList& MRF = MRFfluid[i];
fv::options& fvOptions = fluidFvOptions[i];

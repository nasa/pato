#include "setMotionFluid.H"

if (finalIter)
{
  mesh.data::add("finalIteration", true);
}

parcels.evolve();
surfaceFilm.evolve();

if (oCorr == 0)
{
#include "fireFoam_rhoEqn.H"
}

#include "UEqn.H"
#include "YEEqn.H"

// --- PISO loop
for (int corr=0; corr<nCorr; corr++)
{
#include "pEqn.H"
}

turb.correct();

rho = thermo.rho();

if (finalIter)
{
  mesh.data::remove("finalIteration");
}

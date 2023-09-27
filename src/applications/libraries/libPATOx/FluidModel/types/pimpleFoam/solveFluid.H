#include "setMotionFluid.H"

turb.validate();

if (finalIter)
{
  mesh.data::add("finalIteration", true);
}

// --- Pressure-velocity SIMPLE corrector

#include "UEqn.H"

// --- Pressure corrector loop
for (int corr = 0; corr < nCorr; corr++)
{
#include "pEqn.H"
}

/*  TODO : let the user choose when he wants to correct turbulence */
transp.correct();
turb.correct();

if (finalIter)
{
  mesh.data::remove("finalIteration");
}


// ************************************************************************* //


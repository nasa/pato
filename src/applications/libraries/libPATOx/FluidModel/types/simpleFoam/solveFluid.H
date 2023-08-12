#include "setMotionFluid.H"

turb.validate();

// --- Pressure-velocity SIMPLE corrector

#include "UEqn.H"
#include "pEqn.H"


transp.correct();
turb.correct();


// ************************************************************************* //


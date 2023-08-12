/* TEST VERSION*/

// Do any mesh changes
mesh.update();
// Make the fluxes relative to the mesh-motion
fvc::makeRelative(phi, U);

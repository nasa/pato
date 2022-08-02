{

Tw_n = T.boundaryField()[mesh.boundaryMesh().findPatchID(stitchingFluidPatchName)];
Tw_nless1 = T.oldTime().boundaryField()[mesh.boundaryMesh().findPatchID(stitchingFluidPatchName)];
deltaTw = mag( (Tw_n - Tw_nless1) / (mesh.time()).deltaT().value() );

}

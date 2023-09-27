Info << "Moving mesh using least squares interpolation" << endl;

if (dynamicMesh)   //changer de nom
{
  volPointInterpolation pointInterpolation(mesh);

  // Create point mesh
  pointMesh pMesh(mesh);

  wordList types
  (
      pMesh.boundary().size(),
      calculatedFvPatchVectorField::typeName
  );

  pointVectorField pointDU
  (
      IOobject
      (
          "pointDU",
          runTime.timeName(),
          mesh
      ),
      pMesh,
      dimensionedVector("zero", dimLength, vector::zero),
      types
  );


  tmp<volVectorField> deltaD_tmp = D - D.oldTime();
  const volVectorField& deltaD = deltaD_tmp();

  //Info << "deltaD" << D - D.oldTime() << nl;
  //Info << "D" << D << nl;
  //Info << "D_old" << D.oldTime() << nl;
  pointInterpolation.interpolate(deltaD, pointDU);
  //Info << "pointDU" <<pointDU<< nl;

  const vectorField& pointDUI =
      pointDU.internalField();

  //- Move mesh
  vectorField newPoints = mesh.points();

  forAll (pointDUI, pointI) {
    newPoints[pointI] += pointDUI[pointI];
  }

  // Correct symmetryPlane points

  forAll(mesh.boundaryMesh(), patchI) {
    if (isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchI])) {
      const labelList& meshPoints =
          mesh.boundaryMesh()[patchI].meshPoints();

      vector avgN =
          gAverage(mesh.boundaryMesh()[patchI].pointNormals());

      vector i(1, 0, 0);
      vector j(0, 1, 0);
      vector k(0, 0, 1);

      if (mag(avgN&i) > 0.95) {
        forAll(meshPoints, pI) {
          newPoints[meshPoints[pI]].x() = 0;
        }
      } else if (mag(avgN&j) > 0.95) {
        forAll(meshPoints, pI) {
          newPoints[meshPoints[pI]].y() = 0;
        }
      } else if (mag(avgN&k) > 0.95) {
        forAll(meshPoints, pI) {
          newPoints[meshPoints[pI]].z() = 0;
        }
      }
    }
  }

//#   include "calcUnusedNewPoints.H"

//    twoDPointCorrector twoDCorrector(mesh);
//    twoDCorrector.correctPoints(newPoints);
  mesh.movePoints(newPoints);
//    stressMesh.V00();
  mesh.moving(false);
  /*}
  else
  {
      FatalErrorIn(args.executable())
          << "Negative Jacobian"
          << exit(FatalError);
  } */
}

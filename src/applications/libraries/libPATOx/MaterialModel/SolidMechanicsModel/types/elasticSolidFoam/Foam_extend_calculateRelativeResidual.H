{
  scalar maxDU = gMax(mag(U.internalField() - U.oldTime().internalField()));
  relativeResidual =
  gMax
  (
      mag(U.internalField() - U.prevIter().internalField())
      /(maxDU + SMALL)
  );
}

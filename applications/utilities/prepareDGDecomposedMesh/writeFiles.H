
IOobject dgFacePhysicGaussIO
(
    "dgFacePhysicGauss",
    polyMeshDir,
    runTime,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

OFstream os1(dgFacePhysicGaussIO.objectPath());
dgFacePhysicGaussIO.writeHeader(os1, "vectorListList");
os1 << dgFacePhysicGauss << nl;

IOobject dgFaceJacobianIO
(
    "dgFaceJacobian",
    polyMeshDir,
    runTime,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

OFstream os2(dgFaceJacobianIO.objectPath());
dgFaceJacobianIO.writeHeader(os2, "scalarListList");
os2 << dgFaceJ2D << nl;


#ifndef debugProcessorFaceExchange_H
#define debugProcessorFaceExchange_H

#include "Pstream.H"
#include "vector.H"

namespace debugProcessorFaceExchange
{

inline Foam::label targetFaceID()
{
    const Foam::label procNo = Foam::Pstream::myProcNo();

    if (procNo == 0)
    {
        return 342;
    }

    if (procNo == 1)
    {
        return 341;
    }

    return -1;
}


inline Foam::label targetCellID()
{
    const Foam::label procNo = Foam::Pstream::myProcNo();

    if (procNo == 0)
    {
        return 0;
    }

    if (procNo == 1)
    {
        return 38;
    }

    return -1;
}


inline bool isTargetProc()
{
    return targetFaceID() >= 0;
}


inline void printTargetFaceGeometryOnce
(
    const Foam::dgGeomMesh& dgMesh,
    const char* context
)
{
    static bool printed = false;

    if (printed || !Foam::Pstream::parRun() || !isTargetProc())
    {
        return;
    }

    const Foam::label faceID = targetFaceID();
    const Foam::label cellID = targetCellID();
    const Foam::dgGeomFace& face = *dgMesh.faces()[faceID];

    Foam::Pout
        << "\n[debugProcessorFaceExchange] " << context << Foam::nl
        << "  procNo            : " << Foam::Pstream::myProcNo() << Foam::nl
        << "  target local cell : " << cellID << Foam::nl
        << "  target faceID     : " << faceID << Foam::nl
        << "  patchID           : " << face.patchID() << Foam::nl
        << "  isProcessorPatch  : " << face.isProcessorPatch() << Foam::nl
        << "  nGauss            : " << face.nGauss() << Foam::nl
        << "  physical Gauss / ref-owner / ref-neighbour / J2D:" << Foam::nl;

    const Foam::List<Foam::vector>& physGp = face.physicGaussPoints();
    const Foam::List<Foam::vector>& ownerGp = face.gaussPointsOwner();
    const Foam::List<Foam::vector>& neighGp = face.gaussPointsNeighbor();
    const Foam::List<Foam::scalar>& J2D = face.J2D();

    for (Foam::label g = 0; g < face.nGauss(); ++g)
    {
        Foam::Pout
            << "    gp " << g
            << " : phys=" << physGp[g]
            << ", refOwner=" << ownerGp[g]
            << ", refNeighbour=" << neighGp[g]
            << ", J2D=" << J2D[g] << Foam::nl;
    }

    printed = true;
}


template<class Type>
inline void printFieldPlusValues
(
    const char* context,
    const Foam::dgField<Type>& field
)
{
    if (!Foam::Pstream::parRun() || !isTargetProc())
    {
        return;
    }

    const Foam::label cellID = targetCellID();
    const Foam::label faceID = targetFaceID();
    const Foam::faceGaussField<Type>& faceField =
        field.gaussFields()[cellID].faceField();
    const Foam::label localFaceI = faceField.localFaceIndex(faceID);
    const Foam::dgGeomFace& face = *field.gaussFields()[cellID].dgMesh()->faces()[faceID];

    Foam::Pout
        << "\n[debugProcessorFaceExchange] " << context << Foam::nl
        << "  procNo            : " << Foam::Pstream::myProcNo() << Foam::nl
        << "  field             : " << field.name() << Foam::nl
        << "  target local cell : " << cellID << Foam::nl
        << "  target faceID     : " << faceID << Foam::nl
        << "  localFaceI        : " << localFaceI << Foam::nl;

    for (Foam::label g = 0; g < face.nGauss(); ++g)
    {
        Foam::Pout
            << "    gp " << g
            << " : minus=" << faceField.minusValueOnGlobalFace(faceID, g)
            << ", plus=" << faceField.plusValueOnGlobalFace(faceID, g)
            << Foam::nl;
    }
}


template<class Type>
inline void printFaceFlux
(
    const char* context,
    const Foam::GaussField<Type>& fluxField
)
{
    if
    (
        !Foam::Pstream::parRun()
     || !isTargetProc()
     || fluxField.cellID() != targetCellID()
    )
    {
        return;
    }

    const Foam::label faceID = targetFaceID();
    const Foam::faceGaussField<Type>& faceField = fluxField.faceField();
    const Foam::label localFaceI = faceField.localFaceIndex(faceID);
    const Foam::vector& n = faceField.normals()[localFaceI];

    Foam::Pout
        << "\n[debugProcessorFaceExchange] " << context << Foam::nl
        << "  procNo            : " << Foam::Pstream::myProcNo() << Foam::nl
        << "  target local cell : " << fluxField.cellID() << Foam::nl
        << "  target faceID     : " << faceID << Foam::nl
        << "  localFaceI        : " << localFaceI << Foam::nl
        << "  normal            : " << n << Foam::nl;

    const Foam::dgGeomFace& face = *fluxField.dgMesh()->faces()[faceID];

    for (Foam::label g = 0; g < face.nGauss(); ++g)
    {
        const Foam::vector flux = faceField.fluxOnFace(localFaceI, g);

        Foam::Pout
            << "    gp " << g
            << " : flux=" << flux
            << ", |flux|=" << Foam::mag(flux)
            << ", flux.n=" << (flux & n)
            << Foam::nl;
    }
}


inline void printFaceIntegral
(
    const char* context,
    const Foam::GaussField<Foam::vector>& fluxField
)
{
    if
    (
        !Foam::Pstream::parRun()
     || !isTargetProc()
     || fluxField.cellID() != targetCellID()
    )
    {
        return;
    }

    const Foam::label faceID = targetFaceID();
    const Foam::faceGaussField<Foam::vector>& faceField = fluxField.faceField();
    const Foam::label localFaceI = faceField.localFaceIndex(faceID);
    const Foam::vector& n = faceField.normals()[localFaceI];
    const Foam::dgGeomFace& face = *fluxField.dgMesh()->faces()[faceID];
    const Foam::List<Foam::scalar>& w = face.weights();
    const Foam::List<Foam::scalar>& J2D = face.J2D();

    Foam::scalar faceIntegral = 0.0;

    for (Foam::label g = 0; g < face.nGauss(); ++g)
    {
        faceIntegral += J2D[g] * (faceField.fluxOnFace(localFaceI, g) & n) * w[g];
    }

    Foam::Pout
        << "\n[debugProcessorFaceExchange] " << context << Foam::nl
        << "  procNo            : " << Foam::Pstream::myProcNo() << Foam::nl
        << "  target local cell : " << fluxField.cellID() << Foam::nl
        << "  target faceID     : " << faceID << Foam::nl
        << "  localFaceI        : " << localFaceI << Foam::nl
        << "  normal            : " << n << Foam::nl
        << "  face integral     : " << faceIntegral << Foam::nl;
}


inline void printFaceIntegral
(
    const char* context,
    const Foam::GaussField<Foam::tensor>& fluxField
)
{
    if
    (
        !Foam::Pstream::parRun()
     || !isTargetProc()
     || fluxField.cellID() != targetCellID()
    )
    {
        return;
    }

    const Foam::label faceID = targetFaceID();
    const Foam::faceGaussField<Foam::tensor>& faceField = fluxField.faceField();
    const Foam::label localFaceI = faceField.localFaceIndex(faceID);
    const Foam::vector& n = faceField.normals()[localFaceI];
    const Foam::dgGeomFace& face = *fluxField.dgMesh()->faces()[faceID];
    const Foam::List<Foam::scalar>& w = face.weights();
    const Foam::List<Foam::scalar>& J2D = face.J2D();

    Foam::vector faceIntegral(Foam::Zero);

    for (Foam::label g = 0; g < face.nGauss(); ++g)
    {
        faceIntegral += J2D[g] * faceField.fluxOnFace(localFaceI, g) * w[g];
    }

    Foam::Pout
        << "\n[debugProcessorFaceExchange] " << context << Foam::nl
        << "  procNo            : " << Foam::Pstream::myProcNo() << Foam::nl
        << "  target local cell : " << fluxField.cellID() << Foam::nl
        << "  target faceID     : " << faceID << Foam::nl
        << "  localFaceI        : " << localFaceI << Foam::nl
        << "  normal            : " << n << Foam::nl
        << "  face integral     : " << faceIntegral << Foam::nl
        << "  |face integral|   : " << Foam::mag(faceIntegral) << Foam::nl;
}

} // End namespace debugProcessorFaceExchange

#endif

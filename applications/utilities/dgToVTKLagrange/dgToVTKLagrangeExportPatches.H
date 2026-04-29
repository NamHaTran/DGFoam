/**
 * \brief Export only boundary patches as native high-order VTK Lagrange faces.
 */
void exportPatchOnlyTimeStep
(
    const Time& runTime,
    const fvMesh& mesh,
    const dgGeomMesh& dgMesh,
    dgField<scalar>& rhoField,
    dgField<vector>& rhoUField,
    dgField<scalar>& EField,
    const dgThermoConservative& thermo,
    const scalar limiterEpsilon,
    const bool smoothByClampingTheta,
    const scalar clampTheta,
    const volScalarField* theta1FieldPtr,
    const volScalarField* theta2FieldPtr,
    const wordList& patchNames
)
{
    const label pOrder = dgMesh.pOrder();
    const label hoOrder = dgVtkLagrange::exportOrder(pOrder);

    if (hoOrder > 8)
    {
        FatalErrorInFunction
            << "Current dgToVTKLagrange implementation supports "
            << "high-order export up to order 8. "
            << "Detected DG pOrder = " << pOrder
            << ", which maps to VTK order " << hoOrder
            << exit(FatalError);
    }

    if (smoothByClampingTheta && theta1FieldPtr && theta2FieldPtr)
    {
        clampConservativeFieldsByTheta
        (
            dgMesh,
            rhoField,
            rhoUField,
            EField,
            clampTheta,
            *theta1FieldPtr,
            *theta2FieldPtr
        );
    }

    limitConservativeFieldsAtLagrangeNodes
    (
        dgMesh,
        rhoField,
        rhoUField,
        EField,
        limiterEpsilon
    );

    const fvBoundaryMesh& patches = mesh.boundary();
    List<bool> selectedPatches(patches.size(), patchNames.empty());

    forAll(patchNames, patchNameI)
    {
        const word& patchName = patchNames[patchNameI];
        label patchID = -1;

        forAll(patches, patchI)
        {
            if (patches[patchI].name() == patchName)
            {
                patchID = patchI;
                break;
            }
        }

        if (patchID < 0)
        {
            FatalErrorInFunction
                << "Cannot find boundary patch '" << patchName
                << "' requested by -patchList." << exit(FatalError);
        }

        selectedPatches[patchID] = true;
    }

    label skippedFaces = 0;

    forAll(patches, patchI)
    {
        const fvPatch& patch = patches[patchI];

        if (!selectedPatches[patchI])
        {
            continue;
        }

        if (patch.coupled())
        {
            if (!patchNames.empty())
            {
                WarningInFunction
                    << "Skipping coupled patch '" << patch.name()
                    << "' requested by -patchList." << nl;
            }

            continue;
        }

        dgVtkLagrange::GridBuffers grid;
        ExportPointData pointData;

        const label start = patch.start();
        const label end = start + patch.size();

        for (label faceI = start; faceI < end; ++faceI)
        {
            const face& meshFace = mesh.faces()[faceI];
            const label nFacePoints = meshFace.size();

            if (nFacePoints != 3 && nFacePoints != 4)
            {
                ++skippedFaces;
                continue;
            }

            pointField faceVertices(nFacePoints);

            forAll(meshFace, pointI)
            {
                faceVertices[pointI] = mesh.points()[meshFace[pointI]];
            }

            const label ownerCellI = mesh.faceOwner()[faceI];
            const dgGeomCell* ownerCellPtr = dgMesh.cells()[ownerCellI];

            if (!ownerCellPtr)
            {
                continue;
            }

            const dgCellType ownerCellType = ownerCellPtr->type();
            const pointField& ownerCellVertices = ownerCellPtr->points();
            const label pointOffset = grid.points.size();

            std::vector<point> localPoints;
            std::vector<vector> localEtas;
            std::vector<label> nodemap;
            unsigned char vtkCellType = 0u;

            labelList ownerFaceVertexIDs(nFacePoints, -1);

            forAll(faceVertices, pointI)
            {
                ownerFaceVertexIDs[pointI] =
                    findOwnerVertex(faceVertices[pointI], ownerCellVertices);
            }

            if (nFacePoints == 3)
            {
                localPoints.reserve(trianglePointCount(hoOrder));
                localEtas.reserve(trianglePointCount(hoOrder));

                for (label row = 0; row <= hoOrder; ++row)
                {
                    for (label col = 0; col <= hoOrder - row; ++col)
                    {
                        const scalar r = scalar(col)/scalar(hoOrder);
                        const scalar s = scalar(row)/scalar(hoOrder);
                        const vector2D rs
                        (
                            r,
                            s
                        );

                        localPoints.push_back
                        (
                            mapToPhysicalTriangle(rs, faceVertices)
                        );
                        List<scalar> weights(3);
                        weights[0] = 1.0 - r - s;
                        weights[1] = r;
                        weights[2] = s;
                        localEtas.push_back
                        (
                            referencePointEta
                            (
                                ownerCellType,
                                ownerFaceVertexIDs,
                                weights
                            )
                        );
                    }
                }

                nodemap = vtkLagrangeTriangleNodemap(hoOrder);
                vtkCellType = vtkLagrangeTriangleType;
            }
            else
            {
                localPoints.reserve((hoOrder + 1)*(hoOrder + 1));
                localEtas.reserve((hoOrder + 1)*(hoOrder + 1));

                for (label row = 0; row <= hoOrder; ++row)
                {
                    for (label col = 0; col <= hoOrder; ++col)
                    {
                        const scalar eta1 =
                            -1.0 + 2.0*scalar(col)/scalar(hoOrder);
                        const scalar eta2 =
                            -1.0 + 2.0*scalar(row)/scalar(hoOrder);
                        const vector2D eta
                        (
                            eta1,
                            eta2
                        );
                        const scalar N0 = 0.25*(1 - eta1)*(1 - eta2);
                        const scalar N1 = 0.25*(1 + eta1)*(1 - eta2);
                        const scalar N2 = 0.25*(1 + eta1)*(1 + eta2);
                        const scalar N3 = 0.25*(1 - eta1)*(1 + eta2);

                        localPoints.push_back(mapToPhysicalQuad(eta, faceVertices));
                        List<scalar> weights(4);
                        weights[0] = N0;
                        weights[1] = N1;
                        weights[2] = N2;
                        weights[3] = N3;
                        localEtas.push_back
                        (
                            referencePointEta
                            (
                                ownerCellType,
                                ownerFaceVertexIDs,
                                weights
                            )
                        );
                    }
                }

                nodemap = vtkLagrangeQuadNodemap(hoOrder);
                vtkCellType = vtkLagrangeQuadType;
            }

            for (const label localPointI : nodemap)
            {
                const point& x = localPoints[localPointI];
                const vector& eta = localEtas[localPointI];

                List<scalar> basis;
                computeOwnerBasisAt(ownerCellType, eta, pOrder, basis);

                const scalar rhoVal =
                    evaluateModalField(rhoField, ownerCellI, basis);
                const vector rhoUVal =
                    evaluateModalField(rhoUField, ownerCellI, basis);
                const scalar EVal =
                    evaluateModalField(EField, ownerCellI, basis);

                const scalar rhoSafe = max(rhoVal, SMALL);
                const vector UVal = rhoUVal/rhoSafe;
                const scalar eVal = EVal/rhoSafe - 0.5*magSqr(UVal);

                grid.points.push_back(x);
                pointData.rho.push_back(rhoVal);
                pointData.rhoU.push_back(rhoUVal);
                pointData.E.push_back(EVal);
                pointData.U.push_back(UVal);
                pointData.e.push_back(eVal);
                pointData.T.push_back
                (
                    thermo.calcTemperatureFromRhoHe(ownerCellI, rhoSafe, eVal)
                );
                pointData.p.push_back
                (
                    thermo.calcPressureFromRhoHe(ownerCellI, rhoSafe, eVal)
                );
                pointData.a.push_back
                (
                    thermo.calcSpeedOfSoundFromRhoHe(ownerCellI, rhoSafe, eVal)
                );
            }

            for (label localPointI = 0; localPointI < label(nodemap.size()); ++localPointI)
            {
                grid.connectivity.push_back(pointOffset + localPointI);
            }

            grid.offsets.push_back(grid.connectivity.size());
            grid.cellTypes.push_back(vtkCellType);
            grid.partitions.push_back(Pstream::myProcNo());
        }

        const std::vector<dgVtkLagrange::ScalarPointData> scalarArrays
        {
            {"rho", &pointData.rho},
            {"E", &pointData.E},
            {"e", &pointData.e},
            {"p", &pointData.p},
            {"T", &pointData.T},
            {"a", &pointData.a}
        };

        const std::vector<dgVtkLagrange::VectorPointData> vectorArrays
        {
            {"rhoU", &pointData.rhoU},
            {"U", &pointData.U}
        };

        ensureOutputDir(runTime);

        fileName baseName = patch.name() + "_vtkLagrange";

        if (Pstream::parRun())
        {
            baseName += "_" + Foam::name(Pstream::myProcNo());
        }

        dgVtkLagrange::writeVtuPiece
        (
            outputDir(runTime)/(baseName + ".vtu"),
            runTime.value(),
            grid,
            scalarArrays,
            vectorArrays
        );
    }

    if (skippedFaces > 0)
    {
        WarningInFunction
            << "Skipped " << skippedFaces
            << " boundary face(s) because only triangular and quadrilateral "
            << "patch faces are supported for -no-internal." << nl;
    }
}

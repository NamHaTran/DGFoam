#include "dgGaussFieldVectorMathTmp.H"

{
    Info<< "\n===== testVectorMath (GaussField vector-math tmp) =====\n" << endl;

    if (mesh.C().size() == 0)
    {
        Info<< "Mesh has no cells. Skip vector math tests." << endl;
    }
    else
    {
        const label testCell = 0;

        GaussField<vector> gV1(testCell, &dgMesh);
        GaussField<vector> gV2(testCell, &dgMesh);
        GaussField<tensor> gT(testCell, &dgMesh);
        GaussField<symmTensor> gSym(testCell, &dgMesh);
        GaussField<sphericalTensor> gSph(testCell, &dgMesh);

        gV1 = vector(1.0, 2.0, 3.0);
        gV2 = vector(-1.0, 0.5, 4.0);

        gT = tensor
        (
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0
        );
        gSym = symmTensor(2.0, 0.1, 0.2, 3.0, 0.3, 4.0);
        gSph = sphericalTensor(5.0);

        Info<< "Initial gV1:\n" << gV1 << endl;
        Info<< "Initial gV2:\n" << gV2 << endl;

        // vector - vector
        { tmp<GaussField<vector>> t = gV1 + gV2; Info<< "gV1 + gV2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gV1 - gV2; Info<< "gV1 - gV2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gV1 ^ gV2; Info<< "gV1 ^ gV2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gV1 & gV2; Info<< "gV1 & gV2 (dot):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gV1 * gV2; Info<< "gV1 * gV2 (outer):\n" << t() << endl; t.clear(); }

        // unary
        { tmp<GaussField<scalar>> t = mag(gV1); Info<< "mag(gV1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = magSqr(gV1); Info<< "magSqr(gV1):\n" << t() << endl; t.clear(); }

        // vector & tensor-family
        { tmp<GaussField<vector>> t = gV1 & gT; Info<< "gV1 & gT:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gT & gV1; Info<< "gT & gV1:\n" << t() << endl; t.clear(); }

        { tmp<GaussField<vector>> t = gV1 & gSym; Info<< "gV1 & gSym:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gSym & gV1; Info<< "gSym & gV1:\n" << t() << endl; t.clear(); }

        { tmp<GaussField<vector>> t = gV1 & gSph; Info<< "gV1 & gSph:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gSph & gV1; Info<< "gSph & gV1:\n" << t() << endl; t.clear(); }

        // tmp variants (exercise tmp-aware paths)
        {
            tmp<GaussField<vector>> tA(new GaussField<vector>(gV1));
            tmp<GaussField<vector>> tB(new GaussField<vector>(gV2));
            tmp<GaussField<vector>> tOut = tA + tB;
            Info<< "tmp(gV1) + tmp(gV2):\n" << tOut() << endl;
            tOut.clear();
        }

        {
            tmp<GaussField<vector>> tA(new GaussField<vector>(gV1));
            tmp<GaussField<vector>> tB(new GaussField<vector>(gV2));
            tmp<GaussField<scalar>> tOut = tA & tB;
            Info<< "tmp(gV1) & tmp(gV2):\n" << tOut() << endl;
            tOut.clear();
        }

        {
            tmp<GaussField<vector>> tA(new GaussField<vector>(gV1));
            tmp<GaussField<tensor>> tB(new GaussField<tensor>(gT));
            tmp<GaussField<vector>> tOut = tA & tB;
            Info<< "tmp(gV1) & tmp(gT):\n" << tOut() << endl;
            tOut.clear();
        }

        {
            tmp<GaussField<vector>> tA(new GaussField<vector>(gV1));
            tmp<GaussField<scalar>> tOut = mag(tA);
            Info<< "mag(tmp(gV1)):\n" << tOut() << endl;
            tOut.clear();
        }
    }

    Info<< "===== end testVectorMath =====\n" << endl;
}

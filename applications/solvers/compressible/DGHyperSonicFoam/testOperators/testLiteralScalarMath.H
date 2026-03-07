#include "dgGaussFieldLiteralScalarMathTmp.H"

{
    Info<< "\n===== testLiteralScalarMath (GaussField literal-scalar tmp) =====\n" << endl;

    if (mesh.C().size() == 0)
    {
        Info<< "Mesh has no cells. Skip literal-scalar math tests." << endl;
    }
    else
    {
        const label testCell = 0;

        GaussField<scalar> gS(testCell, &dgMesh);
        GaussField<vector> gV(testCell, &dgMesh);
        GaussField<tensor> gT(testCell, &dgMesh);
        GaussField<symmTensor> gSym(testCell, &dgMesh);
        GaussField<sphericalTensor> gSph(testCell, &dgMesh);

        gS = 4.0;
        gV = vector(1.0, 2.0, 3.0);
        gT = tensor
        (
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0
        );
        gSym = symmTensor(2.0, 0.1, 0.2, 3.0, 0.3, 4.0);
        gSph = sphericalTensor(5.0);

        Info<< "Initial gS:\n" << gS << endl;
        Info<< "Initial gV:\n" << gV << endl;
        Info<< "Initial gT:\n" << gT << endl;
        Info<< "Initial gSym:\n" << gSym << endl;
        Info<< "Initial gSph:\n" << gSph << endl;

        // * and / : scalar with scalar/vector/tensor/symmTensor/sphericalTensor
        { tmp<GaussField<scalar>> t = 2.0 * gS; Info<< "2.0 * gS:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS * 2.0; Info<< "gS * 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = 2.0 * gV; Info<< "2.0 * gV:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gV * 2.0; Info<< "gV * 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = 2.0 * gT; Info<< "2.0 * gT:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT * 2.0; Info<< "gT * 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = 2.0 * gSym; Info<< "2.0 * gSym:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = gSym * 2.0; Info<< "gSym * 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<sphericalTensor>> t = 2.0 * gSph; Info<< "2.0 * gSph:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<sphericalTensor>> t = gSph * 2.0; Info<< "gSph * 2.0:\n" << t() << endl; t.clear(); }

        { tmp<GaussField<scalar>> t = 8.0 / gS; Info<< "8.0 / gS:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS / 2.0; Info<< "gS / 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gV / 2.0; Info<< "gV / 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT / 2.0; Info<< "gT / 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = gSym / 2.0; Info<< "gSym / 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<sphericalTensor>> t = gSph / 2.0; Info<< "gSph / 2.0:\n" << t() << endl; t.clear(); }

        // + and - : scalar with GaussField<scalar>
        { tmp<GaussField<scalar>> t = 2.0 + gS; Info<< "2.0 + gS:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS + 2.0; Info<< "gS + 2.0:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = 2.0 - gS; Info<< "2.0 - gS:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS - 2.0; Info<< "gS - 2.0:\n" << t() << endl; t.clear(); }

        // pow, sqrt
        { tmp<GaussField<scalar>> t = pow(gS, 2.0); Info<< "pow(gS, 2.0):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = sqrt(gS); Info<< "sqrt(gS):\n" << t() << endl; t.clear(); }

        // tmp-input variants
        {
            tmp<GaussField<scalar>> tIn(new GaussField<scalar>(gS));
            tmp<GaussField<scalar>> tOut = 2.0 * tIn;
            Info<< "2.0 * (tmp gS):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<vector>> tIn(new GaussField<vector>(gV));
            tmp<GaussField<vector>> tOut = tIn * 2.0;
            Info<< "(tmp gV) * 2.0:\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<scalar>> tIn(new GaussField<scalar>(gS));
            tmp<GaussField<scalar>> tOut = 2.0 + tIn;
            Info<< "2.0 + (tmp gS):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<scalar>> tIn(new GaussField<scalar>(gS));
            tmp<GaussField<scalar>> tOut = tIn - 1.0;
            Info<< "(tmp gS) - 1.0:\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<scalar>> tIn(new GaussField<scalar>(gS));
            tmp<GaussField<scalar>> tOut = pow(tIn, 3.0);
            Info<< "pow(tmp gS, 3.0):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<scalar>> tIn(new GaussField<scalar>(gS));
            tmp<GaussField<scalar>> tOut = sqrt(tIn);
            Info<< "sqrt(tmp gS):\n" << tOut() << endl;
            tOut.clear();
        }

        // compare (field, scalar) and (scalar, field)
        Info<< "gS > 3.0  : " << (gS > 3.0) << endl;
        Info<< "gS >= 4.0 : " << (gS >= 4.0) << endl;
        Info<< "gS < 5.0  : " << (gS < 5.0) << endl;
        Info<< "gS <= 4.0 : " << (gS <= 4.0) << endl;
        Info<< "gS == 4.0 : " << (gS == 4.0) << endl;
        Info<< "gS != 4.0 : " << (gS != 4.0) << endl;

        Info<< "3.0 > gS  : " << (3.0 > gS) << endl;
        Info<< "4.0 >= gS : " << (4.0 >= gS) << endl;
        Info<< "5.0 < gS  : " << (5.0 < gS) << endl;
        Info<< "4.0 <= gS : " << (4.0 <= gS) << endl;
        Info<< "4.0 == gS : " << (4.0 == gS) << endl;
        Info<< "4.0 != gS : " << (4.0 != gS) << endl;

        {
            tmp<GaussField<scalar>> tCmp(new GaussField<scalar>(gS));
            Info<< "tmp(gS) > 3.0 : " << (tCmp > 3.0) << endl;
            tCmp.clear();
        }
        {
            tmp<GaussField<scalar>> tCmp(new GaussField<scalar>(gS));
            Info<< "3.0 < tmp(gS) : " << (3.0 < tCmp) << endl;
            tCmp.clear();
        }
    }

    Info<< "===== end testLiteralScalarMath =====\n" << endl;
}

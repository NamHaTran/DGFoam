
{
    // Tutorial helper: verifies field-field scalar algebra and tmp-aware paths.
    Info<< "\n===== testScalarMath (GaussField scalar-math tmp) =====\n" << endl;

    if (mesh.C().size() == 0)
    {
        Info<< "Mesh has no cells. Skip scalar math tests." << endl;
    }
    else
    {
        const label testCell = 0;

        GaussField<scalar> gS1(testCell, &dgMesh);
        GaussField<scalar> gS2(testCell, &dgMesh);
        GaussField<vector> gV(testCell, &dgMesh);
        GaussField<tensor> gT(testCell, &dgMesh);
        GaussField<symmTensor> gSym(testCell, &dgMesh);
        GaussField<sphericalTensor> gSph(testCell, &dgMesh);

        gS1 = 4.0;
        gS2 = 2.0;
        gV = vector(1.0, 2.0, 3.0);
        gT = tensor
        (
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0
        );
        gSym = symmTensor(2.0, 0.1, 0.2, 3.0, 0.3, 4.0);
        gSph = sphericalTensor(5.0);

        Info<< "Initial gS1:\n" << gS1 << endl;
        Info<< "Initial gS2:\n" << gS2 << endl;

        // scalar - scalar
        { tmp<GaussField<scalar>> t = gS1 + gS2; Info<< "gS1 + gS2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS1 - gS2; Info<< "gS1 - gS2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS1 * gS2; Info<< "gS1 * gS2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = gS1 / gS2; Info<< "gS1 / gS2:\n" << t() << endl; t.clear(); }

        // scalar - vector
        { tmp<GaussField<vector>> t = gS1 * gV; Info<< "gS1 * gV:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gV * gS1; Info<< "gV * gS1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<vector>> t = gV / gS1; Info<< "gV / gS1:\n" << t() << endl; t.clear(); }

        // scalar - tensor
        { tmp<GaussField<tensor>> t = gS1 * gT; Info<< "gS1 * gT:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT * gS1; Info<< "gT * gS1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT / gS1; Info<< "gT / gS1:\n" << t() << endl; t.clear(); }

        // scalar - symmTensor
        { tmp<GaussField<symmTensor>> t = gS1 * gSym; Info<< "gS1 * gSym:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = gSym * gS1; Info<< "gSym * gS1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = gSym / gS1; Info<< "gSym / gS1:\n" << t() << endl; t.clear(); }

        // scalar - sphericalTensor
        { tmp<GaussField<sphericalTensor>> t = gS1 * gSph; Info<< "gS1 * gSph:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<sphericalTensor>> t = gSph * gS1; Info<< "gSph * gS1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<sphericalTensor>> t = gSph / gS1; Info<< "gSph / gS1:\n" << t() << endl; t.clear(); }

        // tmp variants (exercise tmp-aware paths)
        {
            tmp<GaussField<scalar>> tA(new GaussField<scalar>(gS1));
            tmp<GaussField<scalar>> tB(new GaussField<scalar>(gS2));
            tmp<GaussField<scalar>> tOut = tA + tB;
            Info<< "tmp(gS1) + tmp(gS2):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<vector>> tV(new GaussField<vector>(gV));
            tmp<GaussField<vector>> tOut = gS1 * tV;
            Info<< "gS1 * tmp(gV):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<vector>> tV(new GaussField<vector>(gV));
            tmp<GaussField<scalar>> tS(new GaussField<scalar>(gS1));
            tmp<GaussField<vector>> tOut = tV / tS;
            Info<< "tmp(gV) / tmp(gS1):\n" << tOut() << endl;
            tOut.clear();
        }

        // comparisons: scalar fields
        Info<< "gS1 <  gS2 : " << (gS1 <  gS2) << endl;
        Info<< "gS1 <= gS2 : " << (gS1 <= gS2) << endl;
        Info<< "gS1 >  gS2 : " << (gS1 >  gS2) << endl;
        Info<< "gS1 >= gS2 : " << (gS1 >= gS2) << endl;
        Info<< "gS1 == gS2 : " << (gS1 == gS2) << endl;
        Info<< "gS1 != gS2 : " << (gS1 != gS2) << endl;

        {
            tmp<GaussField<scalar>> tA(new GaussField<scalar>(gS1));
            tmp<GaussField<scalar>> tB(new GaussField<scalar>(gS2));
            Info<< "tmp(gS1) > tmp(gS2) : " << (tA > tB) << endl;
            tA.clear();
            tB.clear();
        }
    }

    Info<< "===== end testScalarMath =====\n" << endl;
}

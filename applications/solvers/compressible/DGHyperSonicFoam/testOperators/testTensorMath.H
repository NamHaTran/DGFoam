#include "dgGaussFieldTensorMathTmp.H"

{
    // Tutorial helper: demonstrates tensor/symmTensor algebra on GaussField.
    Info<< "\n===== testTensorMath (GaussField tensor-math tmp) =====\n" << endl;

    if (mesh.C().size() == 0)
    {
        Info<< "Mesh has no cells. Skip tensor math tests." << endl;
    }
    else
    {
        const label testCell = 0;

        GaussField<tensor> gT1(testCell, &dgMesh);
        GaussField<tensor> gT2(testCell, &dgMesh);
        GaussField<symmTensor> gSym1(testCell, &dgMesh);
        GaussField<symmTensor> gSym2(testCell, &dgMesh);

        gT1 = tensor
        (
            1.0, 2.0, 3.0,
            0.0, 4.0, 5.0,
            0.0, 0.0, 6.0
        );
        gT2 = tensor
        (
            2.0, -1.0, 0.0,
            1.0,  3.0, 1.0,
            4.0,  2.0, 5.0
        );

        gSym1 = symmTensor(3.0, 0.2, 0.1, 4.0, 0.3, 5.0);
        gSym2 = symmTensor(1.0, 0.0, 0.2, 2.0, 0.1, 3.0);

        Info<< "Initial gT1:\n" << gT1 << endl;
        Info<< "Initial gT2:\n" << gT2 << endl;
        Info<< "Initial gSym1:\n" << gSym1 << endl;
        Info<< "Initial gSym2:\n" << gSym2 << endl;

        // tensor unary utilities
        { tmp<GaussField<tensor>> t = dev(gT1); Info<< "dev(gT1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = symm(gT1); Info<< "symm(gT1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = skew(gT1); Info<< "skew(gT1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = T(gT1); Info<< "T(gT1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = tr(gT1); Info<< "tr(gT1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = det(gT1); Info<< "det(gT1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = inv(gT1); Info<< "inv(gT1):\n" << t() << endl; t.clear(); }

        // symmTensor unary utilities
        { tmp<GaussField<scalar>> t = tr(gSym1); Info<< "tr(gSym1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = dev(gSym1); Info<< "dev(gSym1):\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = inv(gSym1); Info<< "inv(gSym1):\n" << t() << endl; t.clear(); }

        // tensor - tensor
        { tmp<GaussField<tensor>> t = gT1 + gT2; Info<< "gT1 + gT2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT1 - gT2; Info<< "gT1 - gT2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT1 & gT2; Info<< "gT1 & gT2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = (gT1 && gT2); Info<< "gT1 && gT2:\n" << t() << endl; t.clear(); }

        // tensor - symmTensor
        { tmp<GaussField<tensor>> t = gT1 + gSym1; Info<< "gT1 + gSym1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT1 - gSym1; Info<< "gT1 - gSym1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gT1 & gSym1; Info<< "gT1 & gSym1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = (gT1 && gSym1); Info<< "gT1 && gSym1:\n" << t() << endl; t.clear(); }

        // symmTensor - tensor
        { tmp<GaussField<tensor>> t = gSym1 + gT1; Info<< "gSym1 + gT1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gSym1 - gT1; Info<< "gSym1 - gT1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gSym1 & gT1; Info<< "gSym1 & gT1:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = (gSym1 && gT1); Info<< "gSym1 && gT1:\n" << t() << endl; t.clear(); }

        // symmTensor - symmTensor
        { tmp<GaussField<symmTensor>> t = gSym1 + gSym2; Info<< "gSym1 + gSym2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<symmTensor>> t = gSym1 - gSym2; Info<< "gSym1 - gSym2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<tensor>> t = gSym1 & gSym2; Info<< "gSym1 & gSym2:\n" << t() << endl; t.clear(); }
        { tmp<GaussField<scalar>> t = (gSym1 && gSym2); Info<< "gSym1 && gSym2:\n" << t() << endl; t.clear(); }

        // tmp variants (exercise tmp-aware paths)
        {
            tmp<GaussField<tensor>> tA(new GaussField<tensor>(gT1));
            tmp<GaussField<tensor>> tB(new GaussField<tensor>(gT2));
            tmp<GaussField<tensor>> tOut = tA + tB;
            Info<< "tmp(gT1) + tmp(gT2):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<tensor>> tA(new GaussField<tensor>(gT1));
            tmp<GaussField<tensor>> tOut = dev(tA);
            Info<< "dev(tmp(gT1)):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<symmTensor>> tA(new GaussField<symmTensor>(gSym1));
            tmp<GaussField<tensor>> tB(new GaussField<tensor>(gT1));
            tmp<GaussField<tensor>> tOut = tA & tB;
            Info<< "tmp(gSym1) & tmp(gT1):\n" << tOut() << endl;
            tOut.clear();
        }
        {
            tmp<GaussField<tensor>> tA(new GaussField<tensor>(gT1));
            tmp<GaussField<scalar>> tOut = det(tA);
            Info<< "det(tmp(gT1)):\n" << tOut() << endl;
            tOut.clear();
        }
    }

    Info<< "===== end testTensorMath =====\n" << endl;
}

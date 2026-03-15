#include "IO.h"
#include <string>
#include "VarDeclaration.h"
#include "DGIOLib.h"
#include "DurstModel.h"
#include <iostream>

namespace IO_Durst {
    void readMassDiffCoefModel()
    {
        //Input do dai cua cac array o day
        const int numOfInt(1),
                numOfDouble(1),
                numOfBool(1),
                numOfStr(1);

        /*Read FlowProperties*/
        std::string fileName("DurstModelProperties.txt");
        std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant");
        std::string keyWordsDouble[numOfDouble] = {"blendingFactor"}, keyWordsInt[numOfInt] = {}, keyWordsBool[numOfBool] = {}, keyWordsStr[numOfStr] = {"massDiffusionCoeff"};
        double outDB[numOfDouble] = {};
        int outInt[numOfInt] = {};
        bool outBool[numOfBool] = {};
        std::string outStr[numOfStr] = {};

        /*NOTE: trong ham readDataFile, 4 argument cuoi cung la so luong variable trong file can phai doc tuong ung voi kieu du lieu:
         * double, int, bool, string.
         * Neu kieu du lieu nao khong co data can doc, them dau '-' vao phia truoc ten bien, vi du:
         * readDataFile(..., numOfDouble, -numOfInt, -numOfBool, -numOfStr) ----> chi doc bien co kieu du lieu double
        */
        IO::readDataFile(fileName, Loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, numOfDouble, -numOfInt, -numOfBool, numOfStr);

        if (outStr[0].compare("ChapmanEnskog")==0)
        {
            extNSF_Durst::massDiffModel_ChapmanEnskog=true;
            extNSF_Durst::massDiffModel_constant=false;

            //Read parameter
            IO_Durst::massDiffCoef_ChapmanEnskog();
        }
        else if (outStr[0].compare("constant")==0)
        {
            extNSF_Durst::massDiffModel_ChapmanEnskog=false;
            extNSF_Durst::massDiffModel_constant=true;

            //Read parameter
            IO_Durst::massDiffCoef_constant();
        }
        else {
            std::string str0("Model '"+outStr[0]+"' of DurstModelProperties>massDiffusionCoeff is not a available.");
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, str0);
        }

        extNSF_Durst::blending=outDB[0];
    }

    void massDiffCoef_ChapmanEnskog()
    {
        //Input do dai cua cac array o day
        const int numOfInt(1),
                numOfDouble(2),
                numOfBool(1),
                numOfStr(1);

        /*Read FlowProperties*/
        std::string fileName("DurstModelProperties.txt");
        std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant");
        std::string keyWordsDouble[numOfDouble] = {"molMass", "collisionDiameter"}, keyWordsInt[numOfInt] = {}, keyWordsBool[numOfBool] = {}, keyWordsStr[numOfStr] = {};
        double outDB[numOfDouble] = {};
        int outInt[numOfInt] = {};
        bool outBool[numOfBool] = {};
        std::string outStr[numOfStr] = {};

        /*NOTE: trong ham readDataFile, 4 argument cuoi cung la so luong variable trong file can phai doc tuong ung voi kieu du lieu:
         * double, int, bool, string.
         * Neu kieu du lieu nao khong co data can doc, them dau '-' vao phia truoc ten bien, vi du:
         * readDataFile(..., numOfDouble, -numOfInt, -numOfBool, -numOfStr) ----> chi doc bien co kieu du lieu double
        */
        IO::readDataFile(fileName, Loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, numOfDouble, -numOfInt, -numOfBool, -numOfStr);

        //Chua add model Chapman-Enskog
    }

    void massDiffCoef_constant()
    {
        //Input do dai cua cac array o day
        const int numOfInt(1),
                numOfDouble(1),
                numOfBool(1),
                numOfStr(1);

        /*Read FlowProperties*/
        std::string fileName("DurstModelProperties.txt");
        std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant");
        std::string keyWordsDouble[numOfDouble] = {"Dm"}, keyWordsInt[numOfInt] = {}, keyWordsBool[numOfBool] = {}, keyWordsStr[numOfStr] = {};
        double outDB[numOfDouble] = {};
        int outInt[numOfInt] = {};
        bool outBool[numOfBool] = {};
        std::string outStr[numOfStr] = {};

        /*NOTE: trong ham readDataFile, 4 argument cuoi cung la so luong variable trong file can phai doc tuong ung voi kieu du lieu:
         * double, int, bool, string.
         * Neu kieu du lieu nao khong co data can doc, them dau '-' vao phia truoc ten bien, vi du:
         * readDataFile(..., numOfDouble, -numOfInt, -numOfBool, -numOfStr) ----> chi doc bien co kieu du lieu double
        */
        IO::readDataFile(fileName, Loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, numOfDouble, -numOfInt, -numOfBool, -numOfStr);

        extNSF_Durst::Dm=outDB[0];
    }

    void showSettings()
    {
        if (extNSF_Durst::enable)
        {
            std::cout<<"Durst model settings:\n";
            std::cout<<"       + Allow Self-Diffusion at Wall: ";
            if (extNSF_Durst::diffusionAtWall)
                std::cout<<"Yes.\n";
            else
                std::cout<<"No.\n";

            std::cout<<"       + Mass diffusion coefficient model: ";
            if (extNSF_Durst::massDiffModel_ChapmanEnskog)
            {
                std::cout<<"Chapman-Enskog.\n";
                std::cout<<"       + Chapman-Enskog's parameters:\n"
                        <<"           molMass           "<<40<<" (g/mol)\n"
                        <<"           collisionDiameter "<<20<<" (angstrom)\n";
            }

            else if (extNSF_Durst::massDiffModel_constant)
            {
                std::cout<<"Constant.\n";
                std::cout<<"       + Constant coefficient:\n"
                         <<"           Dm                "<<extNSF_Durst::Dm<<".\n";
            }
            std::cout<<"       + Blending factor: "<<extNSF_Durst::blending<<".\n";
        }
    }
}

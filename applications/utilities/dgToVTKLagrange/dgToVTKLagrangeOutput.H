/**
 * \brief Return the output directory for the current export time.
 */
fileName outputDir(const Time& runTime)
{
    return runTime.globalPath()/"VTK_Lagrange"/runTime.timeName();
}


/**
 * \brief Return the VTU piece file written by the current MPI rank.
 */
fileName pieceFileName(const Time& runTime)
{
    const fileName dir = outputDir(runTime);

    if (Pstream::parRun())
    {
        return dir/("dgToVTKLagrange_" + Foam::name(Pstream::myProcNo()) + ".vtu");
    }

    return dir/"dgToVTKLagrange.vtu";
}


/**
 * \brief Return the parallel PVTU manifest file for the current time.
 */
fileName pvtuFileName(const Time& runTime)
{
    return outputDir(runTime)/"dgToVTKLagrange.pvtu";
}


/**
 * \brief Return the top-level PVD timeline manifest for all exported times.
 */
fileName pvdFileName(const Time& runTime)
{
    return runTime.globalPath()/"VTK_Lagrange"/"dgToVTKLagrange.pvd";
}


/**
 * \brief Create the output directory for the current time if needed.
 */
void ensureOutputDir(const Time& runTime)
{
    mkDir(outputDir(runTime));
}


/**
 * \brief Write a ParaView PVD manifest that references all exported times.
 */
void writePvdManifest(const Time& runTime, const instantList& timeDirs)
{
    const fileName manifestPath = pvdFileName(runTime);
    mkDir(manifestPath.path());

    OFstream os(manifestPath);
    os.precision(16);

    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <Collection>" << nl;

    forAll(timeDirs, timeI)
    {
        const instant& timeInst = timeDirs[timeI];
        const fileName dataSet =
            fileName(timeInst.name())
           /(Pstream::parRun()
                ? fileName("dgToVTKLagrange.pvtu")
                : fileName("dgToVTKLagrange.vtu"));

        os  << "    <DataSet timestep=\"" << timeInst.value()
            << "\" group=\"\" part=\"0\" file=\"" << dataSet.c_str()
            << "\"/>" << nl;
    }

    os  << "  </Collection>" << nl
        << "</VTKFile>" << nl;
}

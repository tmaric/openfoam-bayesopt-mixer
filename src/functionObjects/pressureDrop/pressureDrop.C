#include "pressureDrop.H"

#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "surfaceFields.H"
#include "volFields.H"

#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pressureDrop, 0);
    addToRunTimeSelectionTable(functionObject, pressureDrop, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pressureDrop::pressureDrop
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pressureFieldName_("p"),
    patch1Name_("inlet"),
    patch2Name_("outlet"),
    csvFileName_("pressureDrop.csv"),
    csvHeaderWritten_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pressureDrop::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("field", pressureFieldName_);
    dict.lookup("patch1") >> patch1Name_;
    dict.lookup("patch2") >> patch2Name_;
    dict.readIfPresent("csvFile", csvFileName_);

    csvHeaderWritten_ = false;

    return true;
}


Foam::label Foam::functionObjects::pressureDrop::patchIndex
(
    const word& patchName
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label patchi = patches.findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch '" << patchName << "' not found." << nl
            << "Available patches: " << patches.names() << nl
            << exit(FatalError);
    }

    return patchi;
}


bool Foam::functionObjects::pressureDrop::writeCsvHeader()
{
    if (!Pstream::master())
    {
        return true;
    }

    const fileName csvPath = time_.globalPath()/csvFileName_;

    if (isFile(csvPath))
    {
        csvHeaderWritten_ = true;
        return true;
    }

    mkDir(csvPath.path());

    std::ofstream os(csvPath.c_str());

    if (!os.good())
    {
        WarningInFunction
            << "Cannot create CSV file at " << csvPath << nl;
        return false;
    }

    os
        << "time"
        << ",field"
        << ",patch1"
        << ",patch2"
        << ",patch1_average_m2_s2"
        << ",patch2_average_m2_s2"
        << ",pressure_drop_m2_s2"
        << '\n';

    csvHeaderWritten_ = true;
    return true;
}


bool Foam::functionObjects::pressureDrop::execute()
{
    return true;
}


bool Foam::functionObjects::pressureDrop::write()
{
    if (!csvHeaderWritten_ && !writeCsvHeader())
    {
        return false;
    }

    const label patch1i = patchIndex(patch1Name_);
    const label patch2i = patchIndex(patch2Name_);

    const volScalarField& p = mesh_.lookupObject<volScalarField>(pressureFieldName_);

    const fvPatchScalarField& pPatch1 = p.boundaryField()[patch1i];
    const fvPatchScalarField& pPatch2 = p.boundaryField()[patch2i];

    const fvsPatchScalarField& magSf1 = mesh_.magSf().boundaryField()[patch1i];
    const fvsPatchScalarField& magSf2 = mesh_.magSf().boundaryField()[patch2i];

    scalar area1 = 0.0;
    scalar areaWeightedP1 = 0.0;

    forAll(pPatch1, facei)
    {
        area1 += magSf1[facei];
        areaWeightedP1 += magSf1[facei]*pPatch1[facei];
    }

    scalar area2 = 0.0;
    scalar areaWeightedP2 = 0.0;

    forAll(pPatch2, facei)
    {
        area2 += magSf2[facei];
        areaWeightedP2 += magSf2[facei]*pPatch2[facei];
    }

    reduce(area1, sumOp<scalar>());
    reduce(areaWeightedP1, sumOp<scalar>());
    reduce(area2, sumOp<scalar>());
    reduce(areaWeightedP2, sumOp<scalar>());

    if (area1 <= SMALL || area2 <= SMALL)
    {
        WarningInFunction
            << "Patch areas are too small for averaging: "
            << patch1Name_ << " area=" << area1 << ", "
            << patch2Name_ << " area=" << area2 << nl;
        return true;
    }

    const scalar p1Average = areaWeightedP1/area1;
    const scalar p2Average = areaWeightedP2/area2;
    const scalar deltaP = p1Average - p2Average;

    Info<< type() << " " << name() << ": "
        << "avg(" << pressureFieldName_ << "," << patch1Name_ << ")=" << p1Average
        << " m2/s2, avg(" << pressureFieldName_ << "," << patch2Name_ << ")=" << p2Average
        << " m2/s2, deltaP=" << deltaP << " m2/s2"
        << nl;

    if (Pstream::master())
    {
        const fileName csvPath = time_.globalPath()/csvFileName_;
        std::ofstream os(csvPath.c_str(), std::ios::app);

        if (!os.good())
        {
            WarningInFunction
                << "Cannot append to CSV file at " << csvPath << nl;
            return false;
        }

        os
            << time_.value()
            << ',' << pressureFieldName_
            << ',' << patch1Name_
            << ',' << patch2Name_
            << ',' << p1Average
            << ',' << p2Average
            << ',' << deltaP
            << '\n';
    }

    return true;
}


// ************************************************************************* //

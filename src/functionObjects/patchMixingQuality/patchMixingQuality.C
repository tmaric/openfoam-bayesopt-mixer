#include "patchMixingQuality.H"

#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"
#include "OFstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "volFields.H"

#include <fstream>
#include <cmath>
#include <limits>

namespace Foam
{
namespace functionObjects
{

defineTypeNameAndDebug(patchMixingQuality, 0);
addToRunTimeSelectionTable(functionObject, patchMixingQuality, dictionary);


patchMixingQuality::patchMixingQuality
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    scalarFieldName_("T"),
    patchName_("outlet"),
    weighting_("none"),
    meanMode_("fromField"),
    aMean_(0.5),
    epsilon_(SMALL),
    reportMixingIndex_(true),
    reportAbsoluteDeviations_(false),
    csvFileName_("mixing.csv"),
    csvHeaderWritten_(false)
{
    read(dict);
}


bool patchMixingQuality::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("scalarField") >> scalarFieldName_;
    dict.lookup("patch") >> patchName_;

    dict.readIfPresent("weighting", weighting_);
    dict.readIfPresent("meanMode", meanMode_);
    dict.readIfPresent("aMean", aMean_);
    dict.readIfPresent("epsilon", epsilon_);
    dict.readIfPresent("reportMixingIndex", reportMixingIndex_);
    dict.readIfPresent("reportAbsoluteDeviations", reportAbsoluteDeviations_);
    dict.readIfPresent("csvFile", csvFileName_);

    if (weighting_ != "none" && weighting_ != "phi" && weighting_ != "Un")
    {
        FatalIOErrorInFunction(dict)
            << "Invalid weighting '" << weighting_
            << "'. Valid values: none | phi | Un" << nl
            << exit(FatalIOError);
    }

    if (meanMode_ != "fromField" && meanMode_ != "fromInletRatio")
    {
        FatalIOErrorInFunction(dict)
            << "Invalid meanMode '" << meanMode_
            << "'. Valid values: fromField | fromInletRatio" << nl
            << exit(FatalIOError);
    }

    epsilon_ = max(epsilon_, SMALL);

    csvHeaderWritten_ = false;

    return true;
}


label patchMixingQuality::patchIndex() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label patchi = patches.findPatchID(patchName_);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch '" << patchName_ << "' not found." << nl
            << "Available patches: " << patches.names() << nl
            << exit(FatalError);
    }

    return patchi;
}


scalar patchMixingQuality::safeDiv
(
    const scalar numerator,
    const scalar denominator,
    const scalar eps
)
{
    return numerator/(mag(denominator) > eps ? denominator : eps);
}


bool patchMixingQuality::writeCsvHeader()
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
        << ",scalar_field"
        << ",patch"
        << ",weighting"
        << ",mean_mode"
        << ",mean_concentration"
        << ",standard_deviation"
        << ",coefficient_of_variation"
        << ",mixing_coefficient"
        << ",flux_weighted_mean_concentration"
        << ",flux_weighted_standard_deviation"
        << ",flux_weighted_coefficient_of_variation"
        << ",flux_weighted_mixing_coefficient"
        << ",intensity_of_segregation"
        << ",flux_weighted_intensity_of_segregation"
        << ",relative_standard_deviation"
        << ",flux_weighted_relative_standard_deviation"
        << ",mixing_index_rsd"
        << ",mixing_index_intensity"
        << ",max_absolute_relative_deviation"
        << ",mean_absolute_relative_deviation"
        << ",delta_x_min"
        << ",delta_x_max"
        << ",delta_x_mean"
        << ",delta_x_average"
        << '\n';

    csvHeaderWritten_ = true;
    return true;
}


bool patchMixingQuality::execute()
{
    return true;
}


bool patchMixingQuality::write()
{
    const scalar nanValue = std::numeric_limits<scalar>::quiet_NaN();

    const label patchi = patchIndex();

    const volScalarField& a = mesh_.lookupObject<volScalarField>(scalarFieldName_);
    const fvPatchScalarField& aPatch = a.boundaryField()[patchi];

    const label nFacesLocal = aPatch.size();
    label nFacesGlobal = nFacesLocal;
    reduce(nFacesGlobal, sumOp<label>());

    if (nFacesGlobal == 0)
    {
        WarningInFunction
            << "Patch '" << patchName_ << "' has no faces." << nl;
        return true;
    }

    const vectorField& Sf = mesh_.Sf().boundaryField()[patchi];

    scalar sumArea = 0.0;
    scalar sumA = 0.0;
    scalar sumA2 = 0.0;
    scalar aMin = GREAT;
    scalar aMax = -GREAT;

    scalar sumW = 0.0;
    scalar sumWA = 0.0;
    scalar sumWA2 = 0.0;

    scalar sumAbsDev = 0.0;

    const surfaceScalarField* phiPtr = nullptr;
    const volVectorField* UPtr = nullptr;

    if (weighting_ == "phi")
    {
        phiPtr = &mesh_.lookupObject<surfaceScalarField>("phi");
    }
    else if (weighting_ == "Un")
    {
        UPtr = &mesh_.lookupObject<volVectorField>("U");
    }

    forAll(aPatch, facei)
    {
        const scalar ai = aPatch[facei];
        const scalar area = mag(Sf[facei]);

        sumArea += area;
        sumA += area*ai;
        sumA2 += area*sqr(ai);

        aMin = min(aMin, ai);
        aMax = max(aMax, ai);

        scalar wi = 1.0;

        if (weighting_ == "phi")
        {
            const scalar phiFace = phiPtr->boundaryField()[patchi][facei];
            wi = max(phiFace, scalar(0));
        }
        else if (weighting_ == "Un")
        {
            const vector Ui = UPtr->boundaryField()[patchi][facei];
            const scalar UnA = Ui & Sf[facei];
            wi = max(UnA, scalar(0));
        }

        sumW += wi;
        sumWA += wi*ai;
        sumWA2 += wi*sqr(ai);
    }

    reduce(sumArea, sumOp<scalar>());
    reduce(sumA, sumOp<scalar>());
    reduce(sumA2, sumOp<scalar>());
    reduce(aMin, minOp<scalar>());
    reduce(aMax, maxOp<scalar>());

    reduce(sumW, sumOp<scalar>());
    reduce(sumWA, sumOp<scalar>());
    reduce(sumWA2, sumOp<scalar>());

    if (sumArea <= VSMALL)
    {
        WarningInFunction
            << "Patch '" << patchName_ << "' has negligible area." << nl;
        return true;
    }

    const scalar meanA = sumA/sumArea;
    const scalar varA = max(sumA2/sumArea - sqr(meanA), scalar(0));
    const scalar sigmaA = sqrt(varA);

    const scalar covRefMean = (meanMode_ == "fromInletRatio") ? aMean_ : meanA;
    const scalar covA = safeDiv(sigmaA, covRefMean, epsilon_);

    scalar meanAF = nanValue;
    scalar varAF = nanValue;
    scalar sigmaAF = nanValue;
    scalar covAF = nanValue;

    if (sumW > epsilon_)
    {
        meanAF = sumWA/sumW;
        varAF = max(sumWA2/sumW - sqr(meanAF), scalar(0));
        sigmaAF = sqrt(varAF);
        covAF = safeDiv(sigmaAF, meanAF, epsilon_);
    }

    const scalar sigma0sq = covRefMean*(1.0 - covRefMean);

    scalar intensity = nanValue;
    scalar intensityF = nanValue;
    scalar rsd = nanValue;
    scalar rsdF = nanValue;

    if (sigma0sq > epsilon_)
    {
        intensity = varA/sigma0sq;
        rsd = sigmaA/sqrt(sigma0sq);

        if (sumW > epsilon_)
        {
            intensityF = varAF/sigma0sq;
            rsdF = sigmaAF/sqrt(sigma0sq);
        }
    }

    scalar mixingIndexRsd = nanValue;
    scalar mixingIndexIntensity = nanValue;

    if (reportMixingIndex_)
    {
        if (!std::isnan(rsd))
        {
            mixingIndexRsd = 1.0 - rsd;
        }

        if (!std::isnan(intensity))
        {
            mixingIndexIntensity = 1.0 - intensity;
        }
    }

    scalar deltaMax = nanValue;
    scalar delta = nanValue;

    if (reportAbsoluteDeviations_)
    {
        forAll(aPatch, facei)
        {
            sumAbsDev += mag(Sf[facei])*mag(aPatch[facei] - covRefMean);
        }
        reduce(sumAbsDev, sumOp<scalar>());

        deltaMax = safeDiv
        (
            max(aMax - covRefMean, covRefMean - aMin),
            covRefMean,
            epsilon_
        );

        delta = safeDiv(sumAbsDev/sumArea, covRefMean, epsilon_);
    }

    const tmp<surfaceScalarField> tDeltaX = Foam::pow(mesh_.deltaCoeffs(), -1.0);
    const surfaceScalarField& deltaX = tDeltaX();

    const scalar deltaXMin = gMin(deltaX);
    const scalar deltaXMax = gMax(deltaX);
    const scalar deltaXMean = gAverage(deltaX);
    const scalar deltaXAverage = average(deltaX).value();

    Log
        << type() << ' ' << name() << " on patch=" << patchName_
        << ": mean=" << meanA
        << ", std_dev=" << sigmaA
        << ", mixing_coefficient=" << covA
        << ", intensity_of_segregation=" << intensity
        << ", delta_x_mean=" << deltaXMean
        << nl;

    if (!csvHeaderWritten_)
    {
        if (!writeCsvHeader())
        {
            return true;
        }
    }

    if (Pstream::master())
    {
        const fileName csvPath = time_.globalPath()/csvFileName_;
        std::ofstream os(csvPath.c_str(), std::ios::app);

        if (!os.good())
        {
            WarningInFunction
                << "Cannot append CSV file at " << csvPath << nl;
            return true;
        }

        os
            << time_.value()
            << ',' << scalarFieldName_
            << ',' << patchName_
            << ',' << weighting_
            << ',' << meanMode_
            << ',' << meanA
            << ',' << sigmaA
            << ',' << covA
            << ',' << covA
            << ',' << meanAF
            << ',' << sigmaAF
            << ',' << covAF
            << ',' << covAF
            << ',' << intensity
            << ',' << intensityF
            << ',' << rsd
            << ',' << rsdF
            << ',' << mixingIndexRsd
            << ',' << mixingIndexIntensity
            << ',' << deltaMax
            << ',' << delta
            << ',' << deltaXMin
            << ',' << deltaXMax
            << ',' << deltaXMean
            << ',' << deltaXAverage
            << '\n';
    }

    return true;
}


} // End namespace functionObjects
} // End namespace Foam

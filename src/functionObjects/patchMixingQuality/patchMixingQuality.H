#ifndef patchMixingQuality_H
#define patchMixingQuality_H

#include "fvMeshFunctionObject.H"
#include "word.H"

namespace Foam
{
namespace functionObjects
{

class patchMixingQuality
:
    public fvMeshFunctionObject
{
    // Private data

    word scalarFieldName_;
    word patchName_;
    word weighting_;
    word meanMode_;

    scalar aMean_;
    scalar epsilon_;

    bool reportMixingIndex_;
    bool reportAbsoluteDeviations_;

    fileName csvFileName_;
    bool csvHeaderWritten_;


    // Private Member Functions

    label patchIndex() const;

    bool writeCsvHeader();

    static scalar safeDiv
    (
        const scalar numerator,
        const scalar denominator,
        const scalar eps
    );


public:

    TypeName("patchMixingQuality");


    // Constructors

    patchMixingQuality
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );


    //- Destructor
    virtual ~patchMixingQuality() = default;


    // Member Functions

    virtual bool read(const dictionary& dict);

    virtual bool execute();

    virtual bool write();
};


} // End namespace functionObjects
} // End namespace Foam

#endif

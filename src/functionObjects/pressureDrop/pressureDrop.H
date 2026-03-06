#ifndef pressureDrop_H
#define pressureDrop_H

#include "fvMeshFunctionObject.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class pressureDrop Declaration
\*---------------------------------------------------------------------------*/

class pressureDrop
:
    public fvMeshFunctionObject
{
    // Private data

    word pressureFieldName_;
    word patch1Name_;
    word patch2Name_;

    fileName csvFileName_;
    bool csvHeaderWritten_;


    // Private member functions

    label patchIndex(const word& patchName) const;

    bool writeCsvHeader();


public:

    //- Runtime type information
    TypeName("pressureDrop");


    // Constructors

    pressureDrop
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    pressureDrop(const pressureDrop&) = delete;
    void operator=(const pressureDrop&) = delete;


    //- Destructor
    virtual ~pressureDrop() = default;


    // Member functions

    virtual bool read(const dictionary& dict);

    virtual bool execute();

    virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

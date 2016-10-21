/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "improvedNoSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::improvedNoSlipFvPatchVectorField::improvedNoSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::improvedNoSlipFvPatchVectorField::improvedNoSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::improvedNoSlipFvPatchVectorField::improvedNoSlipFvPatchVectorField
(
    const improvedNoSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::improvedNoSlipFvPatchVectorField::improvedNoSlipFvPatchVectorField
(
    const improvedNoSlipFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::improvedNoSlipFvPatchVectorField::improvedNoSlipFvPatchVectorField
(
    const improvedNoSlipFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::Field<Foam::vector> > Foam::improvedNoSlipFvPatchVectorField::gradientInternalCoeffs() const
{
    vectorField n = this->patch().nf();
    vectorField impCoeffs(this->size(), pTraits<vector>::zero);

    forAll(impCoeffs, faceI)
    {
	impCoeffs[faceI][0] = scalar(1.0) - n[faceI][0] * n[faceI][0];
	impCoeffs[faceI][1] = scalar(1.0) - n[faceI][1] * n[faceI][1];
	impCoeffs[faceI][2] = scalar(1.0) - n[faceI][2] * n[faceI][2];
    }

    return -impCoeffs*this->patch().deltaCoeffs();
}

Foam::tmp<Foam::Field<Foam::vector> > Foam::improvedNoSlipFvPatchVectorField::gradientBoundaryCoeffs() const
{
    vectorField n = this->patch().nf();
    vectorField expCoeffs(this->size(), pTraits<vector>::zero);
    vectorField vt = *this - ((*this & n)) * n;
    vectorField iF = this->patchInternalField();

    forAll(expCoeffs, faceI)
    {
	expCoeffs[faceI][0] = vt[faceI][0] + 
	                      n[faceI][0]*n[faceI][1]*iF[faceI][1] +
			      n[faceI][0]*n[faceI][2]*iF[faceI][2];
	expCoeffs[faceI][1] = vt[faceI][1] + 
	                      n[faceI][0]*n[faceI][1]*iF[faceI][0] +
			      n[faceI][1]*n[faceI][2]*iF[faceI][2];
	expCoeffs[faceI][2] = vt[faceI][2] + 
	                      n[faceI][0]*n[faceI][2]*iF[faceI][0] +
			      n[faceI][1]*n[faceI][2]*iF[faceI][1];
    }

    return expCoeffs*this->patch().deltaCoeffs();
}



void Foam::improvedNoSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        improvedNoSlipFvPatchVectorField
    );
}

// ************************************************************************* //

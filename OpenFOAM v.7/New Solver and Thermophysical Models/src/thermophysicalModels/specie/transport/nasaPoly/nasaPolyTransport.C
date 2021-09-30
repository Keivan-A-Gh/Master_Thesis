/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Keivan A.Gh / 28, Oct, 2020:

\*---------------------------------------------------------------------------*/

#include "nasaPolyTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::nasaPolyTransport<Thermo>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorInFunction
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::nasaPolyTransport<Thermo>::nasaPolyTransport(const dictionary& dict)
:
    Thermo(dict),
    Tlow_(readScalar(dict.subDict("transport").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("transport").lookup("Thigh"))),
    Tcommon_(readScalar(dict.subDict("transport").lookup("Tcommon"))),
    highMuKappaCoeffs_(dict.subDict("transport").lookup("highMuKappaCoeffs")),
    lowMuKappaCoeffs_(dict.subDict("transport").lookup("lowMuKappaCoeffs"))
{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::nasaPolyTransport<Thermo>::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;
	
    Thermo::write(os);

    dictionary dict("transport");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
    dict.add("Tcommon", Tcommon_);
    dict.add("highMuKappaCoeffs",highMuKappaCoeffs_);
    dict.add("lowMuKappaCoeffs",lowMuKappaCoeffs_);
    
    os  << indent << dict.dictName() << dict;
    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const nasaPolyTransport<Thermo>& npt
)
{
    npt.write(os);
    return os;
}


// ************************************************************************* //

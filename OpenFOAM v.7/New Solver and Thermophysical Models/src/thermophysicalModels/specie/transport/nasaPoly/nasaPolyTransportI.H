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
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline Foam::nasaPolyTransport<Thermo>::nasaPolyTransport
(
    const Thermo& t,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename nasaPolyTransport<Thermo>::coeffArray& highMuKappaCoeffs,
    const typename nasaPolyTransport<Thermo>::coeffArray& lowMuKappaCoeffs
)
:
    Thermo(t),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
	for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highMuKappaCoeffs_[coefLabel] = highMuKappaCoeffs[coefLabel];
            lowMuKappaCoeffs_[coefLabel] = lowMuKappaCoeffs[coefLabel];
        }
}


template<class Thermo>
inline const typename Foam::nasaPolyTransport<Thermo>::coeffArray&
Foam::nasaPolyTransport<Thermo>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowMuKappaCoeffs_;
    }
    else
    {
        return highMuKappaCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::nasaPolyTransport<Thermo>::nasaPolyTransport
(
    const word& name,
    const nasaPolyTransport& npt
)
:
    Thermo(name, npt),
    Tlow_(npt.Tlow_),
    Thigh_(npt.Thigh_),
    Tcommon_(npt.Tcommon_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highMuKappaCoeffs_[coefLabel] = npt.highMuKappaCoeffs_[coefLabel];
        lowMuKappaCoeffs_[coefLabel] = npt.lowMuKappaCoeffs_[coefLabel];
    }
}

template<class Thermo>
inline Foam::autoPtr<Foam::nasaPolyTransport<Thermo>>
Foam::nasaPolyTransport<Thermo>::clone() const
{
    return autoPtr<nasaPolyTransport<Thermo>>
    (
        new nasaPolyTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::nasaPolyTransport<Thermo>>
Foam::nasaPolyTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<nasaPolyTransport<Thermo>>
    (
        new nasaPolyTransport<Thermo>(dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "attempt to use nasaPolyTransport<Thermo>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::Tlow() const
{
    return Tlow_;
}


template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::Thigh() const
{
    return Thigh_;
}


template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::Tcommon() const
{
    return Tcommon_;
}


template<class Thermo>
inline const typename Foam::nasaPolyTransport<Thermo>::coeffArray&
Foam::nasaPolyTransport<Thermo>::highMuKappaCoeffs() const
{
    return highMuKappaCoeffs_;
}


template<class Thermo>
inline const typename Foam::nasaPolyTransport<Thermo>::coeffArray&
Foam::nasaPolyTransport<Thermo>::lowMuKappaCoeffs() const
{
    return lowMuKappaCoeffs_;
}


template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return (exp(a[0]*log(T) + (a[1] + a[2]/T)/T + a[3]))*1e-07;
}


template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return (exp(a[4]*log(T) + (a[5] + a[6]/T)/T + a[7]))*1e-04;
}


template<class Thermo>
inline Foam::scalar Foam::nasaPolyTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::nasaPolyTransport<Thermo>::operator+=
(
    const nasaPolyTransport<Thermo>& npt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(npt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = npt.Y()/this->Y();

        Tlow_ = max(Tlow_, npt.Tlow_);
        Thigh_ = min(Thigh_, npt.Thigh_);

        if
        (
            nasaPolyTransport<Thermo>::debug
         && notEqual(Tcommon_, npt.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << Tcommon_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << npt.Tcommon_ << " for "
                << (npt.name().size() ? npt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<nasaPolyTransport<Thermo>::nCoeffs_;
            coefLabel++
        )
        {
            highMuKappaCoeffs_[coefLabel] =
                Y1*highMuKappaCoeffs_[coefLabel]
              + Y2*npt.highMuKappaCoeffs_[coefLabel];

            lowMuKappaCoeffs_[coefLabel] =
                Y1*lowMuKappaCoeffs_[coefLabel]
              + Y2*npt.lowMuKappaCoeffs_[coefLabel];
        }
    }
}


template<class Thermo>
inline void Foam::nasaPolyTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::nasaPolyTransport<Thermo> Foam::operator+
(
    const nasaPolyTransport<Thermo>& npt1,
    const nasaPolyTransport<Thermo>& npt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(npt1)
      + static_cast<const Thermo&>(npt2)
    );

    if (mag(t.Y()) < small)
    {
        return nasaPolyTransport<Thermo>
        (
            t,
	    0,
            npt1.Tlow_,
            npt1.Thigh_,
            npt1.Tcommon_,
            npt1.highMuKappaCoeffs_,
            npt1.lowMuKappaCoeffs_
        );
    }
    else
    {
        scalar Y1 = npt1.Y()/t.Y();
        scalar Y2 = npt2.Y()/t.Y();

        typename nasaPolyTransport<Thermo>::coeffArray highMuKappaCoeffs;
        typename nasaPolyTransport<Thermo>::coeffArray lowMuKappaCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<nasaPolyTransport<Thermo>::nCoeffs_;
            coefLabel++
        )
        {
            highMuKappaCoeffs[coefLabel] =
                Y1*npt1.highMuKappaCoeffs_[coefLabel]
              + Y2*npt2.highMuKappaCoeffs_[coefLabel];

            lowMuKappaCoeffs[coefLabel] =
                Y1*npt1.lowMuKappaCoeffs_[coefLabel]
              + Y2*npt2.lowMuKappaCoeffs_[coefLabel];
        }

        if
        (
            nasaPolyTransport<Thermo>::debug
         && notEqual(npt1.Tcommon_, npt2.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << npt1.Tcommon_ << " for "
                << (npt1.name().size() ? npt1.name() : "others")
                << " != " << npt2.Tcommon_ << " for "
                << (npt2.name().size() ? npt2.name() : "others")
                << exit(FatalError);
        }

        return nasaPolyTransport<Thermo>
        (
            t,
            max(npt1.Tlow_, npt2.Tlow_),
            min(npt1.Thigh_, npt2.Thigh_),
            npt1.Tcommon_,
            highMuKappaCoeffs,
            lowMuKappaCoeffs
        );
    }
}


template<class Thermo>
inline Foam::nasaPolyTransport<Thermo> Foam::operator*
(
    const scalar s,
    const nasaPolyTransport<Thermo>& npt
)
{
    return nasaPolyTransport<Thermo>
    (
        s*static_cast<const Thermo&>(npt),
        npt.Tlow_,
        npt.Thigh_,
        npt.Tcommon_,
        npt.highMuKappaCoeffs_,
        npt.lowMuKappaCoeffs_
    );
}


// ************************************************************************* //

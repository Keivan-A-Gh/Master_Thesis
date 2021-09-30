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

#include "nasa9PolyThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::nasa9PolyThermo<EquationOfState>::nasa9PolyThermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon1,
    const scalar Tcommon2,
    const typename nasa9PolyThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename nasa9PolyThermo<EquationOfState>::coeffArray& midCpCoeffs,
    const typename nasa9PolyThermo<EquationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon1_(Tcommon1),
    Tcommon2_(Tcommon2)
{
    if (convertCoeffs)
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel]*this->R();
	    midCpCoeffs_[coefLabel] = midCpCoeffs[coefLabel]*this->R();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel]*this->R();
        }
    }
    else
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
	    midCpCoeffs_[coefLabel] = midCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}


template<class EquationOfState>
inline const typename Foam::nasa9PolyThermo<EquationOfState>::coeffArray&
Foam::nasa9PolyThermo<EquationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon1_)
    {
        return lowCpCoeffs_;
    }
    else if (T < Tcommon2_)
    {
	return midCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::nasa9PolyThermo<EquationOfState>::nasa9PolyThermo
(
    const word& name,
    const nasa9PolyThermo& npt
)
:
    EquationOfState(name, npt),
    Tlow_(npt.Tlow_),
    Thigh_(npt.Thigh_),
    Tcommon1_(npt.Tcommon1_),
    Tcommon2_(npt.Tcommon2_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = npt.highCpCoeffs_[coefLabel];
	midCpCoeffs_[coefLabel] = npt.midCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = npt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "attempt to use nasa9PolyThermo<EquationOfState>"
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


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Tcommon1() const
{
    return Tcommon1_;
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Tcommon2() const
{
    return Tcommon2_;
}


template<class EquationOfState>
inline const typename Foam::nasa9PolyThermo<EquationOfState>::coeffArray&
Foam::nasa9PolyThermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::nasa9PolyThermo<EquationOfState>::coeffArray&
Foam::nasa9PolyThermo<EquationOfState>::midCpCoeffs() const
{
    return midCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::nasa9PolyThermo<EquationOfState>::coeffArray&
Foam::nasa9PolyThermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        ((((a[6]*T + a[5])*T + a[4])*T + a[3])*T + a[2])
	  + (a[0]/T + a[1])/T
      + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        ((((a[6]/5.0*T + a[5]/4.0)*T + a[4]/3.0)*T + a[3]/2.0)*T + a[2])*T
	  + (-a[0]/T + a[1]*log(T))
	  + a[7]
    ) + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::Hc() const
{
    const coeffArray& a = lowCpCoeffs_;
    return
    (
        ((((a[6]/5.0*Tstd + a[5]/4.0)*Tstd + a[4]/3.0)*Tstd + a[3]/2.0)*Tstd + a[2])*Tstd 
	  + (-a[0]/Tstd + a[1]*log(Tstd))
	  + a[7]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        (((a[6]/4.0*T + a[5]/3.0)*T + a[4]/2.0)*T + a[3])*T + a[2]*log(T)
	+ (-a[0]/2.0/T - a[1])/T
        + a[8]
    ) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return -(
		((-a[0]/T + a[1]*log(T))/T + a[2] + a[7]/T)/T 
		+ a[3]/2 + T*(a[4]/3 + T*(a[5]/4 + T*a[6]/5))
	    );
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9PolyThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
	(((4*a[6]*T + 3*a[5])*T + 2*a[4])*T + a[3])
	+ ((-2*a[0]/T - a[1])/T)/T
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::nasa9PolyThermo<EquationOfState>::operator+=
(
    const nasa9PolyThermo<EquationOfState>& npt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(npt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = npt.Y()/this->Y();

        Tlow_ = max(Tlow_, npt.Tlow_);
        Thigh_ = min(Thigh_, npt.Thigh_);

        if
        (
            nasa9PolyThermo<EquationOfState>::debug
         && notEqual(Tcommon1_, npt.Tcommon1_)
        )
        {
            FatalErrorInFunction
                << "Tcommon1 " << Tcommon1_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << npt.Tcommon1_ << " for "
                << (npt.name().size() ? npt.name() : "others")
                << exit(FatalError);
        }
		
		if
        (
            nasa9PolyThermo<EquationOfState>::debug
         && notEqual(Tcommon2_, npt.Tcommon2_)
        )
        {
            FatalErrorInFunction
                << "Tcommon2 " << Tcommon2_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << npt.Tcommon2_ << " for "
                << (npt.name().size() ? npt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<nasa9PolyThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs_[coefLabel] =
                Y1*highCpCoeffs_[coefLabel]
              + Y2*npt.highCpCoeffs_[coefLabel];

	    midCpCoeffs_[coefLabel] =
                Y1*midCpCoeffs_[coefLabel]
              + Y2*npt.midCpCoeffs_[coefLabel];
			  
            lowCpCoeffs_[coefLabel] =
                Y1*lowCpCoeffs_[coefLabel]
              + Y2*npt.lowCpCoeffs_[coefLabel];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::nasa9PolyThermo<EquationOfState> Foam::operator+
(
    const nasa9PolyThermo<EquationOfState>& npt1,
    const nasa9PolyThermo<EquationOfState>& npt2
)
{
    EquationOfState eofs = npt1;
    eofs += npt2;

    if (mag(eofs.Y()) < small)
    {
        return nasa9PolyThermo<EquationOfState>
        (
            eofs,
            npt1.Tlow_,
            npt1.Thigh_,
            npt1.Tcommon1_,
	    npt1.Tcommon2_,
            npt1.highCpCoeffs_,
	    npt1.midCpCoeffs_,
            npt1.lowCpCoeffs_
        );
    }
    else
    {
        const scalar Y1 = npt1.Y()/eofs.Y();
        const scalar Y2 = npt2.Y()/eofs.Y();

        typename nasa9PolyThermo<EquationOfState>::coeffArray highCpCoeffs;
	typename nasa9PolyThermo<EquationOfState>::coeffArray midCpCoeffs;
        typename nasa9PolyThermo<EquationOfState>::coeffArray lowCpCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<nasa9PolyThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs[coefLabel] =
                Y1*npt1.highCpCoeffs_[coefLabel]
              + Y2*npt2.highCpCoeffs_[coefLabel];
			  
	    midCpCoeffs[coefLabel] =
                Y1*npt1.midCpCoeffs_[coefLabel]
              + Y2*npt2.midCpCoeffs_[coefLabel];

            lowCpCoeffs[coefLabel] =
                Y1*npt1.lowCpCoeffs_[coefLabel]
              + Y2*npt2.lowCpCoeffs_[coefLabel];
        }

        if
        (
            nasa9PolyThermo<EquationOfState>::debug
         && notEqual(npt1.Tcommon1_, npt2.Tcommon1_)
        )
        {
            FatalErrorInFunction
                << "Tcommon1 " << npt1.Tcommon1_ << " for "
                << (npt1.name().size() ? npt1.name() : "others")
                << " != " << npt2.Tcommon1_ << " for "
                << (npt2.name().size() ? npt2.name() : "others")
                << exit(FatalError);
        }
		
	if
        (
            nasa9PolyThermo<EquationOfState>::debug
         && notEqual(npt1.Tcommon2_, npt2.Tcommon2_)
        )
        {
            FatalErrorInFunction
                << "Tcommon2 " << npt1.Tcommon2_ << " for "
                << (npt1.name().size() ? npt1.name() : "others")
                << " != " << npt2.Tcommon2_ << " for "
                << (npt2.name().size() ? npt2.name() : "others")
                << exit(FatalError);
        }

        return nasa9PolyThermo<EquationOfState>
        (
            eofs,
            max(npt1.Tlow_, npt2.Tlow_),
            min(npt1.Thigh_, npt2.Thigh_),
            npt1.Tcommon1_,
	    npt1.Tcommon2_,
            highCpCoeffs,
	    midCpCoeffs,
            lowCpCoeffs
        );
    }
}


template<class EquationOfState>
inline Foam::nasa9PolyThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const nasa9PolyThermo<EquationOfState>& npt
)
{
    return nasa9PolyThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(npt),
        npt.Tlow_,
        npt.Thigh_,
        npt.Tcommon1_,
	npt.Tcommon2_,
        npt.highCpCoeffs_,
	npt.midCpCoeffs_,
        npt.lowCpCoeffs_
    );
}


template<class EquationOfState>
inline Foam::nasa9PolyThermo<EquationOfState> Foam::operator==
(
    const nasa9PolyThermo<EquationOfState>& npt1,
    const nasa9PolyThermo<EquationOfState>& npt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(npt1)
     == static_cast<const EquationOfState&>(npt2)
    );

    const scalar Y1 = npt2.Y()/eofs.Y();
    const scalar Y2 = npt1.Y()/eofs.Y();

    typename nasa9PolyThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename nasa9PolyThermo<EquationOfState>::coeffArray midCpCoeffs;
    typename nasa9PolyThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        label coefLabel=0;
        coefLabel<nasa9PolyThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            Y1*npt2.highCpCoeffs_[coefLabel]
          - Y2*npt1.highCpCoeffs_[coefLabel];
		  
	midCpCoeffs[coefLabel] =
            Y1*npt2.midCpCoeffs_[coefLabel]
          - Y2*npt1.midCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            Y1*npt2.lowCpCoeffs_[coefLabel]
          - Y2*npt1.lowCpCoeffs_[coefLabel];
    }

    if
    (
        nasa9PolyThermo<EquationOfState>::debug
     && notEqual(npt2.Tcommon1_, npt1.Tcommon1_)
    )
    {
        FatalErrorInFunction
            << "Tcommon1 " << npt2.Tcommon1_ << " for "
            << (npt2.name().size() ? npt2.name() : "others")
            << " != " << npt1.Tcommon1_ << " for "
            << (npt1.name().size() ? npt1.name() : "others")
            << exit(FatalError);
    }
	
    if
    (
        nasa9PolyThermo<EquationOfState>::debug
     && notEqual(npt2.Tcommon2_, npt1.Tcommon2_)
    )
    {
        FatalErrorInFunction
            << "Tcommon2 " << npt2.Tcommon2_ << " for "
            << (npt2.name().size() ? npt2.name() : "others")
            << " != " << npt1.Tcommon2_ << " for "
            << (npt1.name().size() ? npt1.name() : "others")
            << exit(FatalError);
    }

    return nasa9PolyThermo<EquationOfState>
    (
        eofs,
        max(npt2.Tlow_, npt1.Tlow_),
        min(npt2.Thigh_, npt1.Thigh_),
        npt2.Tcommon1_,
	npt2.Tcommon2_,
        highCpCoeffs,
	midCpCoeffs,
        lowCpCoeffs
    );
}


// ************************************************************************* //

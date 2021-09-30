/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
Class
    Foam::nasaPolyTransport

Description
    Transport package using NASA-7 transport polynomial for
    calculating \mu and \kappa.

Usage

    Property			| Description
    highMuKappaCoeffs<8>	| Dynamic viscosity + Thermal Conductivity polynomial coefficients from Tcommon to Thigh
    lowMuKappaCoeffs<8>		| Dynamic viscosity + Thermal Conductivity polynomial coefficients from Tlow to Tcommon


    Example of the specification of the transport properties:
    transport
    {
        highMuKappaCoeffs<8>	( 0.70504381, 0.36287686e+02, -0.72255550e+04, 0.41921607, highkappaCoeff, highkappaCoeff, highkappaCoeff, highkappaCoeff );
        lowMuKappaCoeffs<8>	( 0.68887644, 0.48727168e+01, -0.59565053e+03, 0.55569577, lowkappaCoeff, lowkappaCoeff, lowkappaCoeff, lowkappaCoeff );
    }

    The polynomial expressions are evaluated as so:
    
    if (T < Tcommon)
    {
       \mu    = exp(0.68887644*log(T) + 0.48727168e+01/T + -0.59565053e+03/T^2 + 0.55569577)*1e-7
    }
    else
    {
       \mu    = exp(0.70504381*log(T) + 0.36287686e+02/T + -0.72255550e+04/T^2 + 0.41921607)*1e-7
    }

    if (T < Tcommon)
    {
       \kappa    = exp(...)*1e-4
    }
    else
    {
       \kappa    = exp(...)*1e-4
    }

Note
    - Dynamic viscosity polynomial coefficients evaluate to an expression in
      [Pa.s].
    - Thermal conductivity polynomial coefficients evaluate to an expression in
      [W/m/K].

SourceFiles
    nasaPolyTransportI.H
    nasaPolyTransport.C

\*---------------------------------------------------------------------------*/

#ifndef nasaPolyTransport_H
#define nasaPolyTransport_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class nasaPolyTransport;

template<class Thermo>
inline nasaPolyTransport<Thermo> operator+
(
    const nasaPolyTransport<Thermo>&,
    const nasaPolyTransport<Thermo>&
);

template<class Thermo>
inline nasaPolyTransport<Thermo> operator*
(
    const scalar,
    const nasaPolyTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const nasaPolyTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class nasaPolyTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class nasaPolyTransport
:
    public Thermo
{

public:

    // Public data

        static const int nCoeffs_ = 8;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private Data

        // Temperature limits of applicability of functions [K]
        scalar Tlow_, Thigh_, Tcommon_;
		
	// Coefficient arrays consist of both dynamic viscosity and
	// thermal consuctivity coefficients for different temperature range
        coeffArray highMuKappaCoeffs_;
	coeffArray lowMuKappaCoeffs_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline nasaPolyTransport
        (
            const Thermo& t,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon,
            const coeffArray& highMuKappaCoeffs,
	    const coeffArray& lowMuKappaCoeffs
        );

        //- Construct as a named copy
        inline nasaPolyTransport(const word&, const nasaPolyTransport&);
		
	//- Construct from dictionary
        nasaPolyTransport(const dictionary& dict);

	//- Construct and return a clone
        inline autoPtr<nasaPolyTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<nasaPolyTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "nasaPoly<" + Thermo::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the common temperature
            inline scalar Tcommon() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& highMuKappaCoeffs() const;
			
	    //- Return const access to the low temperature poly coefficients
            inline const coeffArray& lowMuKappaCoeffs() const;


        // Fundamental properties

            //- Dynamic viscosity [kg/m/s]
            inline scalar mu(const scalar p, const scalar T) const;

            //- Thermal conductivity [W/m/K]
	    inline scalar kappa(const scalar p, const scalar T) const;

	    //- Thermal diffusivity of enthalpy [kg/m/s]
	    inline scalar alphah(const scalar p, const scalar T) const;

	    // Species diffusivity
	    // inline scalar D(const scalar p, const scalar T) const;

        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const nasaPolyTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend nasaPolyTransport operator+ <Thermo>
        (
            const nasaPolyTransport&,
            const nasaPolyTransport&
        );

        friend nasaPolyTransport operator* <Thermo>
        (
            const scalar,
            const nasaPolyTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const nasaPolyTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nasaPolyTransportI.H"

#ifdef NoRepository
#include "nasaPolyTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

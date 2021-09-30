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
    Foam::nasa9PolyThermo

Description
    JANAF based thermodynamic package (with 9 Coefficients) templated
    into the equation of state.

SourceFiles
    nasa9PolyThermoI.H
    nasa9PolyThermo.C

\*---------------------------------------------------------------------------*/

#ifndef nasa9PolyThermo_H
#define nasa9PolyThermo_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class nasa9PolyThermo;

template<class EquationOfState>
inline nasa9PolyThermo<EquationOfState> operator+
(
    const nasa9PolyThermo<EquationOfState>&,
    const nasa9PolyThermo<EquationOfState>&
);

template<class EquationOfState>
inline nasa9PolyThermo<EquationOfState> operator*
(
    const scalar,
    const nasa9PolyThermo<EquationOfState>&
);

template<class EquationOfState>
inline nasa9PolyThermo<EquationOfState> operator==
(
    const nasa9PolyThermo<EquationOfState>&,
    const nasa9PolyThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const nasa9PolyThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class nasa9PolyThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class nasa9PolyThermo
:
    public EquationOfState
{

public:

    // Public data

        static const int nCoeffs_ = 9;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private Data

        //- Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon1_, Tcommon2_;

        coeffArray highCpCoeffs_;
	coeffArray midCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline nasa9PolyThermo
        (
            const EquationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon1,
	    const scalar Tcommon2,
            const coeffArray& highCpCoeffs,
	    const coeffArray& midCpCoeffs,
            const coeffArray& lowCpCoeffs,
            const bool convertCoeffs = false
        );

        //- Construct from dictionary
        nasa9PolyThermo(const dictionary& dict);

        //- Construct as a named copy
        inline nasa9PolyThermo(const word&, const nasa9PolyThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "nasa9Poly<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;

    // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the common1 temperature
            inline scalar Tcommon1() const;
			
	    //- Return const access to the common2 temperature
            inline scalar Tcommon2() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& highCpCoeffs() const;
			
	    //- Return const access to the middle temperature poly coefficients
            inline const coeffArray& midCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& lowCpCoeffs() const;

    // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            #include "HtoEthermo.H"

    // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


    // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const nasa9PolyThermo&);


    // Friend operators

        friend nasa9PolyThermo operator+ <EquationOfState>
        (
            const nasa9PolyThermo&,
            const nasa9PolyThermo&
        );

        friend nasa9PolyThermo operator* <EquationOfState>
        (
            const scalar,
            const nasa9PolyThermo&
        );

        friend nasa9PolyThermo operator== <EquationOfState>
        (
            const nasa9PolyThermo&,
            const nasa9PolyThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const nasa9PolyThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nasa9PolyThermoI.H"

#ifdef NoRepository
#include "nasa9PolyThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::prescribedSLERPMotionSolver

Description
    Prescribed mesh motion solver for an fvMesh.

    Applies prescribed non-uniform deformations using SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    prescribedSLERPMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedSLERPMotionSolver_H
#define prescribedSLERPMotionSolver_H

#include "displacementMotionSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                     Class prescribedSLERPMotionSolver Declaration
    \*---------------------------------------------------------------------------*/

    class prescribedSLERPMotionSolver
        : public displacementMotionSolver
    {
        // Private data

        wordRes patches_;

        //- Moving patches
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Coefficient for computing the radius/height in the X-direction
        const scalar coeffX_;

        //- Coefficient for computing the radius/height in the Y-direction
        const scalar coeffY_;

        //- Coefficient for computing the radius/height in the Z-direction
        const scalar coeffZ_;

        //- Choice of harmonic or non-harmonic motion
        const word motionType_;

        //- Choice of dimensionless amplitude for harmonic motion (if not included, the default amplitude is set to 1)
        const scalar amplitude_;

        //- Choice of frequency in rad/s for harmonic motion (if not included, the default frequency is set to 1 rad/s)
        const scalar omega_;

        //- Choice of path where the motion data is stored (if not included, the default path is set to $FOAM_CASE)
        const fileName path_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Index to associate each point to a specific height
        pointScalarField p2heightIdx_;

        //- Fraction to calculate the motion of each point from the one specified in the motion data files
        pointScalarField p2heightFraction_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Motion data at the current saved time-step in timeArray
        mutable List<List<double>> motionData;

        //- Motion data at the next saved time-step in timeArray
        mutable List<List<double>> nextMotionData;

        //- Interpolated motion data in case of non-harmonic motion
        mutable List<List<double>> motionDataInterpolated;

        //- Saved time values for the non-harmonic motion
        mutable List<double> timeArray;

        //- Current time index for reading the motion data in case of non-harmonic motion
        mutable int curTimeIndexRead_;

        // Private Member Functions

        //- No copy construct
        prescribedSLERPMotionSolver(
            const prescribedSLERPMotionSolver &) = delete;

        //- No copy assignment
        void operator=(const prescribedSLERPMotionSolver &) = delete;

    public:
        //- Runtime type information
        TypeName("prescribedSLERPMeshMotion");

        // Constructors

        //- Construct from polyMesh and IOdictionary
        prescribedSLERPMotionSolver(
            const polyMesh &,
            const IOdictionary &dict);

        //- Destructor
        ~prescribedSLERPMotionSolver() = default;

        // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        // Read and return the motion state
        virtual void readMotionData(std::string &filename, List<List<double>> &motion) const;

        //- Return the height index and fraction to associate each point to a specific height
        virtual std::pair<int, double> p2HeightAssociation(Vector<double> &point, List<List<double>> &motion, double N_heights);

        //- Return the current time state to make sure the current motion state is up to date
        virtual double timeFractionCalculation(double timeValue) const;

        //- Return the transformation septernion
        virtual septernion computePointTransformation(int index, double fraction, double mult, List<List<double>> &motion) const;

        //- Solve for motion
        virtual void solve();

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

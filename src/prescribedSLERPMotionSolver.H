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

    Applies prescribed non-uniform deformations using SLERP interpolation
    of movement as a function of distance to the object surface. Two modes
    are supported for computing the height or radius of each mesh point:

    - axial: signed projection onto the body axis direction
        h = axisDirection . p
      Suitable for bodies whose surrounding mesh is arranged in planes 
      perpendicular to the body axis or is unstructured.
      The axisDirection vector is normalised to a unit vector.

    - radial: signed weighted radial distance from a Cartesian rotation
      axis
        h = sign(p_dom) * sqrt((cx*x)^2 + (cy*y)^2 + (cz*z)^2)
      Suitable for bodies whose surrounding mesh is arranged azimuthally.
      Set exactly two radialCoeffs components to 1 (perpendicular to the
      rotation axis) and one to 0 (along the rotation axis). Section
      positions in motionData must be stored as radial distances from the
      rotation axis.

SourceFiles
    prescribedSLERPMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedSLERPMotionSolver_H
#define prescribedSLERPMotionSolver_H

#include "displacementMotionSolver.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class prescribedSLERPMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class prescribedSLERPMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Moving patch names
        wordRes patches_;

        //- Moving patches
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of deforming region)
        const scalar do_;

        //- Height mode: axial (dot product) or radial (weighted magnitude)
        const word heightMode_;

        //- Unit direction vector along the body span (axial mode only).
        //  Normalised to unit length in the constructor.
        vector axisDirection_;

        //- Weighting coefficients (cx, cy, cz) for the radial distance
        //  computation (radial mode only). Set non-zero components for
        //  coordinates perpendicular to the rotation axis.
        vector radialCoeffs_;

        //- Choice of harmonic or non-harmonic motion
        const word motionType_;

        //- Dimensionless amplitude for harmonic motion (default: 1)
        const scalar amplitude_;

        //- Angular frequency in rad/s for harmonic motion (default: 1)
        const scalar omega_;

        //- Path where motion data files are stored (default: $FOAM_CASE)
        const fileName path_;

        //- Current interpolation scale (1 at patches, 0 at outerDistance)
        pointScalarField scale_;

        //- Index associating each point to a specific height/radius section
        pointScalarField p2heightIdx_;

        //- Fraction for interpolating point motion between sections
        pointScalarField p2heightFraction_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Motion data at the current saved time step in timeArray
        mutable List<List<scalar>> motionData;

        //- Motion data at the next saved time step in timeArray
        mutable List<List<scalar>> nextMotionData;

        //- Interpolated motion data for non-harmonic motion
        mutable List<List<scalar>> motionDataInterpolated;

        //- Saved time values for non-harmonic motion
        mutable List<scalar> timeArray;

        //- Current time index for reading motion data files
        mutable label curTimeIndexRead_;


    // Private Member Functions

        //- No copy construct
        prescribedSLERPMotionSolver
        (
            const prescribedSLERPMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const prescribedSLERPMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("prescribedSLERPMeshMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        prescribedSLERPMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~prescribedSLERPMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Read motion data from file into motion list
        virtual void readMotionData
        (
            const fileName& filename,
            List<List<scalar>>& motion
        ) const;

        //- Return height index and fraction associating each point
        //  to a specific height/radius section
        virtual std::pair<label, scalar> p2HeightAssociation
        (
            const point& pt,
            const List<List<scalar>>& motion,
            const label N_heights
        ) const;

        //- Return time fraction for interpolation between motion data files
        virtual scalar timeFractionCalculation(scalar timeValue) const;

        //- Return the transformation septernion for a given point
        virtual septernion computePointTransformation
        (
            const label index,
            const scalar fraction,
            const scalar mult,
            const List<List<scalar>>& motion
        ) const;

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

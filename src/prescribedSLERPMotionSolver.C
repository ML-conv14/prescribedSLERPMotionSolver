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

\*---------------------------------------------------------------------------*/

#include "prescribedSLERPMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"
#include "septernion.H"
#include "quaternion.H"
#include "interpolateSplineXY.H"
#include "unitConversion.H"
#include "Time.H"
#include "IFstream.H"
#include <iomanip>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(prescribedSLERPMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        prescribedSLERPMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prescribedSLERPMotionSolver::prescribedSLERPMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    heightMode_(coeffDict().get<word>("heightMode")),
    axisDirection_(Zero),
    radialCoeffs_(Zero),
    motionType_(coeffDict().get<word>("motionType")),
    amplitude_(coeffDict().getOrDefault("amplitude", scalar(1))),
    omega_(coeffDict().getOrDefault("omega", scalar(1))),
    path_
    (
        coeffDict().found("path")
            ? coeffDict().get<fileName>("path").expand()
            : fileName("$FOAM_CASE").expand()
    ),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    p2heightIdx_
    (
        IOobject
        (
            "heightIdxMatching",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    p2heightFraction_
    (
        IOobject
        (
            "heightFractionMatching",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    curTimeIndexRead_(0)
{
    // Validate heightMode and read the relevant axis/coefficient vector
    if (heightMode_ == "axial")
    {
        axisDirection_ = coeffDict().get<vector>("axisDirection");
        const scalar magAxis = mag(axisDirection_);
        if (magAxis < SMALL)
        {
            FatalErrorInFunction
                << "axisDirection has zero magnitude"
                << exit(FatalError);
        }
        axisDirection_ /= magAxis;
    }
    else if (heightMode_ == "radial")
    {
        radialCoeffs_ = coeffDict().get<vector>("radialCoeffs");
    }
    else
    {
        FatalErrorInFunction
            << "Unknown heightMode: " << heightMode_ << nl
            << "Valid options are: axial, radial"
            << exit(FatalError);
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const pointMesh& pMesh = pointMesh::New(mesh);

    pointPatchDist pDist(pMesh, patchSet_, points0());

    // Scaling: 1 up to di then linear down to 0 at do away from patches
    scale_.primitiveFieldRef() =
        min
        (
            max
            (
                (do_ - pDist.primitiveField()) / (do_ - di_),
                scalar(0)
            ),
            scalar(1)
        );

    // Convert the scale function to a cosine
    scale_.primitiveFieldRef() =
        min
        (
            max
            (
                scalar(0.5) - scalar(0.5)
                    * cos
                      (
                          scale_.primitiveField()
                        * Foam::constant::mathematical::pi
                      ),
                scalar(0)
            ),
            scalar(1)
        );

    pointConstraints::New(pMesh).constrain(scale_);
    scale_.write();

    // Set the current time index to 0
    curTimeIndexRead_ = 0;

    // Read the motion data and compute height index and fraction
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (motionType_ == "harmonic")
    {
        const fileName filename = path_ + "/motionData.dat";
        readMotionData(filename, motionData);

        const label N_heights = motionData.size();

        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                auto result = p2HeightAssociation
                (
                    points0()[pointi],
                    motionData,
                    N_heights
                );

                p2heightIdx_.primitiveFieldRef()[pointi] = result.first;
                p2heightFraction_.primitiveFieldRef()[pointi] = result.second;
            }
            else
            {
                p2heightIdx_.primitiveFieldRef()[pointi] = 0;
                p2heightFraction_.primitiveFieldRef()[pointi] = 0;
            }
        }

        pointConstraints::New(pMesh).constrain(p2heightIdx_);
        p2heightIdx_.write();
        pointConstraints::New(pMesh).constrain(p2heightFraction_);
        p2heightFraction_.write();
    }
    else if (motionType_ == "nonHarmonic")
    {
        const fileName time_filename = path_ + "/timeArray.dat";
        IFstream timeFile(time_filename);

        if (!timeFile.good())
        {
            FatalErrorInFunction
                << "Error opening file: " << time_filename
                << exit(FatalError);
        }

        label N_times;
        timeFile >> N_times;

        scalar tmp_timeStep;
        for (label i = 0; i < N_times; ++i)
        {
            timeFile >> tmp_timeStep;
            timeArray.append(tmp_timeStep);
        }

        // Load first two motion data files
        {
            std::ostringstream stream1;
            stream1 << std::fixed << std::setprecision(9)
                    << timeArray[0];
            const fileName filename1 =
                path_ + "/motionData_" + stream1.str() + ".dat";
            readMotionData(filename1, motionData);
        }
        {
            std::ostringstream stream2;
            stream2 << std::fixed << std::setprecision(9)
                    << timeArray[1];
            const fileName filename2 =
                path_ + "/motionData_" + stream2.str() + ".dat";
            readMotionData(filename2, nextMotionData);
        }

        const label N_heights = motionData.size();

        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                auto result = p2HeightAssociation
                (
                    points0()[pointi],
                    motionData,
                    N_heights
                );

                p2heightIdx_.primitiveFieldRef()[pointi] = result.first;
                p2heightFraction_.primitiveFieldRef()[pointi] = result.second;
            }
            else
            {
                p2heightIdx_.primitiveFieldRef()[pointi] = 0;
                p2heightFraction_.primitiveFieldRef()[pointi] = 0;
            }
        }

        pointConstraints::New(pMesh).constrain(p2heightIdx_);
        p2heightIdx_.write();
        pointConstraints::New(pMesh).constrain(p2heightFraction_);
        p2heightFraction_.write();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown motion type: " << motionType_ << nl
            << "Valid options are: harmonic and nonHarmonic"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prescribedSLERPMotionSolver::readMotionData
(
    const fileName& filename,
    List<List<scalar>>& motion
) const
{
    IFstream inputFile(filename);

    if (!inputFile.good())
    {
        FatalErrorInFunction
            << "Error opening file: " << filename
            << exit(FatalError);
    }

    label N_heights;
    inputFile >> N_heights;

    const label N_columns = 10;

    for (label heightIdx = 0; heightIdx < N_heights; ++heightIdx)
    {
        List<scalar> tmp_height(N_columns);
        for (label i = 0; i < N_columns; ++i)
        {
            inputFile >> tmp_height[i];
        }
        motion.append(tmp_height);
    }
}


std::pair<Foam::label, Foam::scalar>
Foam::prescribedSLERPMotionSolver::p2HeightAssociation
(
    const point& pt,
    const List<List<scalar>>& motion,
    const label N_heights
) const
{
    scalar height;

    if (heightMode_ == "axial")
    {
        // Axial mode: signed projection onto body axis direction
        // h = axisDirection . p  (axisDirection_ normalised in constructor)
        height = (axisDirection_ & pt);
    }
    else // radial
    {
        // Radial mode: signed weighted radial distance
        // h = sign(p_dom) * sqrt((cx*x)^2 + (cy*y)^2 + (cz*z)^2)
        // where p_dom is the weighted coordinate with largest absolute
        // value, used to recover the sign of the radial distance.
        const vector weighted
        (
            radialCoeffs_.x() * pt.x(),
            radialCoeffs_.y() * pt.y(),
            radialCoeffs_.z() * pt.z()
        );

        const scalar radius = Foam::sqrt
        (
            max
            (
                sqr(weighted.x())
              + sqr(weighted.y())
              + sqr(weighted.z()),
                scalar(0)
            )
        );

        // Recover sign from dominant weighted component
        label dominantIndex = 0;
        scalar maxVal = Foam::mag(weighted.x());

        if (Foam::mag(weighted.y()) > maxVal)
        {
            maxVal = Foam::mag(weighted.y());
            dominantIndex = 1;
        }
        if (Foam::mag(weighted.z()) > maxVal)
        {
            dominantIndex = 2;
        }

        const scalar signVal =
            (weighted[dominantIndex] >= 0) ? scalar(1) : scalar(-1);

        height = signVal * radius;
    }

    // Assign point to height interval
    label index = 0;
    scalar fraction = scalar(0);

    if (height < motion[0][0])
    {
        index = 0;
        fraction = scalar(0);
    }
    else if (height > motion[N_heights - 1][0])
    {
        index = N_heights - 2;
        fraction = scalar(1);
    }
    else
    {
        for (label i = 0; i < N_heights - 1; i++)
        {
            if (height >= motion[i][0] && height <= motion[i + 1][0])
            {
                index = i;
                fraction =
                    (height - motion[i][0])
                  / (motion[i + 1][0] - motion[i][0]);
                break;
            }
        }
    }

    return std::pair<label, scalar>(index, fraction);
}


Foam::scalar
Foam::prescribedSLERPMotionSolver::timeFractionCalculation
(
    scalar timeValue
) const
{
    if (timeValue <= timeArray[curTimeIndexRead_ + 1])
    {
        const scalar t1 = timeArray[curTimeIndexRead_];
        const scalar t2 = timeArray[curTimeIndexRead_ + 1];
        return (timeValue - t1) / (t2 - t1);
    }
    else
    {
        curTimeIndexRead_++;

        if (curTimeIndexRead_ >= timeArray.size() - 1)
        {
            curTimeIndexRead_ = timeArray.size() - 2;
        }

        const scalar t1 = timeArray[curTimeIndexRead_];
        const scalar t2 = timeArray[curTimeIndexRead_ + 1];
        const scalar timeFraction = (timeValue - t1) / (t2 - t1);

        motionData = nextMotionData;
        nextMotionData.clear();

        std::ostringstream stream;
        stream << std::fixed << std::setprecision(9)
               << timeArray[curTimeIndexRead_ + 1];
        const fileName filename =
            path_ + "/motionData_" + stream.str() + ".dat";

        readMotionData(filename, nextMotionData);

        return timeFraction;
    }
}


Foam::septernion
Foam::prescribedSLERPMotionSolver::computePointTransformation
(
    const label index,
    const scalar fraction,
    const scalar mult,
    const List<List<scalar>>& motion
) const
{
    // Interpolate and scale translation
    const vector translation1
    (
        motion[index][4], motion[index][5], motion[index][6]
    );
    const vector translation2
    (
        motion[index + 1][4], motion[index + 1][5], motion[index + 1][6]
    );
    const vector interpolatedTranslation =
        ((scalar(1) - fraction) * translation1
      + fraction * translation2) * mult;

    // Extract and interpolate rotation (Euler angles)
    const vector rotation1
    (
        motion[index][7], motion[index][8], motion[index][9]
    );
    const vector rotation2
    (
        motion[index + 1][7], motion[index + 1][8], motion[index + 1][9]
    );
    const vector interpolatedRotation =
        ((scalar(1) - fraction) * rotation1
      + fraction * rotation2) * mult;

    // Convert rotation to quaternion
    const quaternion rotationQuat
    (
        quaternion::XYZ, interpolatedRotation * degToRad()
    );

    // Extract and interpolate centre of rotation
    const vector centre1
    (
        motion[index][1], motion[index][2], motion[index][3]
    );
    const vector centre2
    (
        motion[index + 1][1], motion[index + 1][2], motion[index + 1][3]
    );
    const vector interpolatedCentre =
        (scalar(1) - fraction) * centre1 + fraction * centre2;

    // Construct septernion for combined translation and rotation
    return septernion
    (
        septernion(-interpolatedCentre + -interpolatedTranslation)
      * rotationQuat
      * septernion(interpolatedCentre)
    );
}


Foam::tmp<Foam::pointField>
Foam::prescribedSLERPMotionSolver::curPoints() const
{
    const Time& t = this->db().time();

    if (motionType_ == "harmonic")
    {
        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                const label index = p2heightIdx_[pointi];
                const scalar fraction = p2heightFraction_[pointi];

                // Motion multiplier for harmonic motion
                const scalar m = amplitude_ * sin(omega_ * t.value());

                const septernion transformation = computePointTransformation
                (
                    index, fraction, m, motionData
                );

                // Use solid-body motion where scale = 1
                if (scale_[pointi] > scalar(1) - SMALL)
                {
                    pointDisplacement_.primitiveFieldRef()[pointi] =
                        transformation.transformPoint(points0_[pointi])
                      - points0_[pointi];
                }
                // Slerp septernion interpolation: scaled septernion
                else
                {
                    septernion ss
                    (
                        slerp(septernion::I, transformation, scale_[pointi])
                    );
                    pointDisplacement_.primitiveFieldRef()[pointi] =
                        ss.transformPoint(points0_[pointi])
                      - points0_[pointi];
                }
            }
            else
            {
                pointDisplacement_.primitiveFieldRef()[pointi] = vector::zero;
            }
        }

        return tmp<pointField>
        (
            new pointField(points0() + pointDisplacement_.primitiveField())
        );
    }
    else if (motionType_ == "nonHarmonic")
    {
        const scalar timeFraction = timeFractionCalculation(t.value());

        motionDataInterpolated = motionData;

        // Interpolate motion data at the current time step from the
        // closest two time steps saved in timeArray
        forAll(motionDataInterpolated, i)
        {
            motionDataInterpolated[i] =
                (scalar(1) - timeFraction) * motionData[i]
              + timeFraction * nextMotionData[i];
        }

        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                const label index = p2heightIdx_[pointi];
                const scalar fraction = p2heightFraction_[pointi];

                // For non-harmonic motion, the multiplier is always 1
                const scalar m = scalar(1);

                const septernion transformation = computePointTransformation
                (
                    index, fraction, m, motionDataInterpolated
                );

                // Use solid-body motion where scale = 1
                if (scale_[pointi] > scalar(1) - SMALL)
                {
                    pointDisplacement_.primitiveFieldRef()[pointi] =
                        transformation.transformPoint(points0_[pointi])
                      - points0_[pointi];
                }
                // Slerp septernion interpolation: scaled septernion
                else
                {
                    septernion ss
                    (
                        slerp(septernion::I, transformation, scale_[pointi])
                    );
                    pointDisplacement_.primitiveFieldRef()[pointi] =
                        ss.transformPoint(points0_[pointi])
                      - points0_[pointi];
                }
            }
            else
            {
                pointDisplacement_.primitiveFieldRef()[pointi] = vector::zero;
            }
        }

        return tmp<pointField>
        (
            new pointField(points0() + pointDisplacement_.primitiveField())
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unknown motion type: " << motionType_
            << exit(FatalError);

        return tmp<pointField>(nullptr);
    }
}


void Foam::prescribedSLERPMotionSolver::solve()
{
    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << nl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


// ************************************************************************* //

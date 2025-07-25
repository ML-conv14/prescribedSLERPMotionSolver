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
#include <iomanip>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(prescribedSLERPMotionSolver, 0);

    addToRunTimeSelectionTable(
        motionSolver,
        prescribedSLERPMotionSolver,
        dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prescribedSLERPMotionSolver::prescribedSLERPMotionSolver(
    const polyMesh &mesh,
    const IOdictionary &dict)
    : displacementMotionSolver(mesh, dict, typeName),
      patches_(coeffDict().get<wordRes>("patches")),
      patchSet_(mesh.boundaryMesh().patchSet(patches_)),
      di_(coeffDict().get<scalar>("innerDistance")),
      do_(coeffDict().get<scalar>("outerDistance")),
      coeffX_(coeffDict().get<scalar>("coefficientX")),
      coeffY_(coeffDict().get<scalar>("coefficientY")),
      coeffZ_(coeffDict().get<scalar>("coefficientZ")),
      motionType_(coeffDict().get<word>("motionType")),
      amplitude_(coeffDict().getOrDefault("amplitude", 1.0)),
      omega_(coeffDict().getOrDefault("omega", 1.0)),
      path_(
          coeffDict().found("path")
              ? coeffDict().get<fileName>("path").expand()
              : fileName("$FOAM_CASE").expand()
           ),

      scale_(
          IOobject(
              "motionScale",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false),
          pointMesh::New(mesh),
          dimensionedScalar(dimless, Zero)
          ),

      p2heightIdx_(
          IOobject(
              "heightIdxMatching",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false),
          pointMesh::New(mesh),
          dimensionedScalar(dimless, Zero)
          ),

      p2heightFraction_(
          IOobject(
              "heightFractionMatching",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false),
          pointMesh::New(mesh),
          dimensionedScalar(dimless, Zero)
          ),

      curTimeIndex_(-1)

{

// Calculate scaling factor everywhere
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const pointMesh &pMesh = pointMesh::New(mesh);

    pointPatchDist pDist(pMesh, patchSet_, points0());

    // Scaling: 1 up to di then linear down to 0 at do away from patches
    scale_.primitiveFieldRef() =
        min(
            max(
                (do_ - pDist.primitiveField()) / (do_ - di_),
                scalar(0)
                ),
            scalar(1)
            );

    // Convert the scale function to a cosine
    scale_.primitiveFieldRef() =
        min(
            max(
                0.5 - 0.5 * cos(scale_.primitiveField() * Foam::constant::mathematical::pi),
                scalar(0)
                ),
            scalar(1)
            );

    pointConstraints::New(pMesh).constrain(scale_);
    scale_.write();
    
    // Set the current time index to 0
    curTimeIndexRead_ = 0;


// Read the motion data from the file and compute height index and fraction
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (motionType_ == "harmonic")
    {
        // Read the motion data from the file
        std::string filename = path_ + "/motionData.dat";
        readMotionData(filename, motionData);

        double N_heights = motionData.size();

        // Calculate the height index and fraction for each point
        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                // Calculate the height index and fraction for the point
                auto result = p2HeightAssociation(points0()[pointi], motionData, N_heights);

                int index = result.first;
                double fraction = result.second;

                // Assign the height index and fraction to the point
                p2heightIdx_.primitiveFieldRef()[pointi] = index;
                p2heightFraction_.primitiveFieldRef()[pointi] = fraction;
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
        // Read time values from file
        std::string time_filename = path_ + "/timeArray.dat";
        std::ifstream timeFile(time_filename);
        if (!timeFile.is_open())
        {
            FatalErrorInFunction << "Error opening file: " << time_filename << exit(FatalError);
            return;
        }

        // Read number of time values
        double N_times;
        timeFile >> N_times;

        double tmp_timeStep;
        for (int i = 0; i < N_times; ++i)
        {
            timeFile >> tmp_timeStep;
            timeArray.append(tmp_timeStep);
        }

        timeFile.close();

        double timeDouble = timeArray[0];

        // Use a stringstream to format the time value with fixed precision
        std::ostringstream stream1;
        stream1 << std::fixed << std::setprecision(9) << timeDouble; // Set precision to 9 decimal places
        std::string filename1 = path_ + "/motionData_" + stream1.str() + ".dat";

        readMotionData(filename1, motionData);

        timeDouble = timeArray[1];

        // Use a stringstream to format the time value with fixed precision
        std::ostringstream stream2;
        stream2 << std::fixed << std::setprecision(9) << timeDouble; // Set precision to 9 decimal places
        std::string filename2 = path_ + "/motionData_" + stream2.str() + ".dat";

        readMotionData(filename2, nextMotionData);

        double N_heights = motionData.size();

        // Calculate the height index and fraction for each point
        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                // Calculate the height index and fraction for the point
                auto result = p2HeightAssociation(points0()[pointi], motionData, N_heights);

                int index = result.first;
                double fraction = result.second;

                // Assign the height index and fraction to the point
                p2heightIdx_.primitiveFieldRef()[pointi] = index;
                p2heightFraction_.primitiveFieldRef()[pointi] = fraction;
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
            << "Unknown motion type: " << motionType_ << endl
            << "Valid options are: harmonic and nonHarmonic" << exit(FatalError);
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Function to read the motion data from a file
void Foam::prescribedSLERPMotionSolver::readMotionData(std::string &filename, List<List<double>> &motion) const
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        FatalErrorInFunction << "Error opening file: " << filename << exit(FatalError);
    }

    // Read number of heights in the mesh motion file
    double N_heights;
    inputFile >> N_heights;

    int N_columns = 10; // Number of columns in the motion data file

    // Read centre of rotation, translation, and rotation values for each height
    for (int heightIdx = 0; heightIdx < N_heights; ++heightIdx)
    {
        List<double> tmp_height(N_columns);
        for (int i = 0; i < N_columns; ++i)
        {
            inputFile >> tmp_height[i];
        }
        motion.append(tmp_height);
    }

    inputFile.close();
}

// Function to associate each point in the mesh to a specific height/radius by calculating height index and fraction of each point
std::pair<int, double> Foam::prescribedSLERPMotionSolver::p2HeightAssociation(Vector<double> &point, List<List<double>> &motion, double N_heights)
{
    // Multiply each coordinate by the corresponding coefficient
    double x = coeffX_ * point[0];
    double y = coeffY_ * point[1];
    double z = coeffZ_ * point[2];

    // Calculate the height/radius of each point
    double magnitude = sqrt(x * x + y * y + z * z);

    // Determine the sign of the height/radius based on the dominant coefficient
    double values[] = {x, y, z};
    double absValues[] = {std::abs(x), std::abs(y), std::abs(z)};
    int dominantIndex = std::distance(absValues, std::max_element(absValues, absValues + 3));
    double sign = (values[dominantIndex] >= 0) ? 1.0 : -1.0;

    // Apply the sign to the magnitude to get the signed height
    double height = sign * magnitude;

    // Assign each points in the mesh to a height index based on the height/radius
    int index = 0;
    double fraction = 0.0;

    // Handle edge cases
    if (height < motion[0][0])
    {
        index = 0;
        fraction = 0.0;
    }
    else if (height > motion[N_heights - 1][0])
    {
        index = N_heights - 2;
        fraction = 1.0;
    }
    else
    {
        // Find the heights interval
        for (int i = 0; i < N_heights - 1; i++)
        {
            if (height >= motion[i][0] && height <= motion[i + 1][0])
            {
                index = i;
                fraction = (height - motion[i][0]) / (motion[i + 1][0] - motion[i][0]);
                break;
            }
        }
    }

    return std::pair<int, double>(index, fraction);
}

// Function to calculate the time fraction needed for interpolation between motion data at two consecutive time-steps
// This function will read the next motion data file if the current time is greater than the last time-step saved
double Foam::prescribedSLERPMotionSolver::timeFractionCalculation(double timeValue) const
{
    if (timeValue <= timeArray[curTimeIndexRead_ + 1])
    {
        // Interpolate between the two closest time values
        double t1 = timeArray[curTimeIndexRead_];
        double t2 = timeArray[curTimeIndexRead_ + 1];
        double timeFraction = (timeValue - t1) / (t2 - t1);

        return timeFraction;
    }
    else
    {
        // Read the next file
        curTimeIndexRead_++;

        // Safeguard against going out of bounds
        if (curTimeIndexRead_ >= timeArray.size() - 1)
        {
            curTimeIndexRead_ = timeArray.size() - 2;
        }

        // Interpolate between the two closest time values
        double t1 = timeArray[curTimeIndexRead_];
        double t2 = timeArray[curTimeIndexRead_ + 1];
        double timeFraction = (timeValue - t1) / (t2 - t1);

        // Read the next file
        motionData = nextMotionData;

        // Clear the nextMotionData list
        nextMotionData.clear();

        double timeDouble = timeArray[curTimeIndexRead_ + 1];

        // Use a stringstream to format the time value with fixed precision
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(9) << timeDouble; // Set precision to 9 decimal places
        std::string filename = path_ + "/motionData_" + stream.str() + ".dat";

        readMotionData(filename, nextMotionData);

        return timeFraction;
    }
}

// Function to compute the transformation septernion for a point in the mesh based on its height index and fraction
Foam::septernion Foam::prescribedSLERPMotionSolver::computePointTransformation(int index, double fraction, double m, List<List<double>> &motion) const
{
    // Interpolate translation
    vector translation1(motion[index][4], motion[index][5], motion[index][6]);
    vector translation2(motion[index + 1][4], motion[index + 1][5], motion[index + 1][6]);
    vector interpolatedTranslation = ((1 - fraction) * translation1 + fraction * translation2) * m;

    // Extract and interpolate rotation (Euler angles)
    vector rotation1(motion[index][7], motion[index][8], motion[index][9]);
    vector rotation2(motion[index + 1][7], motion[index + 1][8], motion[index + 1][9]);
    vector interpolatedRotation = ((1 - fraction) * rotation1 + fraction * rotation2) * m;

    // Convert rotation to quaternion
    quaternion rotationQuat(quaternion::XYZ, interpolatedRotation * degToRad());

    // Extract and interpolate centre of rotation
    vector centre1(motion[index][1], motion[index][2], motion[index][3]);
    vector centre2(motion[index + 1][1], motion[index + 1][2], motion[index + 1][3]);
    vector interpolatedCentre = (1 - fraction) * centre1 + fraction * centre2;

    // Construct septernion for combined translation and rotation
    septernion transformation(septernion(-interpolatedCentre + -interpolatedTranslation) * rotationQuat * septernion(interpolatedCentre));

    return transformation;
}

Foam::tmp<Foam::pointField>
Foam::prescribedSLERPMotionSolver::curPoints() const
{
    const Time &t = this->db().time();

    if (motionType_ == "harmonic")
    {
        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                int index = p2heightIdx_[pointi];
                double fraction = p2heightFraction_[pointi];

                double m = amplitude_ * sin(omega_ * t.value());  // Calculate the motion multiplier for harmonic motion

                // Calculate the transformation septernion for the point
                septernion transformation = computePointTransformation(index, fraction, m, motionData);

                // Calculate the scaled point displacement from the point transformation

                // Use solid-body motion where scale = 1
                if (scale_[pointi] > 1 - SMALL)
                {
                    pointDisplacement_.primitiveFieldRef()[pointi] = transformation.transformPoint(points0_[pointi]) - points0_[pointi];
                }

                // Slerp septernion interpolation: scaled septernion calculation
                else
                {
                    septernion ss(slerp(septernion::I, transformation, scale_[pointi]));

                    pointDisplacement_.primitiveFieldRef()[pointi] = ss.transformPoint(points0_[pointi]) - points0_[pointi];
                }
            }
            else
            {
                pointDisplacement_.primitiveFieldRef()[pointi] = Foam::vector(0, 0, 0);
            }
        }

        // Calculate the updated points location
        tmp<pointField> newPoints(
            points0() + pointDisplacement_.primitiveField());
            
        return newPoints;
    }

    else if (motionType_ == "nonHarmonic")
    {
        double timeFraction = timeFractionCalculation(t.value());

        motionDataInterpolated = motionData;

        // Calculate the interpolated motion data at the current time-step from the clostest two time-steps saved in timeArray
        forAll(motionDataInterpolated, i)
        {
            motionDataInterpolated[i] = (1 - timeFraction) * motionData[i] + timeFraction * nextMotionData[i];
        }

        forAll(points0(), pointi)
        {
            if (scale_[pointi] > SMALL)
            {
                int index = p2heightIdx_[pointi];
                double fraction = p2heightFraction_[pointi];

                double m = 1.0;  // For non-harmonic motion, the multiplier is always 1

                // Calculate the transformation septernion for the point
                septernion transformation = computePointTransformation(index, fraction, m, motionDataInterpolated);

                // Calculate the scaled point displacement from the point transformation

                // Use solid-body motion where scale = 1
                if (scale_[pointi] > 1 - SMALL)
                {
                    pointDisplacement_.primitiveFieldRef()[pointi] = transformation.transformPoint(points0_[pointi]) - points0_[pointi];
                }
                
                // Slerp septernion interpolation: scaled septernion calculation
                else
                {
                    septernion ss(slerp(septernion::I, transformation, scale_[pointi]));

                    pointDisplacement_.primitiveFieldRef()[pointi] = ss.transformPoint(points0_[pointi]) - points0_[pointi];
                }
            }
            else
            {
                pointDisplacement_.primitiveFieldRef()[pointi] = Foam::vector(0, 0, 0);
            }
        }

        // Calculate the updated points location
        tmp<pointField> newPoints(
            points0() + pointDisplacement_.primitiveField());
            
        return newPoints;
    }

    else
    {
        return 0;
    }
}

void Foam::prescribedSLERPMotionSolver::solve()
{

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Displacement has changed. Update boundary conditions
    pointConstraints::New(
        pointDisplacement_.mesh())
        .constrainDisplacement(pointDisplacement_);
        
}


// ************************************************************************* //

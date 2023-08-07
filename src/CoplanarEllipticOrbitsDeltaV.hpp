//
// Created by julia on 28.07.23.
//

#ifndef VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP
#define VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP
#include <cmath>
#include <array>
#include <optional>
#include "Types.h"
namespace Maneuvers {

    /** angleBetweenImpulseAndTransversal1 -- angle Between 1st Impulse and Transversal Velocity [0;360] deg
     *  angleBetweenImpulseAndTransversal2 -- angle Between 2nd Impulse and Transversal Velocity [0;360] deg
     *  polarEq1 -- polar Equation of the start Orbit (meters)
     *  polarEq2 -- polar Equation of the final Orbit (meters)
     *  periapseA1 -- periapse of the start Orbit (1 / sqrt(meters))
     *  periapseA2 -- periapse of the final Orbit (1 / sqrt(meters))
     *  semimajorAxisTransferOrbit -- semimajor Axis of the Transfer Orbit (meters)
     */

    struct parametersForTransfer {
        scalar polarEq1;
        scalar polarEq2;
        scalar angleBetweenImpulseAndTransversal1;
        scalar angleBetweenImpulseAndTransversal2;
        scalar periapseA1;
        scalar periapseA2;
        scalar semimajorAxisTransferOrbit;
    };

    bool isCorrectValues(const scalar cosDeltaTrueAnomaly1AscendNode, const scalar cosDeltaTrueAnomaly2AscendNode,
                         const scalar angleBetweenImpulseAndTransversalVelocity1, const scalar angleBetweenImpulseAndTransversalVelocity2,
                         const scalar polarEq1, const scalar polarEq2){

        if (std::abs(cosDeltaTrueAnomaly1AscendNode - cosDeltaTrueAnomaly2AscendNode)
           < std::numeric_limits<scalar>::epsilon() and std::abs(polarEq1 - polarEq2)
                                                        > std::numeric_limits<scalar>::epsilon()){
             return false;
        }
        auto is_correct {[](const scalar angleBetweenImpulseAndTransversalVelocity)
        {return std::abs(angleBetweenImpulseAndTransversalVelocity - static_cast<scalar>(0)) < std::numeric_limits<scalar>::epsilon() or
            std::abs(angleBetweenImpulseAndTransversalVelocity - M_PI) < std::numeric_limits<scalar>::epsilon() or
            std::abs(angleBetweenImpulseAndTransversalVelocity - 2 * M_PI) < std::numeric_limits<scalar>::epsilon(); }};

        if (!is_correct(angleBetweenImpulseAndTransversalVelocity1)
            or !is_correct(angleBetweenImpulseAndTransversalVelocity2))
        {return false;}
    }

    [[nodiscard]] std::optional<parametersForTransfer> computeParamsForTransfer(const EllipticOrbit &first, const EllipticOrbit &final,
                                                                 const scalar trueAnomaly1, const scalar trueAnomaly2,
                                                                 const scalar ascendNodeTransferOrbit, const scalar angleBetweenImpulseAndTransversal1,
                                                                 const scalar angleBetweenImpulseAndTransversal2) {
        const scalar semiminorAxis1 = (1 - first.e) * first.a;
        const scalar semiminorAxis2 = (1 - final.e) * final.a;

        const scalar cosDeltaTrueAnomaly1AscendNode = std::cos(trueAnomaly1 - ascendNodeTransferOrbit);
        const scalar cosDeltaTrueAnomaly2AscendNode = std::cos(trueAnomaly2 - ascendNodeTransferOrbit);

        // Polar equation of an orbit
        // u = a + b cos(trueAnomaly - ascendNode)
        const scalar polarEq1 = first.a + first.b * cosDeltaTrueAnomaly1AscendNode;
        const scalar polarEq2 = final.a + final.b * cosDeltaTrueAnomaly2AscendNode;

        if(!isCorrectValues(cosDeltaTrueAnomaly1AscendNode, cosDeltaTrueAnomaly2AscendNode,
                            angleBetweenImpulseAndTransversal1, angleBetweenImpulseAndTransversal2,
                            polarEq1, polarEq2)){
            return std::nullopt;
        }

        const scalar semiminorAxisTransferOrbit = (polarEq1 - polarEq2) /
                (cosDeltaTrueAnomaly1AscendNode - cosDeltaTrueAnomaly2AscendNode);
        const scalar semimajorAxisTransferOrbit = polarEq1 -
                semiminorAxisTransferOrbit * cosDeltaTrueAnomaly1AscendNode;

        // A1 -- the periapse of 1st Orbit
        // angleBetweenImpulseAndTransversal1 is measured counterclockwise between 0 and 360 deg
        const scalar sinDeltaTrueAnomaly1AscendNode1 = semiminorAxis1 * std::sin(trueAnomaly1 - first.w);
        const scalar periapseA1 = polarEq1 / std::sqrt(first.a) -
                sinDeltaTrueAnomaly1AscendNode1 * std::atan(angleBetweenImpulseAndTransversal1);

        // A2 -- the periapse of 2nd Orbit
        // angleBetweenImpulseAndTransversal2 is measured counterclockwise between 0 and 360 deg
        const scalar sinDeltaTrueAnomaly2AscendNode2 = semiminorAxis2 * std::sin(trueAnomaly2 - final.w);
        const scalar periapseA2 = polarEq2 / std::sqrt(final.a) -
                sinDeltaTrueAnomaly2AscendNode2 * std::atan(angleBetweenImpulseAndTransversal2);

        return parametersForTransfer{polarEq1, polarEq2,
                                     angleBetweenImpulseAndTransversal1, angleBetweenImpulseAndTransversal2,
                                     periapseA1, periapseA2,
                                     semimajorAxisTransferOrbit};
    }

    scalar getFunctionForMinimum(const parametersForTransfer &paramsForTransfer){

       const scalar H = ((paramsForTransfer.polarEq1 + paramsForTransfer.semimajorAxisTransferOrbit) /
               paramsForTransfer.periapseA1 / std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit) + 1) * std::cos(paramsForTransfer.angleBetweenImpulseAndTransversal1) -
               ((paramsForTransfer.polarEq2 + paramsForTransfer.semimajorAxisTransferOrbit) /
                paramsForTransfer.periapseA2 / std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit) + 1) * std::cos(paramsForTransfer.angleBetweenImpulseAndTransversal2);

       const scalar G = (paramsForTransfer.polarEq1 / paramsForTransfer.periapseA1 - paramsForTransfer.periapseA1) * std::sin(paramsForTransfer.angleBetweenImpulseAndTransversal1) -
               (paramsForTransfer.polarEq2 / paramsForTransfer.periapseA2 - paramsForTransfer.periapseA2) * std::sin(paramsForTransfer.angleBetweenImpulseAndTransversal2);

       const scalar F = (paramsForTransfer.polarEq1 - paramsForTransfer.semimajorAxisTransferOrbit) *
               (1 + std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit) / paramsForTransfer.periapseA1) * std::cos(paramsForTransfer.angleBetweenImpulseAndTransversal1) +
               (paramsForTransfer.polarEq1 - paramsForTransfer.periapseA1 * std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit)) * std::sin(paramsForTransfer.angleBetweenImpulseAndTransversal1) -
               (paramsForTransfer.polarEq2 - paramsForTransfer.semimajorAxisTransferOrbit) *
               (1 + std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit) / paramsForTransfer.periapseA2) * std::cos(paramsForTransfer.angleBetweenImpulseAndTransversal2) +
               (paramsForTransfer.polarEq2 - paramsForTransfer.periapseA2 * std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit)) * std::sin(paramsForTransfer.angleBetweenImpulseAndTransversal2);

       return std::sqrt(H * H + G * G + F * F);

    }

    std::array<scalar, 3> computeParametersForVelocity(const EllipticOrbit &first, const EllipticOrbit &final,
                                                       const scalar hTrueAnomaly1, const scalar hTrueAnomaly2,
                                                       const scalar hAscendNode, const scalar hAngle1, const scalar hAngle2){
        auto params = computeParamsForTransfer(first, final, hTrueAnomaly1,
                                                                hTrueAnomaly2, hAscendNode,
                                                                hAngle1, hAngle2).value();

        scalar functionForMin = getFunctionForMinimum(params);



    }



}

#endif //VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP

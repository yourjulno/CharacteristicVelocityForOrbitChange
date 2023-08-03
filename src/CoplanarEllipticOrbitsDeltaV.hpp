//
// Created by julia on 28.07.23.
//

#ifndef VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP
#define VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP
#include <cmath>
#include <array>
#include "Types.h"
namespace Maneuvers {

    /** angleBetweenImpulseAndTransversal1 -- angle Between 1st Impulse and Transversal Velocity [0;360] deg
     *  angleBetweenImpulseAndTransversal2 -- angle Between 2nd Impulse and Transversal Velocity [0;360] deg
     *  polarEq1 -- polar Equation of the start Orbit (meters)
     *  polarEq2 -- polar Equation of the final Orbit (meters)
     *  periapseA1 -- periapse of the start Orbit (meters)
     *  periapseA2 -- periapse of the final Orbit (meters)
     *  semimajorAxisTransferOrbit -- semimajor Axis of the Transfer Orbit (meters)
     */

    struct parametersForTransfer {
        scalar angleBetweenImpulseAndTransversal1;
        scalar angleBetweenImpulseAndTransversal2;
        scalar polarEq1;
        scalar polarEq2;
        scalar periapseA1;
        scalar periapseA2;
        scalar semimajorAxisTransferOrbit;
    };


    [[nodiscard]] parametersForTransfer computeParamsForTransfer(const EllipticOrbit &first, const EllipticOrbit &final,
                                   const scalar trueAnomaly1, const scalar trueAnomaly2,
                                   const scalar ascendNodeTransferOrbit) {
        const scalar semiminorAxis1 = (1 - first.e) * first.a;
        const scalar semiminorAxis2 = (1 - final.e) * final.a;
        const scalar cosDeltaTrueAnomaly1AscendNode = std::cos(trueAnomaly1 - ascendNodeTransferOrbit);
        const scalar cosDeltaTrueAnomaly2AscendNode = std::cos(trueAnomaly2 - ascendNodeTransferOrbit);


        // Polar equation of an orbit
        // u = a + b cos(trueAnomaly - ascendNode)
        const scalar polarEq1 = first.a + first.b * cosDeltaTrueAnomaly1AscendNode;
        const scalar polarEq2 = final.a + final.b * cosDeltaTrueAnomaly2AscendNode;

        const scalar semiminorAxisTransferOrbit = (polarEq1 - polarEq2) /
                                                  (cosDeltaTrueAnomaly1AscendNode - cosDeltaTrueAnomaly2AscendNode);
        const scalar semimajorAxisTransferOrbit = polarEq1 -
                                                  semiminorAxisTransferOrbit * cosDeltaTrueAnomaly1AscendNode;

        // A1 -- the periapse of 1st Orbit
        // angleBetweenImpulseAndTransversal1 is measured counterclockwise between 0 and 360 deg
        const scalar sinDeltaTrueAnomaly1AscendNode1 = semiminorAxis1 * std::sin(trueAnomaly1 - first.w);
        const scalar sinDeltaTrueAnomaly1AscendNode = semiminorAxisTransferOrbit * std::sin(trueAnomaly1 - ascendNodeTransferOrbit);

        const scalar periapseA1 = polarEq1 * (sinDeltaTrueAnomaly1AscendNode1 - sinDeltaTrueAnomaly1AscendNode) /
                                  (std::sqrt(semimajorAxisTransferOrbit) * sinDeltaTrueAnomaly1AscendNode1 - std::sqrt(first.a));
        const scalar angleBetweenImpulseAndTransversal1 = std::atan(sinDeltaTrueAnomaly1AscendNode1) /
                                                          (polarEq1 - periapseA1 * std::sqrt(first.a));

        // A2 -- the periapse of 2nd Orbit
        // angleBetweenImpulseAndTransversal2 is measured counterclockwise between 0 and 360 deg
        const scalar sinDeltaTrueAnomaly2AscendNode = semiminorAxisTransferOrbit * std::sin(trueAnomaly2 - ascendNodeTransferOrbit);
        const scalar sinDeltaTrueAnomaly2AscendNode2 = semiminorAxis2 * std::sin(trueAnomaly2 - final.w);

        const scalar periapseA2 = polarEq2 * (sinDeltaTrueAnomaly2AscendNode2 - sinDeltaTrueAnomaly2AscendNode) /
                                  (std::sqrt(semimajorAxisTransferOrbit) * sinDeltaTrueAnomaly2AscendNode2 - std::sqrt(final.a));
        const scalar angleBetweenImpulseAndTransversal2 = std::atan(sinDeltaTrueAnomaly2AscendNode2) /
                                                          (polarEq2 - periapseA2 * std::sqrt(final.a));
        return parametersForTransfer{angleBetweenImpulseAndTransversal1, angleBetweenImpulseAndTransversal2,
                                     polarEq1, polarEq2,
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
                                                       const scalar hAscendNode){
        auto params = computeParamsForTransfer(first, final, hTrueAnomaly1,
                                               hTrueAnomaly2, hAscendNode);
        scalar functionForMin = getFunctionForMinimum(params);



    }



}

#endif //VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP

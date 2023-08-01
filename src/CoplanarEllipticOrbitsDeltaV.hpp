//
// Created by julia on 28.07.23.
//

#ifndef VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP
#define VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP
#include <cmath>

#include "Types.h"
namespace Maneuvers {

    scalar simpleIterationMethodOneIter(scalar digit, const scalar semimajorAxis1, const scalar semiminorAxis1,
                                 const scalar semimajorAxis2, const scalar semiminorAxis2){
        const scalar digitPow2 = digit * digit;
        const scalar const1 = digitPow2 * digitPow2 - 2 * digitPow2 * semimajorAxis1
                + semimajorAxis1 * semimajorAxis1 - semiminorAxis1 * semiminorAxis1;
        const scalar const2 = digitPow2 * digitPow2 - 2 * digitPow2 * semimajorAxis2
                + semimajorAxis2 * semimajorAxis2 - semiminorAxis2 * semiminorAxis2;
        const scalar deltaConsts = const1 - const2;
        digit = (digitPow2 * deltaConsts + semiminorAxis2 * const1 - semimajorAxis1 * const2) /
                (2 * (std::sqrt(semimajorAxis2) * const1 - std::sqrt(semimajorAxis1) * const2));
        return digit;
    }
    scalar simpleIterationMethod( const scalar semimajorAxis1, const scalar semiminorAxis1,
                                 const scalar semimajorAxis2, const scalar semiminorAxis2){
        scalar startPeriapse = semimajorAxis1;
        startPeriapse = simpleIterationMethodOneIter(startPeriapse,semimajorAxis1, semiminorAxis1,
                    semimajorAxis2, semiminorAxis2);
        while (std::abs(startPeriapse - simpleIterationMethodOneIter(startPeriapse,semimajorAxis1, semiminorAxis1,
                                                                     semimajorAxis2, semiminorAxis2)) > std::numeric_limits<scalar>::epsilon()){
            startPeriapse = simpleIterationMethodOneIter(startPeriapse,semimajorAxis1, semiminorAxis1,
                                                         semimajorAxis2, semiminorAxis2);
        }
        return startPeriapse;
    }
    scalar computeSemimajorAxisTransferOrbit(const Orbit2 &first, const Orbit2 &final,
                                             const scalar trueAnomaly1, const scalar trueAnomaly2,
                                             const scalar ascendNodeTransferOrbit){

        const scalar cosDeltaTrueAnomaly1AscendNode = std::cos(trueAnomaly1 - ascendNodeTransferOrbit);
        const scalar cosDeltaTrueAnomaly2AscendNode = std::cos(trueAnomaly2 - ascendNodeTransferOrbit);

        if (std::abs(cosDeltaTrueAnomaly1AscendNode - cosDeltaTrueAnomaly2AscendNode) < std::numeric_limits<scalar>::epsilon()){

        }
        // Polar equation of an orbit
        // u = a + b cos(trueAnomaly - ascendNode)
        const scalar polarEq1 = first.a + first.b * cosDeltaTrueAnomaly1AscendNode;
        const scalar polarEq2 = final.a + final.b * cosDeltaTrueAnomaly2AscendNode;

        const scalar semiminorAxisTransferOrbit = (polarEq1 - polarEq2) /
                                                  (cosDeltaTrueAnomaly1AscendNode - cosDeltaTrueAnomaly2AscendNode);
        const scalar semimajorAxisTransferOrbit = polarEq1 -
                                                  semiminorAxisTransferOrbit * cosDeltaTrueAnomaly1AscendNode;
        return semimajorAxisTransferOrbit;
    }

    void coplanarEllipticOrbitsDeltaV(const Orbit2 &first, const Orbit2 &final,
                                   const scalar trueAnomaly1, const scalar trueAnomaly2,
                                   const scalar ascendNodeTransferOrbit, const scalar gravParameter){
        const scalar semiminorAxis1 = (1 - first.e) * first.a;
        const scalar semiminorAxis2 = (1 - final.e) * final.a;
        const scalar cosDeltaTrueAnomalyAscendNode1 = std::cos(trueAnomaly1 - ascendNodeTransferOrbit);
        const scalar cosDeltaTrueAnomalyAscendNode2 = std::cos(trueAnomaly2 - ascendNodeTransferOrbit);

        // Polar equation of an orbit
        // u = a + b cos(trueAnomaly - ascendNode)
        const scalar polarEq1 = first.a + first.b * cosDeltaTrueAnomalyAscendNode1;
        const scalar polarEq2 = final.a + final.b * cosDeltaTrueAnomalyAscendNode2;

        const scalar semiminorAxisTransferOrbit = (polarEq1 - polarEq2) /
                                                  (cosDeltaTrueAnomalyAscendNode1 - cosDeltaTrueAnomalyAscendNode2);
        const scalar semimajorAxisTransferOrbit = polarEq1 -
                                                  semiminorAxisTransferOrbit * cosDeltaTrueAnomalyAscendNode1;

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

    }

}
#endif //VELOCITY_COPLANARELLIPTICORBITSDELTAV_HPP

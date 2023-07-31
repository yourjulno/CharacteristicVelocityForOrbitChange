//
// Created by julia on 20.07.23.
//

#ifndef VELOCITY_COMPLANARCIRCALORBITSDELTAV_HPP
#define VELOCITY_COMPLANARCIRCALORBITSDELTAV_HPP
#include <cmath>
#include "Types.h"
namespace Maneuvers {


    scalar coplanarManeuverDeltaV(const scalar velocityInApogree,
                                                const Orbit &start,
                                                const Orbit &final){
        const scalar deltaE_x = final.e * std::cos(final.periapsisArg)
                - start.e * std::cos(start.periapsisArg);
        const scalar deltaE_y = final.e * std::sin(final.periapsisArg)
                - start.e * std::sin(start.e);
        const scalar deltaE = std::sqrt(deltaE_x * deltaE_x
                + deltaE_y * deltaE_y);
        const scalar r0 = (final.a + start.a) / 2;

        const scalar deltaSemimajorAxis = (final.a - start.a) / r0;
        const scalar absDeltaSemimajorAxis = std::abs(deltaSemimajorAxis);
        const scalar delta1 = 0.25 * absDeltaSemimajorAxis * velocityInApogree;
        const scalar resultdV = absDeltaSemimajorAxis > deltaE ?
                absDeltaSemimajorAxis / 2 * velocityInApogree :
                std::sqrt(deltaE * deltaE - 0.75 * deltaSemimajorAxis * deltaSemimajorAxis) * velocityInApogree;

        return resultdV;
    }

};
#endif //VELOCITY_COMPLANARCIRCALORBITSDELTAV_HPP

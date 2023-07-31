//
// Created by julia on 15.07.23.
//

#ifndef VELOCITY_PLANECHANGEMANEUVER_HPP
#define VELOCITY_PLANECHANGEMANEUVER_HPP

#include <cmath>

#include "Types.h"

namespace Maneuvers {


    [[nodiscard]] scalar planeRotationDeltaV(const scalar momentum1,
                                        const scalar inclination1, const scalar inclination2,
                                        const scalar ascendNode1, const scalar ascendNode2,
                                        const scalar apocenter) noexcept {


        const scalar cosdeltaAscendNode = std::cos(ascendNode2 - ascendNode1);

        const scalar cosdeltaIncl = std::cos(inclination1) * std::cos(inclination2)
                                    + std::sin(inclination1) * std::sin(inclination2) * cosdeltaAscendNode;

        return std::sqrt(2 * momentum1 * momentum1 * (1 - cosdeltaIncl)) / apocenter;
    }

}  // namespace Maneuvers

#endif // VELOCITY_PLANECHANGEMANEUVER_HPP

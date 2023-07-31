//
// Created by julia on 19.07.23.
//
#include <gtest/gtest.h>
#include "../src/PlaneChangeManeuver.hpp"

using namespace Maneuvers;

TEST(PlaneRot, SECOND) {
    const scalar inclination1 = 0;
    const scalar inclination2 = M_PI / 2;
    const scalar ascendNode1 = M_PI;
    const scalar ascendNode2 = M_PI;
    const scalar apocenter = 1;
    const scalar momentum1 = 1;
    const scalar deltaV = planeRotationDeltaV(momentum1, inclination1, inclination2,
                                        ascendNode1, ascendNode2, apocenter);
    std::cout << deltaV;
    //ASSERT_NEAR(std::sqrt(static_cast<scalar>(2)), deltaV, std::numeric_limits<scalar>::epsilon());
}

TEST(PlaneRot, THIRD) {
    const scalar nu = 3.986028e14;

    const scalar inclination1 = 0;
    const scalar inclination2 = 0;
    const scalar ascendNode1 = 1;
    const scalar ascendNode2 = 1;
    const scalar apocenter = 67;
    const scalar momentum1 = 1;
    const scalar deltaV = planeRotationDeltaV(momentum1, inclination1, inclination2,
                                        ascendNode1, ascendNode2, apocenter);
    //std::cout << deltaV;
    ASSERT_NEAR(100, deltaV, std::numeric_limits<scalar>::epsilon());

}

TEST(PlaneRot, FOURTH) {

    const scalar inclination1 = 0;
    const scalar inclination2 = M_PI;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;
    const scalar apocenter = 1;
    const scalar momentum1 = 1;
    const scalar deltaV = planeRotationDeltaV(momentum1, inclination1,
                                        inclination2, ascendNode1, ascendNode2, apocenter);
    std::cout << deltaV;
    //ASSERT_NEAR(2, deltaV, std::numeric_limits<scalar>::epsilon());


}

TEST(PlaneRot, FIFTH) {

    const scalar inclination1 = M_PI / 2;
    const scalar inclination2 = M_PI / 2;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = M_PI / 2 ;
    const scalar apocenter = 1;
    const scalar momentum1 = 1;
    const scalar deltaV = planeRotationDeltaV(momentum1, inclination1,
                                        inclination2, ascendNode1, ascendNode2, apocenter);
    ASSERT_NEAR(std::sqrt(static_cast<scalar>(2)), deltaV, std::numeric_limits<scalar>::epsilon());

}

TEST(PlaneRot, EXMAPLE) {
    const scalar nu = 3.986028e14;

    const scalar deg = M_PI / 180;
    const scalar inclination1 = 97*deg;
    const scalar inclination2 = 97*deg;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 5 * deg;
    const scalar apocenter = 500e3 + 6731e3;
    const scalar V = std::sqrt(nu / apocenter);
    const scalar momentum1 = apocenter * V;
    const scalar deltaV = planeRotationDeltaV(momentum1, inclination1,
                                              inclination2, ascendNode1, ascendNode2, apocenter);
    //ASSERT_NEAR(std::sqrt(static_cast<scalar>(2)), deltaV, std::numeric_limits<scalar>::epsilon());
    std::cout << deltaV;
}



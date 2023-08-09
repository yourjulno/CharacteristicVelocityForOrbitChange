//
// Created by julia on 09.08.23.
//
#include <gtest/gtest.h>
#include "../src/CoplanarCircularOrbit.hpp"
#include "../src/CoplanarEllipticOrbit.hpp"
using namespace Maneuvers;

constexpr scalar deg = M_PI / 180;

//not intersect
TEST(transferOrbits, FIRST) {
    const scalar nu = 3.986028e14;
    const scalar EarthRadius = 6731e3;
    const scalar semimajorAxis1 = 500e3 + EarthRadius;
    const scalar semimajorAxis2 = 500e3 + EarthRadius;
    const scalar eccentricity1 = 0;
    const scalar eccentricity2 = 0;
    const scalar inclination1 = 97 * deg;
    const scalar inclination2 = 97 * deg;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 5 * deg;
    const scalar periapsisArg1 = M_PI / 2;
    const scalar periapsisArg2 = M_PI / 2;

    const CircularOrbit o1 = {semimajorAxis1, eccentricity1, inclination1,
                              ascendNode1, periapsisArg1};
    const CircularOrbit o2 = {semimajorAxis2, eccentricity2, inclination2,
                              ascendNode2, periapsisArg2};
    const scalar velocityInApogree = std::sqrt(nu / semimajorAxis1);
    const scalar resultDeltaV = coplanarManeuverDeltaV(velocityInApogree,
                                                       o1, o2);
    std::cout << resultDeltaV;
}

TEST(transferOrbits, SECOND) {

    const scalar nu = 3.986028e14;
    const scalar semimajorAxis1 = 6566e3;
    const scalar semimajorAxis2 = 6721e3;
    const scalar eccentricity1 = 0.00228;
    const scalar eccentricity2 = 0.00149;
    const scalar inclination1 = 51.7 * deg;
    const scalar inclination2 = 51.69 * deg;
    const scalar ascendNode1 = 17.49 * deg;
    const scalar ascendNode2 = 17.5 * deg;
    const scalar pericenterArg1 = 20 * deg;
    const scalar pericenterArg2 = 150 * deg;

    const CircularOrbit o1 = {semimajorAxis1, eccentricity1, inclination1,
                              ascendNode1, pericenterArg1};
    const CircularOrbit o2 = {semimajorAxis2, eccentricity2, inclination2,
                              ascendNode2, pericenterArg2};
    const scalar apogree2 = (semimajorAxis2 + semimajorAxis1) / 2;
    const scalar velocityInApogree = std::sqrt(nu / apogree2);
    const scalar resultDeltaV = coplanarManeuverDeltaV(velocityInApogree,
                                                       o1, o2);
    const scalar realResultDeltaV = 90.36;
    ASSERT_NEAR(resultDeltaV, realResultDeltaV, 1.2);
}

TEST(transferOrbits, THIRD) {

    const scalar nu = 3.986028e14;
    const scalar semimajorAxis1 = 590.0e3;
    const scalar semimajorAxis2 = 600.0e3;
    const scalar eccentricity1 = 0;
    const scalar eccentricity2 = 0;
    const scalar inclination1 = 0;
    const scalar inclination2 = 0;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;
    const scalar periapsisArg1 = 0;
    const scalar periapsisArg2 = 0;

    const CircularOrbit o1 = {semimajorAxis1, eccentricity1, inclination1,
                              ascendNode1, periapsisArg1};
    const CircularOrbit o2 = {semimajorAxis2, eccentricity2, inclination2,
                              ascendNode2, periapsisArg2};
    const scalar velocityOnFirstrbit  = std::sqrt(nu / semimajorAxis1);
    const scalar velocityOnSEcondOrbit = std::sqrt(nu / semimajorAxis2);
    const scalar r0 = (semimajorAxis1 + semimajorAxis2) / 2 ;
    const scalar attitude = semimajorAxis2  / semimajorAxis1;
    const scalar ellipticVelocity1 = std::sqrt(nu * 2 * semimajorAxis2 /
                                               semimajorAxis1 / (semimajorAxis1 + semimajorAxis2));
    const scalar ellipticVelocity2 = std::sqrt(nu * 2 * semimajorAxis1 /
                                               semimajorAxis2 / (semimajorAxis1 + semimajorAxis2));
    const scalar firstdelta = ellipticVelocity1 - velocityOnFirstrbit;

    const scalar d = velocityOnFirstrbit - velocityOnSEcondOrbit;
    const scalar seconddelta = velocityOnSEcondOrbit - ellipticVelocity2;
    const scalar deltaFirstV = ellipticVelocity1 - velocityOnFirstrbit;
    const scalar deltaSecondV = velocityOnSEcondOrbit * (1 - std::sqrt(static_cast<scalar>(2) * semimajorAxis1
                                                                       / (semimajorAxis2 + semimajorAxis1)));

    const scalar realResultDeltaV =  firstdelta + seconddelta;

    const scalar velocityInApogree = std::sqrt(nu / r0);
    const scalar resultDeltaV = coplanarManeuverDeltaV(velocityInApogree,
                                                       o1, o2);

    //ASSERT_NEAR(resultDeltaV, realResultDeltaV, 1.2);
    std::cout << realResultDeltaV << " " << seconddelta << " " << resultDeltaV;
}
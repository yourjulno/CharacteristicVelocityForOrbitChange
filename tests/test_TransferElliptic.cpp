//
// Created by julia on 08.08.23.
//
//
// Created by julia on 22.07.23.
//
#include <gtest/gtest.h>
#include "../src/CoplanarCircularOrbit.hpp"
#include "../src/CoplanarEllipticOrbit.hpp"
using namespace Maneuvers;

constexpr scalar deg = M_PI / 180;


TEST(CorrectParametres, FIRST){
    const scalar cosDeltaTrueAnomaly1AscendNode = 0;
    const scalar cosDeltaTrueAnomaly2AscendNode = 0;
    const scalar angleBetweenImpulseAndTransversalVelocity1 = M_PI +
                                                              std::numeric_limits<scalar>::epsilon();
    const scalar angleBetweenImpulseAndTransversalVelocity2 = M_PI;
    const scalar polarEq1 = 0;
    const scalar polarEq2 = 0;
    const auto result = isCorrectValues(cosDeltaTrueAnomaly1AscendNode, cosDeltaTrueAnomaly2AscendNode,
                                 angleBetweenImpulseAndTransversalVelocity1, angleBetweenImpulseAndTransversalVelocity2,
                                 polarEq1, polarEq2);
    ASSERT_EQ(result, false);
}

TEST(CorrectParametres, SECOND){
    const scalar cosDeltaTrueAnomaly1AscendNode = 100;
    const scalar cosDeltaTrueAnomaly2AscendNode = 0;
    const scalar angleBetweenImpulseAndTransversalVelocity1 = M_PI +
                                                              2 * std::numeric_limits<scalar>::epsilon();
    const scalar angleBetweenImpulseAndTransversalVelocity2 = M_PI;
    const scalar polarEq1 = 0;
    const scalar polarEq2 = 0;
    const auto result = isCorrectValues(cosDeltaTrueAnomaly1AscendNode, cosDeltaTrueAnomaly2AscendNode,
                                        angleBetweenImpulseAndTransversalVelocity1, angleBetweenImpulseAndTransversalVelocity2,
                                        polarEq1, polarEq2);
    ASSERT_EQ(result, true);
}

TEST(computeParametersForVelocity, FIRST){
    const scalar nu = 3.986028e14;
    const scalar semimajorAxis1 = 3;
    const scalar semiminorAxis1 = 1;
    const scalar semimajorAxis2 = 2;
    const scalar semiminorAxis2 = 1;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;
    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};

    const std::size_t count = 500;
    const CountsForSteps counts = {count, count,
                                   count, count, count};
    const auto result = computeParametersForVelocity_(first, final, counts);
    std::cout << result.value().semimajorAxisTransferOrbit << " ";
    std::cout << result.value().periapseA1 << " ";
    std::cout << result.value().angleBetweenImpulseAndTransversal1 * 180 / M_PI << " ";
    std::cout << result.value().functional << " ";
    std::cout << result.value().angleBetweenImpulseAndTransversal2 * 180 / M_PI;

}

TEST(computeVelocity, FIRST){
    const scalar nu = 3.986028e14;
    const scalar G =  6.67430e-11;
    const scalar semimajorAxis1 = 3;
    const scalar semiminorAxis1 = 1;
    const scalar semimajorAxis2 = 2;
    const scalar semiminorAxis2 = 1;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;
    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};

    const std::size_t count = 500;
    const CountsForSteps counts = {count, count,
                                   count, count, count};
    const auto result = computeParametersForVelocity_(first, final, counts);
    const auto deltaV = computeVelocity(first, final, result, nu);
    const scalar realVelocity = std::sqrt(nu / semimajorAxis2);

    const scalar ellipticVelocity1 = std::sqrt(nu * 2 * semimajorAxis2 /
                                               semimajorAxis1 / (semimajorAxis1 + semimajorAxis2));
    const scalar ellipticVelocity2 = std::sqrt(nu * 2 * semimajorAxis1 /
                                               semimajorAxis2 / (semimajorAxis1 + semimajorAxis2));
    const scalar velocityOnFirstrbit  = std::sqrt(nu / semimajorAxis1);
    const scalar velocityOnSEcondOrbit = std::sqrt(nu / semimajorAxis2);
    const scalar firstdelta = ellipticVelocity1 - velocityOnFirstrbit;

    const scalar seconddelta = velocityOnSEcondOrbit - ellipticVelocity2;
    const scalar realResultDeltaV =  firstdelta + seconddelta;
    std::cout << deltaV  << " " << realResultDeltaV;
}



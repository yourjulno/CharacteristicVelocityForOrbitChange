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
    const scalar ascendNode2 = 30 * deg;
    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};

    const std::size_t count = 400;
    const CountsForSteps counts = {count, count,
                                   count};
    const scalar trashHold = 0.01;
    const auto result = computeParametersForVelocity_(first, final, counts, trashHold);
    std::cout << result.value().semimajorAxisTransferOrbit << " ";
    std::cout << result.value().periapseA1 << " ";
    std::cout << result.value().periapseA2 << " ";
    std::cout << result.value().angleBetweenImpulseAndTransversal1 * 180 / M_PI << " ";
    std::cout << result.value().functional << " ";


}

TEST(functional, FIRST) {
    const scalar nu = 3.986028e14;
    const scalar G = 6.67430e-11;
    const scalar semimajorAxis1 = 3;
    const scalar semiminorAxis1 = 1;
    const scalar semimajorAxis2 = 2;
    const scalar semiminorAxis2 = 1;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;
    const scalar e1 = 0.5;
    const scalar e2 = semiminorAxis2 / semimajorAxis2;
    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};

    const std::size_t count = 500;
    const CountsForSteps counts = {count, count,
                                   count};
    const scalar trashHold = 0.01;
    const auto result = computeParametersForVelocity_(first, final, counts, trashHold);
    const auto deltaV = computeVelocity(first, final, result);
    std::cout << result.value().semimajorAxisTransferOrbit << " ";
    //std::cout << result.value().periapseA1 << " ";
    std::cout << result.value().angleBetweenImpulseAndTransversal1 * 180 / M_PI << " ";
    std::cout << result.value().trueAnomaly1 * 180 / M_PI << " ";
    std::cout << result.value().trueAnomaly2 * 180 / M_PI << " ";

}

TEST(computeVelocity, FIRST){
    const scalar semimajorAxis1 = 3;
    const scalar semiminorAxis1 = 1;
    const scalar semimajorAxis2 = 2;
    const scalar semiminorAxis2 = 1;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;
    const scalar anomaly1 = 115.196 * deg;
    const scalar anomaly2 = -54.961 * deg;
    const scalar polarEq1 = semimajorAxis1 + semiminorAxis1 * std::cos(anomaly1 - ascendNode1);
    const scalar polarEq2 = semimajorAxis2 + semiminorAxis2 * std::cos(anomaly2 - ascendNode2);
    const scalar angle1 = 9.594 * deg;
    const scalar angle2 = 9.594 * deg;
    const scalar semimajorAxisTransferOrbit = 2.49981;
    const scalar periapse1 = semiminorAxis1 * std::sin(anomaly1 - ascendNode1) /
                             std::tan(angle1) / std::sqrt(semimajorAxis1);
    const scalar periapse2 = semiminorAxis2 * std::sin(anomaly2 - ascendNode2) /
                             std::tan(angle2) / std::sqrt(semimajorAxis2);
    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};
    const CorrectParametresForTransfer params = {angle1, angle2,
                                                 polarEq1, polarEq2,
                                                 periapse1, periapse2,
                                                 semimajorAxisTransferOrbit };
    const auto resultFunctional = getFunctionForMinimum(params);
    std::cout << resultFunctional;
}
TEST(Velocity, FIRST){
    const scalar nu = 3.986028e14;
    const scalar G =  6.67430e-11;
    const scalar semimajorAxis1 = 500e3;
    const scalar semiminorAxis1 = 500e3;
    const scalar semimajorAxis2 = 600e3;
    const scalar semiminorAxis2 = 600e3;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 0;

    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};

    const std::size_t count = 400;
    const CountsForSteps counts = {count, count,
                                   count};
    const scalar trashHold = 0.01;

    const scalar velocityOnFirstrbit  = std::sqrt(nu / semimajorAxis1);
    const scalar velocityOnSEcondOrbit = std::sqrt(nu / semimajorAxis2);
    const scalar firstdelta = std::sqrt(nu / semimajorAxis1) * (std::sqrt(2 * semimajorAxis2 /
                                                                                  (semimajorAxis2 + semimajorAxis1)) -1);


    const scalar seconddelta = std::sqrt(nu / semimajorAxis2) * ( - std::sqrt(2 * semimajorAxis1 /
                                                                           (semimajorAxis2 + semimajorAxis1)) + 1);
    const auto result = computeParametersForVelocity_(first, final, counts, trashHold);
    const auto deltaV = computeVelocity(first, final, result);
    const scalar realResultDeltaV =  firstdelta + seconddelta;
    std::cout << result.value().semimajorAxisTransferOrbit << " ";
//    std::cout << result.value().periapseA1 << " ";
//    std::cout << result.value().periapseA2 << " ";
    std::cout << result.value().angleBetweenImpulseAndTransversal1  << " ";
    //std::cout << result.value().functional << " ";
    std::cout << result.value().angleBetweenImpulseAndTransversal2 << " ";
    std::cout << deltaV   << " " << realResultDeltaV;
}



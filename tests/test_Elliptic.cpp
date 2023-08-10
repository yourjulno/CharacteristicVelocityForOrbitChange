//
// Created by julia on 10.08.23.
//
#include <gtest/gtest.h>
#include "../src/CoplanarCircularOrbit.hpp"
#include "../src/Elliptic.hpp"
using namespace Maneuvers;

constexpr scalar deg = M_PI / 180;
TEST(functional, elliptic) {
    const scalar nu = 3.986028e14;
    const scalar G = 6.67430e-11;
    const scalar semimajorAxis1 = 3;
    const scalar semiminorAxis1 = 1;
    const scalar semimajorAxis2 = 2;
    const scalar semiminorAxis2 = 1;
    const scalar ascendNode1 = 0;
    const scalar ascendNode2 = 30 * deg;
    const EllipticOrbit first = {semimajorAxis1, semiminorAxis1, ascendNode1};
    const EllipticOrbit final = {semimajorAxis2, semiminorAxis2, ascendNode2};

    const std::size_t count = 300;
    const CountsForSteps counts = {count, count,
                                   count};
    const scalar trashHold = 0.009;
    const auto result = computeParametersForVelocity(first, final, counts, trashHold);
    const auto deltaV = computeVelocity(first, final, result);
    std::cout << result.value().semimajorAxisTransferOrbit << " ";
    //std::cout << result.value().periapseA1 << " ";
    //std::cout << std::tan(result.value().impulseAngle1)  << " " << deltaV << " ";
    std:: cout << result.value().periapse1 << " ";
    std::cout << result.value().impulseAngle1 * 180 / M_PI << " ";
   std::cout << result.value().polarEq1 << " ";

}

TEST(Velocity, elliptic){
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
    const scalar firstdelta = std::sqrt(nu / semimajorAxis1) * (std::sqrt(2 * semimajorAxis2 /
                                                                          (semimajorAxis2 + semimajorAxis1)) -1);


    const scalar seconddelta = std::sqrt(nu / semimajorAxis2) * ( - std::sqrt(2 * semimajorAxis1 /
                                                                              (semimajorAxis2 + semimajorAxis1)) + 1);
    const auto result = computeParametersForVelocity(first, final, counts, trashHold);
    const auto deltaV = computeVelocity(first, final, result);
    const scalar realResultDeltaV =  firstdelta + seconddelta;
    std::cout << result.value().semimajorAxisTransferOrbit << std::endl
    << deltaV << std::endl << realResultDeltaV << std::endl << result.value().impulseAngle2 * 180 / M_PI << std::endl
    << result.value().polarEq1;
}
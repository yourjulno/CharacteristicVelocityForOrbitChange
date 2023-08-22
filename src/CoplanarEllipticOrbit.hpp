//
// Created by julia on 08.08.23.
//

#ifndef CHARACTERISTICVELOCITYFORORBITCHANGE_COPLANARELLIPTICORBIT_HPP
#define CHARACTERISTICVELOCITYFORORBITCHANGE_COPLANARELLIPTICORBIT_HPP


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

    struct CorrectParametresForTransfer {
        scalar angleBetweenImpulseAndTransversal1;
        scalar angleBetweenImpulseAndTransversal2;
        scalar trueAnomaly1;
        scalar trueAnomaly2;
        scalar semimajor1;
        scalar semimajor2;
        scalar polarEq1;
        scalar polarEq2;
        scalar periapseA1;
        scalar periapseA2;
        scalar semimajorAxisTransferOrbit;
        scalar semiminorAxisTransferOrbit;
        scalar functional;
    };


    bool isCorrectValues(const scalar cosDeltaTrueAnomaly1AscendNode, const scalar cosDeltaTrueAnomaly2AscendNode,
                         const scalar angleBetweenImpulseAndTransversalVelocity1,
                         const scalar angleBetweenImpulseAndTransversalVelocity2,
                         const scalar polarEq1, const scalar polarEq2) {
        const auto isDeltaCosCorrect = [](const scalar cosDeltaTrueAnomaly1AscendNode,
                                          const scalar cosDeltaTrueAnomaly2AscendNode) {
            return std::abs(cosDeltaTrueAnomaly1AscendNode - cosDeltaTrueAnomaly2AscendNode)
                   > std::numeric_limits<scalar>::epsilon();
        };

        const auto isPolarEqCorrect = [](const scalar polarEq1, const scalar polarEq2) {
            return std::abs(polarEq1 - polarEq2) > std::numeric_limits<scalar>::epsilon();
        };
        const auto isAnglesCorrect = [](const scalar angleBetweenImpulseAndTransversalVelocity) {
            return std::abs(angleBetweenImpulseAndTransversalVelocity - static_cast<scalar>(0)) >
                   std::numeric_limits<scalar>::epsilon() ||
                   std::abs(angleBetweenImpulseAndTransversalVelocity - M_PI) >
                   std::numeric_limits<scalar>::epsilon() ||
                   std::abs(angleBetweenImpulseAndTransversalVelocity - 2 * M_PI) >
                   std::numeric_limits<scalar>::epsilon();
        };

        return (isAnglesCorrect(angleBetweenImpulseAndTransversalVelocity1)
                && isAnglesCorrect(angleBetweenImpulseAndTransversalVelocity2))
               && (isDeltaCosCorrect(cosDeltaTrueAnomaly1AscendNode, cosDeltaTrueAnomaly2AscendNode) ||
                   isPolarEqCorrect(polarEq1, polarEq2));
    }

    struct TrigonometricFunctions {
        scalar cosDeltaTrueAnomaly1AscendNode;
        scalar cosDeltaTrueAnomaly2AscendNode;
        scalar sinDeltaTrueAnomaly1AscendNode;
        scalar sinDeltaTrueAnomaly2AscendNode;
        scalar sinDeltaTrueAnomaly1AscendNode1;
        scalar sinDeltaTrueAnomaly2AscendNode2;
        scalar semimajorAxis1;
        scalar semimajorAxis2;
    };

    TrigonometricFunctions computeCosDeltaTrueAnomaly(const EllipticOrbit &first, const EllipticOrbit &final,
                                                      const scalar trueAnomaly1, const scalar trueAnomaly2,
                                                      const scalar ascendNodeTransferOrbit) {
        const scalar cosDeltaTrueAnomaly1AscendNode = std::cos(trueAnomaly1 - ascendNodeTransferOrbit);
        const scalar cosDeltaTrueAnomaly2AscendNode = std::cos(trueAnomaly2 - ascendNodeTransferOrbit);
        const scalar sinDeltaTrueAnomaly1AscendNode = std::sin(trueAnomaly1 - ascendNodeTransferOrbit);
        const scalar sinDeltaTrueAnomaly2AscendNode = std::sin(trueAnomaly2 - ascendNodeTransferOrbit);
        const scalar sinDeltaTrueAnomaly1AscendNode1 = first.b * std::sin(trueAnomaly1 - first.w);
        const scalar sinDeltaTrueAnomaly2AscendNode2 = final.b * std::sin(trueAnomaly2 - final.w);
        return TrigonometricFunctions{cosDeltaTrueAnomaly1AscendNode, cosDeltaTrueAnomaly2AscendNode,
                                      sinDeltaTrueAnomaly1AscendNode, sinDeltaTrueAnomaly2AscendNode,
                                      sinDeltaTrueAnomaly1AscendNode1, sinDeltaTrueAnomaly2AscendNode2,
                                      first.a, final.a};
    }
    struct trigFunctions{
        scalar cosAnomaly1;
        scalar cosAnomaly2;
        scalar sinAnomaly1;
        scalar sinAnomaly2;
        scalar cosAscend;
        scalar sinAscend;
    };

    trigFunctions computeTrigFunctions(const scalar trueAnomaly1, const scalar trueAnomaly2,
                                       const scalar ascendNodeTransferOrbit){
        const scalar cosAnomaly1 = std::cos(trueAnomaly1);
        const scalar cosAnomaly2 = std::cos(trueAnomaly2);
        const scalar sinAnomaly1 = std::sin(trueAnomaly1);
        const scalar sinAnomaly2 = std::sin(trueAnomaly2);
        const scalar cosAscend = std::cos(ascendNodeTransferOrbit);
        const scalar sinAscend = std::sin(ascendNodeTransferOrbit);
        return trigFunctions{cosAnomaly1, cosAnomaly2, sinAnomaly1, sinAnomaly2,
                             cosAscend, sinAscend};
    }

    struct PolarEquations {
        scalar polarEq1;
        scalar polarEq2;
        scalar semimajor1;
        scalar semimajor2;
        scalar trueAnomaly1;
        scalar trueAnomaly2;
    };

    PolarEquations computeOrbitPolarEquations(const EllipticOrbit &first, const EllipticOrbit &final,
                                              const scalar trueAnomaly1, const scalar trueAnomaly2) {
        const scalar polarEq1 = first.a + first.b * std::cos(trueAnomaly1 - first.w);
        const scalar polarEq2 = final.a + final.b * std::cos(trueAnomaly2 - final.w);
        return PolarEquations{polarEq1, polarEq2, first.a, final.a,
                              trueAnomaly1, trueAnomaly2};
    }

    bool checkOnCorrectParams_(const CorrectParametresForTransfer &params) {

        const auto isParamsCorrect = [](const scalar parameter) -> bool {
            return !std::isinf(parameter) &&
                   !std::isnan(parameter);
        };
        const auto isAxisCorrect = [](const scalar semimajorAxis){
            return semimajorAxis > 0;
        };
        const auto isAngle1Correct = [](const scalar angle, const CorrectParametresForTransfer params){
            return params.polarEq1 * (1 / std::sqrt(params.semimajor1) - 1 / std::sqrt(params.semimajorAxisTransferOrbit)) /
            std::cos(angle) > 0;
        };
        const auto isAngle2Correct = [](const scalar angle, const CorrectParametresForTransfer params){
            return params.polarEq2 * ( - 1 / std::sqrt(params.semimajor2) + 1 / std::sqrt(params.semimajorAxisTransferOrbit)) /
                   std::cos(angle) > 0;
        };
        return isParamsCorrect(params.semimajorAxisTransferOrbit) &&
               isParamsCorrect(params.periapseA1) &&
               isParamsCorrect(params.periapseA2) &&
               isParamsCorrect(params.polarEq1) &&
               isParamsCorrect(params.polarEq2) &&
               isParamsCorrect(params.angleBetweenImpulseAndTransversal1) &&
               isParamsCorrect(params.angleBetweenImpulseAndTransversal2) &&
               isAxisCorrect(params.semimajorAxisTransferOrbit) &&
               isAxisCorrect(params.semiminorAxisTransferOrbit) &&
               isAngle1Correct(params.angleBetweenImpulseAndTransversal1, params) &&
               isAngle2Correct(params.angleBetweenImpulseAndTransversal2, params);
    }


    CorrectParametresForTransfer
    computeParamsForTransfer_(const PolarEquations &polarEq,
                              const trigFunctions &func_, const TrigonometricFunctions &func) {
        const scalar semiminorAxisTransferOrbit = (polarEq.polarEq1 - polarEq.polarEq2) /
                                                  (func_.cosAscend * (func_.cosAnomaly1 - func_.cosAnomaly1) +
                                                  func_.sinAscend * (func_.sinAnomaly2 - func_.sinAnomaly1));

        const scalar semimajorAxisTransferOrbit = polarEq.polarEq1 -
                semiminorAxisTransferOrbit * (func_.cosAnomaly1 * func_.cosAscend - func_.sinAnomaly1 * func_.sinAscend);
        //const scalar semiminorAxisTransferOrbit = polarEq.e * semimajorAxisTransferOrbit;

        const scalar sinDeltaTrueAnomaly1AscendNode =
                semiminorAxisTransferOrbit * (func_.sinAnomaly1 * func_.cosAscend - func_.cosAnomaly1 * func_.sinAscend);
        const scalar periapseA1 =
                polarEq.polarEq1 * (func.sinDeltaTrueAnomaly1AscendNode1 - sinDeltaTrueAnomaly1AscendNode) /
                (func.sinDeltaTrueAnomaly1AscendNode1 * std::sqrt(semimajorAxisTransferOrbit) -
                 sinDeltaTrueAnomaly1AscendNode * std::sqrt(func.semimajorAxis1));
        const scalar angleBetweenImpulseAndTransversal1 = std::atan(semiminorAxisTransferOrbit /
                                                                    (polarEq.polarEq1 -
                                                                     periapseA1 * std::sqrt(semimajorAxisTransferOrbit)));

        const scalar sinDeltaTrueAnomaly2AscendNode =
                semiminorAxisTransferOrbit * (func_.sinAnomaly2 * func_.cosAscend - func_.cosAnomaly2 * func_.sinAscend);
        const scalar periapseA2 =
                polarEq.polarEq2 * (func.sinDeltaTrueAnomaly2AscendNode2 - sinDeltaTrueAnomaly2AscendNode) /
                (func.sinDeltaTrueAnomaly2AscendNode2 * std::sqrt(semimajorAxisTransferOrbit) -
                 sinDeltaTrueAnomaly2AscendNode * std::sqrt(func.semimajorAxis2));
        const scalar angleBetweenImpulseAndTransversal2 = std::atan(sinDeltaTrueAnomaly2AscendNode /
                                                                    (polarEq.polarEq2 -
                                                                     periapseA2 * std::sqrt(semimajorAxisTransferOrbit)));

        return CorrectParametresForTransfer{angleBetweenImpulseAndTransversal1,
                                            angleBetweenImpulseAndTransversal2,
                                            polarEq.trueAnomaly1, polarEq.trueAnomaly2,
                                            polarEq.semimajor1, polarEq.semimajor2,
                                            polarEq.polarEq1, polarEq.polarEq2,
                                            periapseA1, periapseA2,
                                            semimajorAxisTransferOrbit, semiminorAxisTransferOrbit};

    }
    bool isAnglesCorrect(scalar impulseAngle1, scalar impulseAngle2, const CorrectParametresForTransfer params ){
        const auto isAngle1Correct = [](const scalar angle, const CorrectParametresForTransfer params){
            return params.polarEq1 * (1 / std::sqrt(params.semimajor1) - 1 / std::sqrt(params.semimajorAxisTransferOrbit)) /
                   std::cos(angle) > 0;
        };
        return isAngle1Correct(impulseAngle1, params) && isAngle1Correct(impulseAngle2, params);
    }


    scalar getFunctionForMinimum(const CorrectParametresForTransfer &paramsForTransfer) {

        const scalar cosAngleBetweenImpulseAndTransversal1 = std::cos(
                paramsForTransfer.angleBetweenImpulseAndTransversal1);
        const scalar sinAngleBetweenImpulseAndTransversal1 = std::sin(
                paramsForTransfer.angleBetweenImpulseAndTransversal1);
        const scalar cosAngleBetweenImpulseAndTransversal2 = std::cos(
                paramsForTransfer.angleBetweenImpulseAndTransversal2);
        const scalar sinAngleBetweenImpulseAndTransversal2 = std::sin(
                paramsForTransfer.angleBetweenImpulseAndTransversal2);
        const scalar sqrtSemimajorAxisTransferOrbit = std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit);
        const scalar H = ((paramsForTransfer.polarEq1 + paramsForTransfer.semimajorAxisTransferOrbit) /
                          paramsForTransfer.periapseA1 / sqrtSemimajorAxisTransferOrbit + 1) *
                         cosAngleBetweenImpulseAndTransversal1 -
                         ((paramsForTransfer.polarEq2 + paramsForTransfer.semimajorAxisTransferOrbit) /
                          paramsForTransfer.periapseA2 / sqrtSemimajorAxisTransferOrbit + 1) *
                         cosAngleBetweenImpulseAndTransversal2;

        const scalar G = ((paramsForTransfer.polarEq1 / paramsForTransfer.periapseA1 - paramsForTransfer.periapseA1) *
                          sinAngleBetweenImpulseAndTransversal1 -
                          (paramsForTransfer.polarEq2 / paramsForTransfer.periapseA2 - paramsForTransfer.periapseA2) *
                          sinAngleBetweenImpulseAndTransversal2) / sqrtSemimajorAxisTransferOrbit;

        const scalar F = ((paramsForTransfer.polarEq1 - paramsForTransfer.semimajorAxisTransferOrbit) *
                          (1 + sqrtSemimajorAxisTransferOrbit / paramsForTransfer.periapseA1) *
                          cosAngleBetweenImpulseAndTransversal1 +
                          (paramsForTransfer.polarEq1 -
                           paramsForTransfer.periapseA1 * sqrtSemimajorAxisTransferOrbit) *
                          sinAngleBetweenImpulseAndTransversal1 *
                          std::tan(paramsForTransfer.angleBetweenImpulseAndTransversal1) -
                          (paramsForTransfer.polarEq2 - paramsForTransfer.semimajorAxisTransferOrbit) *
                          (1 + sqrtSemimajorAxisTransferOrbit / paramsForTransfer.periapseA2) *
                          cosAngleBetweenImpulseAndTransversal2 +
                          (paramsForTransfer.polarEq2 -
                           paramsForTransfer.periapseA2 * sqrtSemimajorAxisTransferOrbit) *
                          sinAngleBetweenImpulseAndTransversal2 *
                          std::tan(paramsForTransfer.angleBetweenImpulseAndTransversal2))
                         / paramsForTransfer.semimajorAxisTransferOrbit;

        return std::sqrt(H * H + G * G + F * F);

    }
    scalar computeVelocity(const EllipticOrbit &first, const EllipticOrbit &final,
                           const std::optional<CorrectParametresForTransfer> &params) {

        const scalar sqrtSemimajorAxis = std::sqrt(params.value().semimajorAxisTransferOrbit);
        return  params.value().polarEq1 * (1 / sqrtSemimajorAxis -
                                          1 / std::sqrt(first.a)) /
               std::cos(M_PI) +
               params.value().polarEq2 * (-1 / sqrtSemimajorAxis +
                                          1 / std::sqrt(final.a)) /
               std::cos(0);
    }
    struct CountsForSteps {
        std::size_t anomaly1Count;
        std::size_t anomaly2Count;
        std::size_t ascendNodeCount;
    };

    std::optional<CorrectParametresForTransfer>
    computeParametersForVelocity_(const EllipticOrbit &first, const EllipticOrbit &final,
                                  const CountsForSteps &count, const scalar trashHold) {
        const scalar hTrueAnomaly1 = 2 * M_PI / count.anomaly1Count;
        const scalar hTrueAnomaly2 = 2 *  M_PI / count.anomaly2Count;
        const scalar hAscendNode = 2 * M_PI / count.ascendNodeCount;

        scalar minimumDeltaV = std::numeric_limits<scalar>::max();
        scalar minimumValue = std::numeric_limits<scalar>::max();
        std::optional<CorrectParametresForTransfer> params = std::nullopt;
        for (std::size_t hNode = 0; hNode < count.ascendNodeCount; hNode++){
            const scalar ascendingNode = hNode * hAscendNode;
            for (std::size_t hAnomaly2 = 0; hAnomaly2 < count.anomaly2Count; hAnomaly2++) {
                const scalar anomaly2 = hAnomaly2 * hTrueAnomaly2;
                for (std::size_t hAnomaly1 = 0; hAnomaly1 < count.anomaly1Count; hAnomaly1++)  {
                    const scalar anomaly1 = hAnomaly1* hTrueAnomaly1;
                    const auto trigonometricFunc = computeCosDeltaTrueAnomaly(first, final, anomaly1, anomaly2,
                                                                              ascendingNode);
                    const auto trigFunc = computeTrigFunctions(anomaly1, anomaly2, ascendingNode);
                    const auto polarEquations = computeOrbitPolarEquations(first, final, anomaly1, anomaly2);
                    const auto paramsForMinimum = computeParamsForTransfer_(polarEquations, trigFunc, trigonometricFunc);
                    const bool par = isCorrectValues(trigonometricFunc.cosDeltaTrueAnomaly1AscendNode, trigonometricFunc.cosDeltaTrueAnomaly2AscendNode,
                                                     paramsForMinimum.angleBetweenImpulseAndTransversal1, paramsForMinimum.angleBetweenImpulseAndTransversal2,
                                                     paramsForMinimum.polarEq1, paramsForMinimum.polarEq2);
                    const bool checkParams = checkOnCorrectParams_(paramsForMinimum);
                    if (par) {
                        const auto functionValue = getFunctionForMinimum(paramsForMinimum);
                        const auto velocityValue = computeVelocity(first, final, paramsForMinimum);
                        if (functionValue < minimumValue) {
                            minimumValue = functionValue;

                            //const auto velocityValue = computeVelocity(first, final, paramsForMinimum);
                            //if (velocityValue < minimumDeltaV) {
                            minimumDeltaV = velocityValue;
                            params = paramsForMinimum;
                            params->functional = minimumValue;
                            //}

                        }
                    }
                }
            }
        }


            return params;


    }



}


#endif //CHARACTERISTICVELOCITYFORORBITCHANGE_COPLANARELLIPTICORBIT_HPP

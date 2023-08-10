//
// Created by julia on 10.08.23.
//

#ifndef CHARACTERISTICVELOCITYFORORBITCHANGE_ELLIPTIC_HPP
#define CHARACTERISTICVELOCITYFORORBITCHANGE_ELLIPTIC_HPP

#include <cmath>
#include <array>
#include <optional>
#include "Types.h"

namespace Maneuvers {
    struct ParametersForTransfer{
        scalar semimajorAxisTransferOrbit;
        scalar semiminorAxisTransferOrbit;
        scalar polarEq1;
        scalar polarEq2;
        scalar periapse1;
        scalar periapse2;
        scalar impulseAngle1;
        scalar impulseAngle2;
    };

    struct PolarEquations{
        scalar polarEq1;
        scalar polarEq2;
    };

    PolarEquations computePolarEquations(const EllipticOrbit &first, const EllipticOrbit &final,
                                   const scalar trueAnomaly1, const scalar trueAnomaly2){
        const scalar polarEq1 = first.a + first.b * std::cos(trueAnomaly1 - first.w);
        const scalar polarEq2 = final.a + final.b * std::cos(trueAnomaly2 - final.w);
        return PolarEquations{polarEq1, polarEq2};
    }

    struct GeometryParamaters{
        scalar semimajorAxisTransferOrbit;
        scalar semiminorAxisTransferOrbit;
    };

    GeometryParamaters computeGeometryParameters(const PolarEquations &polarEq, const scalar trueAnomaly1,
                                                 const scalar trueAnomaly2, const scalar ascendNode){
        const scalar semiminorAxisTransferOrbit = (polarEq.polarEq1 - polarEq.polarEq2) /
                (std::cos(trueAnomaly1 - ascendNode) - std::cos(trueAnomaly2 - ascendNode));
        const scalar semimajorAxisTransferOrbit = polarEq.polarEq1 - semiminorAxisTransferOrbit *
                std::cos(trueAnomaly1 - ascendNode);
        return GeometryParamaters{semimajorAxisTransferOrbit, semiminorAxisTransferOrbit};
    }

    struct ImpulseAngles{
        scalar periapse1;
        scalar periapse2;
        scalar impulseAngle1;
        scalar impulseAngle2;
    };

    ImpulseAngles computeImpulseAngles(const EllipticOrbit &first, const EllipticOrbit &final,
                                       const scalar trueAnomaly1, const scalar trueAnomaly2,
                                       const scalar ascenNode, const GeometryParamaters &geomParams,
                                       const PolarEquations &polarEq){
        const scalar sinTrueAnomaly1AscendNode1 = first.b * std::sin(trueAnomaly1 - first.w);
        const scalar sinTrueAnomaly1AscendNode = geomParams.semiminorAxisTransferOrbit *
                std::sin(trueAnomaly1 - ascenNode);
        const scalar periapseA1 =
                polarEq.polarEq1 * (sinTrueAnomaly1AscendNode1 - sinTrueAnomaly1AscendNode) /
                (sinTrueAnomaly1AscendNode1 * std::sqrt(geomParams.semimajorAxisTransferOrbit) -
                 sinTrueAnomaly1AscendNode * std::sqrt(first.a));
        const scalar sinTrueAnomaly2AscendNode2 = first.b * std::sin(trueAnomaly2 - first.w);
        const scalar sinTrueAnomaly2AscendNode = geomParams.semiminorAxisTransferOrbit *
                                                 std::sin(trueAnomaly2 - ascenNode);
        const scalar periapseA2 =
                polarEq.polarEq2 * (sinTrueAnomaly2AscendNode2 - sinTrueAnomaly2AscendNode) /
                (sinTrueAnomaly2AscendNode2 * std::sqrt(geomParams.semimajorAxisTransferOrbit) -
                 sinTrueAnomaly2AscendNode * std::sqrt(final.a));

        const scalar impulseAngle1 = std::atan(sinTrueAnomaly1AscendNode1 /
                (polarEq.polarEq1 - periapseA1 * std::sqrt(first.a) ));
        const scalar impulseAngle2 = std::atan(sinTrueAnomaly2AscendNode2 /
                                               (polarEq.polarEq2 - periapseA2 * std::sqrt(final.a) ));
        return ImpulseAngles{periapseA1, periapseA2,
                             impulseAngle1, impulseAngle2};
    }

    bool isCorrectParams(const EllipticOrbit &first, const EllipticOrbit &final,
                         const PolarEquations &polar, const scalar trueAnomaly1, const scalar trueAnomaly2,
                         const scalar ascendNode, const ImpulseAngles &angles, const GeometryParamaters &geomParams){
        const auto isCosCorrect = [](const scalar trueAnomaly1, const scalar trueAnomaly2,
                const scalar ascendNode){
            return std::abs(std::cos(trueAnomaly1 - ascendNode) - std::cos(trueAnomaly2 - ascendNode)) >
            std::numeric_limits<scalar>::epsilon();
        };
//        const auto isVelocityPlus = [](const scalar polar, const scalar transferA,
//                const scalar A, const scalar impulseAngle){
//            return polar * (1 / std::sqrt(transferA) - 1 / std::sqrt(A)) / std::cos(impulseAngle) > 0;
//        };
        const auto isPolarEqCorrect = [](const PolarEquations polar){
            return std::abs(polar.polarEq1 - polar.polarEq2) > std::numeric_limits<scalar>::epsilon();
        };
        const auto isImpulseCorrect = [](const scalar impulseAngle){
            return std::abs(impulseAngle - M_PI / 2) >
                   std::numeric_limits<scalar>::epsilon() ||
                   std::abs(impulseAngle - 3 * M_PI / 2) >
                   std::numeric_limits<scalar>::epsilon();
        };
        const auto isAnglesCorrect = [](const scalar angleBetweenImpulseAndTransversalVelocity) {
            return std::abs(angleBetweenImpulseAndTransversalVelocity - static_cast<scalar>(0)) >
                   std::numeric_limits<scalar>::epsilon() ||
                   std::abs(angleBetweenImpulseAndTransversalVelocity - M_PI) >
                   std::numeric_limits<scalar>::epsilon() ||
                   std::abs(angleBetweenImpulseAndTransversalVelocity - 2 * M_PI) >
                   std::numeric_limits<scalar>::epsilon();
        };
        const auto isCorrectAxis = [](const scalar polar, const scalar periapse, const scalar semimajor){
            return std::abs(polar - periapse * std::sqrt(semimajor)) - static_cast<scalar>(0) >
            std::numeric_limits<scalar>::epsilon();
        };
        return isCosCorrect(trueAnomaly1, trueAnomaly2, ascendNode) &&
                isPolarEqCorrect(polar) &&
                isAnglesCorrect(angles.impulseAngle1) &&
                isAnglesCorrect(angles.impulseAngle2) &&
                isCorrectAxis(polar.polarEq1, angles.periapse1, angles.impulseAngle1) &&
                isCorrectAxis(polar.polarEq2, angles.periapse2, angles.impulseAngle2) &&
               // isVelocityPlus(polar.polarEq1, geomParams.semimajorAxisTransferOrbit, first.a, angles.impulseAngle1)&&
                //isVelocityPlus(polar.polarEq2,  final.a,geomParams.semimajorAxisTransferOrbit, angles.impulseAngle2) &&
                isImpulseCorrect(angles.impulseAngle1) &&
                isImpulseCorrect(angles.impulseAngle2);
    }
    scalar getFunctionForMinimum(const ParametersForTransfer &paramsForTransfer) {

        const scalar cosAngleBetweenImpulseAndTransversal1 = std::cos(
                paramsForTransfer.impulseAngle1);
        const scalar sinAngleBetweenImpulseAndTransversal1 = std::sin(
                paramsForTransfer.impulseAngle1);
        const scalar cosAngleBetweenImpulseAndTransversal2 = std::cos(
                paramsForTransfer.impulseAngle2);
        const scalar sinAngleBetweenImpulseAndTransversal2 = std::sin(
                paramsForTransfer.impulseAngle2);
        const scalar sqrtSemimajorAxisTransferOrbit = std::sqrt(paramsForTransfer.semimajorAxisTransferOrbit);
        const scalar H = ((paramsForTransfer.polarEq1 + paramsForTransfer.semimajorAxisTransferOrbit) /
                          paramsForTransfer.periapse1 / sqrtSemimajorAxisTransferOrbit + 1) *
                         cosAngleBetweenImpulseAndTransversal1 -
                         ((paramsForTransfer.polarEq2 + paramsForTransfer.semimajorAxisTransferOrbit) /
                          paramsForTransfer.periapse2 / sqrtSemimajorAxisTransferOrbit + 1) *
                         cosAngleBetweenImpulseAndTransversal2;

        const scalar G = ((paramsForTransfer.polarEq1 / paramsForTransfer.periapse1 - paramsForTransfer.periapse1) *
                          sinAngleBetweenImpulseAndTransversal1 -
                          (paramsForTransfer.polarEq2 / paramsForTransfer.periapse2 - paramsForTransfer.periapse2) *
                          sinAngleBetweenImpulseAndTransversal2) / sqrtSemimajorAxisTransferOrbit;

        const scalar F = ((paramsForTransfer.polarEq1 - paramsForTransfer.semimajorAxisTransferOrbit) *
                          (1 + sqrtSemimajorAxisTransferOrbit / paramsForTransfer.periapse1) *
                          cosAngleBetweenImpulseAndTransversal1 +
                          (paramsForTransfer.polarEq1 -
                           paramsForTransfer.periapse1 * sqrtSemimajorAxisTransferOrbit) *
                          sinAngleBetweenImpulseAndTransversal1 *
                          std::tan(paramsForTransfer.impulseAngle1) -
                          (paramsForTransfer.polarEq2 - paramsForTransfer.semimajorAxisTransferOrbit) *
                          (1 + sqrtSemimajorAxisTransferOrbit / paramsForTransfer.periapse2) *
                          cosAngleBetweenImpulseAndTransversal2 +
                          (paramsForTransfer.polarEq2 -
                           paramsForTransfer.periapse2 * sqrtSemimajorAxisTransferOrbit) *
                          sinAngleBetweenImpulseAndTransversal2 *
                          std::tan(paramsForTransfer.impulseAngle2))
                         / paramsForTransfer.semimajorAxisTransferOrbit;

        return std::sqrt(H * H + G * G + F * F);

    }
    scalar computeVelocity(const EllipticOrbit &first, const EllipticOrbit &final,
                           const std::optional<ParametersForTransfer> &params) {

        const scalar sqrtSemimajorAxis = std::sqrt(params.value().semimajorAxisTransferOrbit);
        return params.value().polarEq1 * (1 / sqrtSemimajorAxis -
                                          1 / std::sqrt(first.a)) /
               std::cos(params->impulseAngle1) +
               params.value().polarEq2 * (-1 / sqrtSemimajorAxis +
                                          1 / std::sqrt(final.a)) /
               std::cos(params->impulseAngle2);
    }

    struct CountsForSteps {
        std::size_t anomaly1Count;
        std::size_t anomaly2Count;
        std::size_t ascendNodeCount;
    };

    std::optional<ParametersForTransfer>
    computeParametersForVelocity(const EllipticOrbit &first, const EllipticOrbit &final,
                                  const CountsForSteps &count, const scalar trashHold) {
        const scalar hTrueAnomaly1 = 2 * M_PI / count.anomaly1Count;
        const scalar hTrueAnomaly2 = 2 * M_PI / count.anomaly2Count;
        const scalar hAscendNode = 2 * M_PI / count.ascendNodeCount;

        scalar minimumDeltaV = std::numeric_limits<scalar>::max();
        scalar minimumValue = std::numeric_limits<scalar>::max();
        std::optional<ParametersForTransfer> params = std::nullopt;
        for (std::size_t hNode = 0; hNode < count.ascendNodeCount; hNode++) {
            const scalar ascendingNode = hNode * hAscendNode;
            for (std::size_t hAnomaly2 = 0; hAnomaly2 < count.anomaly2Count; hAnomaly2++) {
                const scalar anomaly2 = hAnomaly2 * hTrueAnomaly2;
                for (std::size_t hAnomaly1 = 0; hAnomaly1 < count.anomaly1Count; hAnomaly1++) {
                    const scalar anomaly1 = hAnomaly1 * hTrueAnomaly1;
                    const auto polarEquations = computePolarEquations(first, final, anomaly1, anomaly2);
                    const auto geomParams = computeGeometryParameters(polarEquations, anomaly1, anomaly2,
                                                                      ascendingNode);
                    const auto impulseAgles = computeImpulseAngles(first, final, anomaly1, anomaly2,
                                                                   ascendingNode, geomParams, polarEquations);

                    const bool checkParams = isCorrectParams(first, final, polarEquations, anomaly1, anomaly2,
                                                             ascendingNode, impulseAgles, geomParams);
                    if (checkParams) {
                        const ParametersForTransfer paramsForMinimum = {geomParams.semimajorAxisTransferOrbit,
                                                                        geomParams.semiminorAxisTransferOrbit,
                                                                        polarEquations.polarEq1,
                                                                        polarEquations.polarEq2,
                                                                        impulseAgles.periapse1, impulseAgles.periapse2,
                                                                        impulseAgles.impulseAngle1,
                                                                        impulseAgles.impulseAngle2};
                        const auto functionValue = getFunctionForMinimum(paramsForMinimum);
                        if (functionValue < trashHold) {
                            minimumValue = functionValue;

                            const auto velocityValue = computeVelocity(first, final, paramsForMinimum);
                            if (velocityValue < minimumDeltaV) {
                               // minimumDeltaV = velocityValue;
                                params = paramsForMinimum;
                                //params->functional = minimumValue;
                                }

                            }

                        }
                    }
                }


            }

        return params;
    }


}
#endif //CHARACTERISTICVELOCITYFORORBITCHANGE_ELLIPTIC_HPP

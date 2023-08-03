//
// Created by julia on 17.07.23.
//

#ifndef VELOCITY_TYPES_H
#define VELOCITY_TYPES_H

namespace Maneuvers {
    using scalar = double;

    struct CircularOrbit {
    scalar a; // semimajor axis
    scalar e; // eccentricity
    scalar i; //inclination
    scalar ascendNode; // longitude of the ascending node
    scalar periapsisArg; // pericenter argument
    };

    struct EllipticOrbit{
        scalar a;
        scalar e; // eccentricity
        scalar b; //semiminorAxis
        scalar w;
    };
}  // namespace Maneuvers

#endif //VELOCITY_TYPES_H

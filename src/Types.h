//
// Created by julia on 17.07.23.
//

#ifndef VELOCITY_TYPES_H
#define VELOCITY_TYPES_H

namespace Maneuvers {
    using scalar = double;

    struct Orbit {
    scalar a; // semimajor axis
    scalar e; // eccentricity
    scalar i; //inclination
    scalar ascendNode; // longitude of the ascending node
    scalar periapsisArg; // pericenter argument
    };

    struct Orbit2{
        scalar a;
        scalar e; // eccentricity
        scalar b;
        scalar w;
    };
}  // namespace Maneuvers

#endif //VELOCITY_TYPES_H

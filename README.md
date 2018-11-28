# CoordSys

## Methods:
    * coordinateSystem - create a local coordinate system relative to
                         global coordinate system (0, 0, 0)
    * toLocal          - transforms a set of points (x, y, z) from global coordinate system
                         to local coordinate system
    * toGlobal         - transforms a set of points (x, y, z) from local coordinate system
                         to global coordinate system
## Properties:
    * origin           - a 1x3 array defining coordinate system offset (X, Y, Z)
                         from global coordinate systemLanguage
    * rotation         - a 3x3 transformation (rotation) matrix defining local coordinate system
                         axes relative to global coordinate system

## Notes:
    * When constructing a coordinate system using unit vectors (as oppose to Euler angles),
      CoordSys tests to see if these vectors are appropriate to define roataion matrix by:
      1) checking that they are unit vectors
      2) checking that the rotation matrix is orthogonal
      Since this process is time consuming, the user is advised to construct a coordinates
      system using Euler angles or removing the above tests by using the 'fast' method.
    * When using the 'fast' method to construct a local coordinate system using Euler angles,
      then the worst accuracy of the trigonometric approximation is 3e-3, but they are at least
      twice faster then javascript builtin standard trigonometric functions.


## Example:
An airplane is flying thru the air.
Its true air speed is 86 knots, with 4 degrees angle of attack and 2 degrees side slip angle.
its pitch is 12 degrees, bank angle is 35 degrees and it is headed on bearing 125 degrees.
what is its air speed in body and earth (inertial) coordinate system?

```javascript

// information
var airSpeed = 86,
    alpha = 4,
    beta = 2,
    pitch = 12,
    roll = 35,
    yaw = 125;

// define a coordinate system in which airplane body is global and measurements are local
var bodyWindCS = new coordinateSystem([0, 0, 0], alpha, 0, beta);

// airspeed in body coordinate system 
var airspeedBody = bodyWindCS.toGlobal([86, 0, 0]);

// define airplane body coordinate system relative to earth (inertial)
var bodyEarthCS = new coordinateSystem([0, 0, 0], pitch, roll ,yaw);

// airspeed in earth coordinate system
var airspeedEarth = bodyEarthCS.toGlobal(airspeedBody);
```

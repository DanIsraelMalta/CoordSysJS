/*
CoordSys - 2D/3D Cartesian Coordinate System Class

Methods:
    * coordinateSystem - create a local coordinate system relative to
                         global coordinate system (0, 0, 0)
    * toLocal          - transforms a set of points (x, y, z) from global coordinate system
                         to local coordinate system
    * toGlobal         - transforms a set of points (x, y, z) from local coordinate system
                         to global coordinate system
Properties:
    * origin           - a 1x3 array defining coordinate system offset (X, Y, Z)
                         from global coordinate systemLanguage
    * rotation         - a 3x3 transformation (rotation) matrix defining local coordinate system
                         axes relative to global coordinate system

Notes:
    * When constructing a coordinate system using unit vectors (as oppose to Euler angles),
      CoordSys tests to see if these vectors are appropriate to define roataion matrix by:
      1) checking that they are unit vectors
      2) checking that the rotation matrix is orthogonal
      Since this process is time consuming, the user is advised to construct a coordinates
      system using Euler angles or removing the above tests by using the 'fast' method.
    * When using the 'fast' method to construct a local coordinate system using Euler angles,
      then the worst accuracy of the trigonometric approximation is 3e-3, but they are at least
      twice faster then javascript builtin standard trigonometric functions.


Example:
   An airplane is flying thru the air.
   Its true air speed is 86 knots, with 4 degrees angle of attack and 2 degrees side slip angle.
   its pitch is 12 degrees, roll if 35 degrees and it is headed on bearing 125 degrees.
   what is its air speed in body and earth (inertial) coordinate system?

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


Dan I. Malta (malta.dan@gmail.com)

*/

// @constructor
// @param  {array}          coordinate system origin offset from global coordinate system (0,0,0)
//                          along [X, Y, Z] axis
// @param  {number | array} if an array then it is the coordinate system orthogonal "X" unit vector
//                          (which together with the other two vectors create an orthogonal base)
//                          if a number then it is the coordinate system pitch (elevation) angle
//                          in degrees
// @param  {number | array} if an array then it is the coordinate system orthogonal "Y" unit vector
//                          (which together with the other two vectors create an orthogonal base)
//                          if a number then it is the coordinate system roll angle
//                          in degrees
// @param  {number | array} if an array then it is the coordinate system orthogonal "Z" unit vector
//                          (which together with the other two vectors create an orthogonal base)
//                          if a number then it is the coordinate system psi (bearing / heading) angle
//                          in degrees
// @param  {string}         'fast'    - create coordinates system by "cutting corners":
//                                      * If a coordinate system is constructed using Euler angles,
//                                        then coordinateSystem use very accurate but fast trigonometric
//                                        approximations.
//                                      * If a coordinate system is constructed using unit vectors,
//                                        then coordinateSystem doesn't test unit vectors validity
//                                      This is very usefull when many low accuracy local coordinate
//                                      systems are needed - like for animation purposes.
//                          otherwise - create coordinates system using system trigonometric functions
// @return {}               coordinateSystem instance
function coordinateSystem(xi_offset, xi_U, xi_V, xi_W, xi_type) {

    // @param  {number} angle [degrees]
    // @return {number} sine of input argument
    // worse error is 3.0508965E-3
    var _sin = function(xi_angle) {
        if (xi_angle < 9) {
            return (0.017453292519943295 * xi_angle);
        } else if (xi_angle <= 180) {
            return ((4 * xi_angle * (180 - xi_angle)) / (40500 - xi_angle * (180 - xi_angle)));
        } else {
            xi_angle %= 180;
            return -((4 * xi_angle * (180 - xi_angle)) / (40500 - xi_angle * (180 - xi_angle)));
        }
    };

    // @param  {number} angle [degrees]
    // @return {number} sine of input argument
    // worse error is 3.0508965E-3
    var _cos = function(xi_angle) {
        // angle wrapping
        xi_angle = (xi_angle + 270) % 360;

        // output
        if (xi_angle < 9) {
            return (-0.017453292519943295 * xi_angle);
        } else if (xi_angle <= 180) {
            return -((4 * xi_angle * (180 - xi_angle)) / (40500 - xi_angle * (180 - xi_angle)));
        } else {
            xi_angle %= 180;
            return ((4 * xi_angle * (180 - xi_angle)) / (40500 - xi_angle * (180 - xi_angle)));
        }
    };

    // locals
    var c_deg2rad, _sTheta, _cTheta, _sPhi, _cPhi, _sPsi, _cPsi, _theta, _phi, _psi;

    // housekeeping
    xi_type = xi_type || 'SLOW';

    // coordinate system origin defined as an offset from global coordinate system (0, 0, 0)
    this.origin = [xi_offset[0], xi_offset[1], xi_offset[2]];

    // coordinate system axes definition (right system)
    if ((xi_U.length === 3) && (xi_V.length === 3) && (xi_W.length === 3)) {
        if (xi_type.toUpperCase() !== 'FAST') {
            // test to see if coordinate system is cartesian
            if ((Math.abs(1 - Math.sqrt(xi_U[0] * xi_U[0] + xi_U[1] * xi_U[1] + xi_U[2] * xi_U[2])) > 0.1) ||
                (Math.abs(1 - Math.sqrt(xi_V[0] * xi_V[0] + xi_V[1] * xi_V[1] + xi_V[2] * xi_V[2])) > 0.1) ||
                (Math.abs(1 - Math.sqrt(xi_W[0] * xi_W[0] + xi_W[1] * xi_W[1] + xi_W[2] * xi_W[2])) > 0.1)) {
                throw('coordinateSystem rotation vectors do not seem to be unit vectors.');
            }

            // test to see if rotation matrix is orthogonal
            if ((Math.abs(1 - (xi_U[0] * xi_U[0] + xi_V[0] * xi_V[0] + xi_W[0] * xi_W[0])) > 0.1) ||
                (Math.abs(1 - (xi_U[1] * xi_U[1] + xi_V[1] * xi_V[1] + xi_W[1] * xi_W[1])) > 0.1) ||
                (Math.abs(1 - (xi_U[2] * xi_U[2] + xi_V[2] * xi_V[2] + xi_W[2] * xi_W[2])) > 0.1) ||
                (xi_U[0] * xi_U[1] + xi_V[0] * xi_V[1] + xi_W[0] * xi_W[1] > 0.1) ||
                (xi_U[0] * xi_U[2] + xi_V[0] * xi_V[2] + xi_W[0] * xi_W[2] > 0.1) ||
                (xi_U[1] * xi_U[2] + xi_V[1] * xi_V[2] + xi_W[1] * xi_W[2] > 0.1)) {
                throw('coordinateSystem rotation unit vectors do not seem to create a valid transformation matrix.');
            }
        }

        // rotation matrix
        this.rotation = [[xi_U[0], xi_V[0], xi_W[0]],
                         [xi_U[1], xi_V[1], xi_W[1]],
                         [xi_U[2], xi_V[2], xi_W[2]]];
    } else if ((isFinite(xi_U)) && (isFinite(xi_V)) && (isFinite(xi_W))) {
        if (xi_type.toUpperCase() === 'FAST') {
            _sTheta = _sin(xi_U);
            _cTheta = _cos(xi_U);
            _sPhi = _sin(xi_V);
            _cPhi = _cos(xi_V);
            _sPsi = _sin(xi_W);
            _cPsi = _cos(xi_W);
        } else {
            // transform angles to radians
            c_deg2rad = Math.PI / 180;
            _theta = xi_U * c_deg2rad;
            _phi = xi_V * c_deg2rad;
            _psi = xi_W * c_deg2rad;
     
            // calculation
            _sTheta = Math.sin(_theta);
            _cTheta = Math.cos(_theta);
            _sPhi = Math.sin(_phi);
            _cPhi = Math.cos(_phi);
            _sPsi = Math.sin(_psi);
            _cPsi = Math.cos(_psi);
        }

        // rotation matrix orthogonal base from euler angles
        this.rotation = [[_cTheta * _cPsi, _cPhi * _sPsi + _sPhi * _sTheta * _cPsi,
                          _sPhi * _sPsi - _cPhi * _sTheta * _cPsi],
                         [-_cTheta * _sPsi, _cPhi * _cPsi - _sPhi * _sTheta * _sPsi,
                          _sPhi * _cPsi + _cPhi * _sTheta * _sPsi],
                         [_sTheta, -_sPhi * _cTheta, _cPhi * _cTheta]];
    } else {
        throw('coordinateSystem second, third and fourth input arguments must either be' +
              '1x3 vectors or finite numbers.');
    }
}

// @param  {2D array} array of points in the pattern of:
//                    [[point 1 X value, point 1 Y value, point 1 Z value]
//                     [point 2 X value, point 2 Y value, point 2 Z value]]
// @return {2D array} input points transformed from global coordinate system
//                    to local coordinate system
coordinateSystem.prototype.toLocal = function(xi_points) {

    // @private
    // @param  {number} value
    // @return {number} input value or 0
    var _zero = function(xi_value) {
        return (xi_value || 0);
    };

    // @param  {2D array} 3x3 matrix
    // @return {number}   determinant of input matrix
    var _determinant = function(xi_matrix) {
        return (xi_matrix[0][0] * xi_matrix[1][1] * xi_matrix[2][2] +
                xi_matrix[0][1] * xi_matrix[1][2] * xi_matrix[2][0] +
                xi_matrix[0][2] * xi_matrix[1][0] * xi_matrix[2][1] -
                xi_matrix[0][2] * xi_matrix[1][1] * xi_matrix[2][0] -
                xi_matrix[0][0] * xi_matrix[1][2] * xi_matrix[2][1] -
                xi_matrix[0][1] * xi_matrix[1][0] * xi_matrix[2][2]);
    };

    // @param  {2D array} 3x3 matrix
    // @param  {array} 1x3 array
    // @return {array} solution (X) of <input matrix> * X = <input array>
    var _solve = function(xi_matrix, xi_array) {
        var _detInv, _det = _determinant(xi_matrix), xo_out = [];

        if (Math.abs(_det) > 1e-10) {
            _detInv = 1 / _det;
            xo_out[0] = (xi_array[0] * xi_matrix[1][1] * xi_matrix[2][2] +
                         xi_array[2] * xi_matrix[0][1] * xi_matrix[1][2] +
                         xi_array[1] * xi_matrix[0][2] * xi_matrix[2][1] -
                         xi_array[2] * xi_matrix[0][2] * xi_matrix[1][1] -
                         xi_array[0] * xi_matrix[1][2] * xi_matrix[2][1] -
                         xi_array[1] * xi_matrix[0][1] * xi_matrix[2][2]) * _detInv;
            xo_out[1] = (xi_array[0] * xi_matrix[1][1] * xi_matrix[2][2] +
                         xi_array[2] * xi_matrix[0][1] * xi_matrix[1][2] +
                         xi_array[1] * xi_matrix[0][2] * xi_matrix[2][1] -
                         xi_array[2] * xi_matrix[0][2] * xi_matrix[1][1] -
                         xi_array[0] * xi_matrix[1][2] * xi_matrix[2][1] -
                         xi_array[1] * xi_matrix[0][1] * xi_matrix[2][2]) * _detInv;
            xo_out[2] = (xi_array[2] * xi_matrix[0][0] * xi_matrix[1][1] +
                         xi_array[1] * xi_matrix[0][1] * xi_matrix[2][0] +
                         xi_array[0] * xi_matrix[1][0] * xi_matrix[2][1] -
                         xi_array[0] * xi_matrix[1][1] * xi_matrix[2][0] -
                         xi_array[1] * xi_matrix[0][0] * xi_matrix[2][1] -
                         xi_array[2] * xi_matrix[0][1] * xi_matrix[1][0]) * _detInv;
        }
        return xo_out;
    };

    // locals
    var _detInv, _i, _point = [],
        _columns = _zero(xi_points[0].length),
        _rows = xi_points.length,
        xo_out = [];

    // iterate over all points
    if (_columns > 0) {
        for (_i = 0; _i < _rows; _i++) {

            // remove coordinate system offset from point
            _point = [_zero(xi_points[_i][0]) - this.origin[0],
                      _zero(xi_points[_i][1]) - this.origin[1],
                      _zero(xi_points[_i][2]) - this.origin[2]];

                      // transform to local
            xo_out.push(_solve(this.rotation, _point));
        }
    } else {
        // remove coordinate system offset from point
        _point = [_zero(xi_points[0]) - this.origin[0],
                 _zero(xi_points[1]) - this.origin[1],
                 _zero(xi_points[2]) - this.origin[2]];

        // transform to local
        xo_out = _solve(this.rotation, _point);
    }
    
    // output
    return xo_out;
};

// @param  {2D array} array of points in the pattern of:
//                    [[point 1 X value, point 1 Y value, point 1 Z value]
//                     [point 2 X value, point 2 Y value, point 2 Z value]]
// @return {array}    input points transformed from local coordinate system
//                    to global coordinate system
coordinateSystem.prototype.toGlobal = function(xi_points) {

    // @private
    // @param  {number} value
    // @return {number} input value or 0
    var _zero = function(xi_value) {
        return (xi_value || 0);
    };

    // @param  {2D array} 3x3 matrix
    // @param  {array} 1x3 array
    // @return {array} <1x3 input array> * <3x3 input matrix>
    var _matrixMul = function(xi_matrix, xi_array) {
        var _x = xi_array[0], _y = xi_array[1], _z = xi_array[2];
        return [_x * xi_matrix[0][0] + _y * xi_matrix[0][1] + _z * xi_matrix[0][2],
                _x * xi_matrix[1][0] + _y * xi_matrix[1][1] + _z * xi_matrix[1][2],
                _x * xi_matrix[2][0] + _y * xi_matrix[2][1] + _z * xi_matrix[2][2]];
    };

    // locals
    var _i, _point = [], _out = [],
        _columns = _zero(xi_points[0].length),
        _rows = xi_points.length,
        xo_out = [];

    if (_columns > 0) {
        for (_i = 0; _i < _rows; _i++) {
            _point = [_zero(xi_points[_i][0]),
                      _zero(xi_points[_i][1]),
                      _zero(xi_points[_i][2])];

            // transform point from local to global
            _out = _matrixMul(this.rotation, _point);

            // compensate coordinate system offset
            _out[0] = _out[0] + this.origin[0];
            _out[1] = _out[1] + this.origin[1];
            _out[2] = _out[2] + this.origin[2];
            xo_out.push(_out);
        }
    } else {
        _point = [_zero(xi_points[0]), _zero(xi_points[1]), _zero(xi_points[2])];

        // transform point from local to global
        xo_out = _matrixMul(this.rotation, _point);

        // compensate coordinate system offset
        xo_out[0] = xo_out[0] + this.origin[0];
        xo_out[1] = xo_out[1] + this.origin[1];
        xo_out[2] = xo_out[2] + this.origin[2];
    }

    // output
    return xo_out;
};

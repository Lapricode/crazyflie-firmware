/**
 * ,---------,       ____  _ __
 * |  ,-^-,  |      / __ )(_) /_______________ _____  ___
 * | (  O  ) |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * | / ,--´  |    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *    +------`   /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2019 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "app.h"
#include "FreeRTOS.h"
#include "task.h"

#include "controller.h"
#include "controller_pid.h"

#include "log.h"
#include "param.h"
#include "math.h"
#include "math3d.h"
#include "debug.h"
#include "physicalConstants.h"
#include "power_distribution.h"
#include "platform_defaults.h"
#include "stabilizer_types.h"

#define DEBUG_MODULE "MY_CONTROLLER"
#define NORM_FORCES_CONTROL 0
#define FORCE_TORQUE_CONTROL 1
#define CONTROL_SPEEDS 0
#define CONTROL_THRUSTS 1

static bool isInit = false;
static const double tol = 1e-10f;
static unsigned int debug_print_counter = 0; // a counter for debug printing rate

typedef struct mat33_s
{
  float m[3][3];
} mat33_t; // 3x3 matrix

typedef union vec3_u
{
  float v[3];
  struct
  {
    float x, y, z;
  };
} vec3_t; // 3x1 column vector

typedef struct vec4_u
{
  float v[4];
  struct
  {
    float v1, v2, v3, v4;
  };
} vec4_t; // 4x1 column vector

typedef struct vec12_u
{
  float v[12];
  struct
  {
    float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;
  };
} vec12_t; // 12x1 column vector

typedef union quat_u
{
  float q[4];
  struct
  {
    float w, x, y, z;
  };
} quat_t; // quaternion

typedef struct cf_state_s
{
  vec3_t rw;   // position w.r.t. world frame
  mat33_t Rwb; // body orientation in rotation matrix form w.r.t. world frame
  vec3_t vb;   // linear velocity w.r.t body frame
  vec3_t ob;   // angular velocity w.r.t body frame
} cf_state_t;  // crazyflie's state structure

// define some crazyflie model parameters
// Quadrotor system:
// state:      x = [rw, Rwb, vb, omegab] that has dimension 12 (from a Lie Theory point of view)
//             where:   rw \in R^(3x1) is the position in the world frame
//                      Rwb \in R^(3x3) is the rotation matrix of the body frame w.r.t. the world frame
//                      vb \in R^(3x1) is the linear velocity in the body frame
//                      omegab \in R^(3x1) is the angular velocity in the body frame
//                      qwb \in R^(4x1) is the orientation quaternion of the body frame w.r.t. the world frame
// control:    u = [u1, u2, u3, u4] \in R^(4x1)
//             where: ui is the angular speed of the i-th motor in rad/s
//                    ui is the thrust force produced by the i-th rotor in N
// rotors:     Fi = Kf * wi^2
//             Ti = Kt * wi^2
//             where:   Fi is the thrust force (in N) produced by the i-th rotor
//                      Ti is the torque (in Nm) produced by the i-th rotor
//                      wi is the angular speed (in rad/sec) produced by the i-th rotor
//             the motors are numbered in a clockwise manner, with motor 4 being in xy direction
//             the rotors 1, 3 rotate counter-clockwise, and the rotors 2, 4 rotate clockwise
//             the motor arms form right angles (90 degrees) with each other
// note:   we use the rotation matrix Rwb directly for the dynamics and the jacobians calculations
//         Rwb \in R^(3x3) is the rotation matrix of the body frame w.r.t. the world frame
//         Rwb \in R^(3x3) has dimension 3
static const float g = 9.81f;                              // gravity's acceleration (in m/sec^2)
static const float m_cf = 0.033f;                          // mass (in kg)
static const float l_cf = 0.046f;                          // arm length (in m)
static const float body_yaw0 = -3.0f / 4.0f * (float)M_PI; // assuming body_yaw0 is for the motor 1 at positive y direction, motor 2 at positive x direction and clockwise motor numbers
// static const mat33_t I_cf =
//     {{{16.6e-6f, 0.83e-6f, 0.72e-6f},
//       {0.83e-6f, 16.6e-6f, 1.8e-6f},
//       {0.72e-6f, 1.8e-6f, 29.3e-6f}}};
static const float kf = 2.25e-08f; // the coefficient parameter of the square model: thrust (N) vs. rotor_speed (rad/sec), for a single motor
static const float kt = 1.34e-10f; // the coefficient parameter of the square model: torque (N*m) vs. rotor_speed (rad/sec), for a single motor, kt = 0.00596 * kf

// define the LQR controller variables
static float Kinf[4][12] = {{0.0f}}; // initialize the LQR controller's Kinf matrix
static unsigned int Kinf_choice = 0; // parameter for the choice of the LQR controller's Kinf matrix
// the default Kinf constant LQR controller's matrix
// // for thrusts control
// static const float Kinf_default[4][12] = {
//     {-1.95735086e-02f, 1.56396351e-02f, 2.84411512e-01f, -8.72092347e-02f, -1.09199072e-01f, -1.11927317e-01f, -2.85184121e-02f, 2.27838048e-02f, 2.91126606e-01f, -1.60881220e-02f, -2.01556955e-02f, -1.12611490e-01f},
//     {1.83476168e-02f, 9.15738253e-03f, 2.84411512e-01f, -5.09942802e-02f, 1.02364745e-01f, 1.11001203e-01f, 2.67326348e-02f, 1.33356841e-02f, 2.91126606e-01f, -9.39330745e-03f, 1.88952170e-02f, 1.11678989e-01f},
//     {6.42730717e-03f, -1.03836440e-02f, 2.84411512e-01f, 5.78306154e-02f, 3.57187053e-02f, -1.08670124e-01f, 9.35488598e-03f, -1.51219963e-02f, 2.91126606e-01f, 1.06541458e-02f, 6.56469451e-03f, -1.09332019e-01f},
//     {-5.20141542e-03f, -1.44133736e-02f, 2.84411512e-01f, 8.03728995e-02f, -2.88843781e-02f, 1.09596238e-01f, -7.56910867e-03f, -2.09974926e-02f, 2.91126606e-01f, 1.48272836e-02f, -5.30421594e-03f, 1.10264521e-01f}};
// good
// static const float Kinf_default[4][12] = {
//     {-1.99420103e-02f, 1.57931844e-02f, 3.69101839e-01f, -8.80439216e-02f, -1.11202295e-01f, -1.17492620e-01f, -2.90516608e-02f, 2.30059996e-02f, 3.75420648e-01f, -1.62377284e-02f, -2.05147503e-02f, -1.18149355e-01f},
//     {1.87088023e-02f, 9.02202143e-03f, 3.69101839e-01f, -5.02583888e-02f, 1.04328215e-01f, 1.16487913e-01f, 2.72552978e-02f, 1.31398028e-02f, 3.75420648e-01f, -9.26140430e-03f, 1.92471478e-02f, 1.17138650e-01f},
//     {6.09419392e-03f, -1.02554339e-02f, 3.69101839e-01f, 5.71335790e-02f, 3.39077454e-02f, -1.13967399e-01f, 8.87283818e-03f, -1.49364616e-02f, 3.75420648e-01f, 1.05292058e-02f, 6.24009373e-03f, -1.14603183e-01f},
//     {-4.86098591e-03f, -1.45597719e-02f, 3.69101839e-01f, 8.11687313e-02f, -2.70336658e-02f, 1.14972106e-01f, -7.07647524e-03f, -2.12093408e-02f, 3.75420648e-01f, 1.49699269e-02f, -4.97249122e-03f, 1.15613889e-01f}};
// // for speeds control
// static const float Kinf_default[4][12] = {
//     {-4.34612169e+02f, 2.25644573e+02f, 1.07444819e+03f, -9.41655365e+02f, -1.20221510e+03f, -8.49498402e+02f, -3.56292177e+02f, 2.43187738e+02f, 7.49740794e+02f, -1.74737382e+02f, -1.98468534e+02f, -6.25310833e+02f},
//     {4.04150916e+02f, 1.71361278e+02f, 1.07444819e+03f, -7.09858035e+02f, 1.11832284e+03f, 8.47710105e+02f, 3.31365009e+02f, 1.84191075e+02f, 7.49740794e+02f, -1.30593138e+02f, 1.84689421e+02f, 6.23855890e+02f},
//     {2.77671758e+02f, -1.89076156e+02f, 1.07444819e+03f, 7.83778559e+02f, 7.52894884e+02f, -8.41753850e+02f, 2.25767115e+02f, -2.03282872e+02f, 7.49740794e+02f, 1.44306492e+02f, 1.21498330e+02f, -6.19160783e+02f},
//     {-2.47210505e+02f, -2.07929696e+02f, 1.07444819e+03f, 8.67734841e+02f, -6.69002624e+02f, 8.43542147e+02f, -2.00839947e+02f, -2.24095941e+02f, 7.49740794e+02f, 1.61024028e+02f, -1.07719217e+02f, 6.20615726e+02f}};
// very good Kinf!
static const float Kinf_default[4][12] = {
    {-8.46581487e+02f, 7.36733795e+02f, 4.72268962e+03f, -1.29534439e+03f, -1.49422694e+03f, -1.08286382e+03f, -5.19059177e+02f, 4.50885761e+02f, 1.08974537e+03f, -1.84844968e+02f, -2.14105733e+02f, -7.91565474e+02f},
    {7.87892956e+02f, 5.28998202e+02f, 4.72268962e+03f, -9.22677153e+02f, 1.39111275e+03f, 1.07797145e+03f, 4.83152862e+02f, 3.22528056e+02f, 1.08974537e+03f, -1.30529303e+02f, 1.99404240e+02f, 7.87900017e+02f},
    {4.72748481e+02f, -5.87769873e+02f, 4.72268962e+03f, 1.02593178e+03f, 8.18283109e+02f, -1.06465301e+03f, 2.87199386e+02f, -3.58484283e+02f, 1.08974537e+03f, 1.45250021e+02f, 1.14787654e+02f, -7.77962027e+02f},
    {-4.14059950e+02f, -6.77962125e+02f, 4.72268962e+03f, 1.19208976e+03f, -7.15168913e+02f, 1.06954537e+03f, -2.51293070e+02f, -4.14929533e+02f, 1.08974537e+03f, 1.70124250e+02f, -1.00086161e+02f, 7.81627484e+02f}};
// nice, very steady
// static const float Kinf_default[4][12] = {
//     {-2.05335885e+02f, 1.73245931e+02f, 1.52133651e+03f, -8.29196864e+02f, -9.84448360e+02f, -9.40618194e+02f, -2.14164522e+02f, 1.80558256e+02f, 7.32607133e+02f, -1.86034679e+02f, -2.21273750e+02f, -9.59859729e+02f},
//     {1.91544037e+02f, 1.16217978e+02f, 1.52133651e+03f, -5.54113732e+02f, 9.18470338e+02f, 9.34957957e+02f, 1.99793738e+02f, 1.20915368e+02f, 7.32607133e+02f, -1.23794995e+02f, 2.06479631e+02f, 9.54042342e+02f},
//     {9.70328200e+01f, -1.30024724e+02f, 1.52133651e+03f, 6.20161410e+02f, 4.60689738e+02f, -9.20175833e+02f, 1.00764693e+02f, -1.35301517e+02f, 7.32607133e+02f, 1.38604350e+02f, 1.02441943e+02f, -9.38864198e+02f},
//     {-8.32409726e+01f, -1.59439185e+02f, 1.52133651e+03f, 7.63149186e+02f, -3.94711716e+02f, 9.25836070e+02f, -8.63939084e+01f, -1.66172107e+02f, 7.32607133e+02f, 1.71225324e+02f, -8.76478235e+01f, 9.44681585e+02f}};
// first
// static const float Kinf_default[4][12] = {
//     {-1.59871592e+02f, 1.47942593e+02f, 4.85003523e+02f, -8.50985765e+02f, -9.24240986e+02f, -4.17322852e+02f, -2.36763983e+02f, 2.18805915e+02f, 5.76010628e+02f, -1.62353402e+02f, -1.77350962e+02f, -4.33812011e+02f},
//     {1.49004966e+02f, 1.20465728e+02f, 4.85003523e+02f, -6.86735740e+02f, 8.61592168e+02f, 4.18036045e+02f, 2.20682205e+02f, 1.77746246e+02f, 5.76010628e+02f, -1.29656803e+02f, 1.65374462e+02f, 4.34399595e+02f},
//     {1.18014214e+02f, -1.31355774e+02f, 4.85003523e+02f, 7.49514459e+02f, 6.68382092e+02f, -4.18925278e+02f, 1.73832117e+02f, -1.93862325e+02f, 5.76010628e+02f, 1.41657182e+02f, 1.25208356e+02f, -4.34970655e+02f},
//     {-1.07147588e+02f, -1.37052547e+02f, 4.85003523e+02f, 7.88207046e+02f, -6.05733273e+02f, 4.18212084e+02f, -1.57750340e+02f, -2.02689836e+02f, 5.76010628e+02f, 1.50353022e+02f, -1.13231857e+02f, 4.34383070e+02f}};

// static vec4_t hover_speeds = {{1900.0f, 1900.0f, 1900.0f, 1900.0f}}; // the angular speeds (in rad/sec) of the 4 rotors, for the crazyflie to hover
static float hover_adjust = 0.0f;                   // adjust hover speeds for hover calibration
static vec4_t control_speeds = {{0.0f}};            // the controlled angular speeds (in rad/sec) of the 4 rotors
static const float max_control_speed = 2618.0;      // the maximum angular speed (in rad/sec) of a rotor, approximately 25000.0f * (2.0f * (float)M_PI / 60.0f) = 2618.0f rad/sec
static vec4_t hover_speeds;                         // the angular speeds (in rad/sec) of the 4 rotors, for the crazyflie to hover
static vec4_t control_thrusts = {{0.0f}};           // the controlled thrusts (in N) of the 4 rotors
static const float max_control_thrust = 0.154f;     // the maximum thrust (in N) generated by only 1 motor, approximately kf * powf(max_control_speed, 2.0f)
static vec4_t hover_thrusts;                        // the thrusts (in N) of the 4 rotors, for the crazyflie to hover
static uint32_t update_rate = RATE_HL_COMMANDER;    // the update rate of the control loop (100 Hz default rate)
static uint32_t control_mode = NORM_FORCES_CONTROL; // or FORCE_TORQUE_CONTROL
static uint32_t control_inputs = CONTROL_THRUSTS;   // or CONTROL_SPEEDS

// for the normalized forces control mode (controlModeForce), control_mode = NORM_FORCES_CONTROL
static vec4_t control_norm_thrusts = {{0.0f}}; // the controlled normalized thrusts, in [0, 1], of the 4 rotors

// for the forces-torques control mode (controlModeForceTorque), control_mode = FORCE_TORQUE_CONTROL
static float control_thrust_total = 0.0f;      // the total thrust (in N) generated by the 4 rotors
static vec3_t control_body_torques = {{0.0f}}; // the body torques for each axis (x, y, z)

// functions definitions
float clamp_value(float, float, float);
quat_t qrpy_quat(const vec3_t);
mat33_t Rq_mat(const quat_t);
static mat33_t mat33_transpose(const mat33_t);
static vec3_t mat33_vec3_multiply(const mat33_t, const vec3_t);
static mat33_t mat33_mat33_multiply(const mat33_t, const mat33_t);
static vec4_t mat412_vec12_multiply(const float[4][12], const vec12_t);
static vec3_t SO3_minus_right(const mat33_t, const mat33_t);
static vec12_t compute_state_error(const cf_state_t, const cf_state_t);
void hover_control_init(void);
void Kinf_LQR_init(const unsigned int);

// clamp a float value betweeen two limit values
float clamp_value(float value, float min_lim, float max_lim)
{
  if (value < min_lim)
    return min_lim;
  if (value > max_lim)
    return max_lim;
  return value;
}

// convert roll, pitch, yaw angles to the corresponding quaternion
// there is also the function "struct quat rpy2quat(struct vec rpy)" of "math3d.h"
quat_t qrpy_quat(const vec3_t rpy) // rpy carries the roll, pitch, and yaw angles in radians
{
  quat_t q;
  float cr = cosf(0.5f * rpy.x);
  float sr = sinf(0.5f * rpy.x);
  float cp = cosf(0.5f * rpy.y);
  float sp = sinf(0.5f * rpy.y);
  float cy = cosf(0.5f * rpy.z);
  float sy = sinf(0.5f * rpy.z);

  q.w = cr * cp * cy + sr * sp * sy;
  q.x = sr * cp * cy - cr * sp * sy;
  q.y = cr * sp * cy + sr * cp * sy;
  q.z = cr * cp * sy - sr * sp * cy;

  return q;
}

// convert a quaternion to the corresponding rotation matrix
// there is also the function "struct mat33 quat2rotmat(struct quat q)"" of "math3d.h"
mat33_t Rq_mat(const quat_t q) // q is the quaternion
{
  mat33_t R;
  float w = q.w, x = q.x, y = q.y, z = q.z;

  R.m[0][0] = w * w + x * x - y * y - z * z;
  R.m[0][1] = 2.0f * (x * y - w * z);
  R.m[0][2] = 2.0f * (x * z + w * y);

  R.m[1][0] = 2.0f * (x * y + w * z);
  R.m[1][1] = w * w - x * x + y * y - z * z;
  R.m[1][2] = 2.0f * (y * z - w * x);

  R.m[2][0] = 2.0f * (x * z - w * y);
  R.m[2][1] = 2.0f * (y * z + w * x);
  R.m[2][2] = w * w - x * x - y * y + z * z;

  return R;
}

// compute the transpose of a 3x3 matrix
static mat33_t mat33_transpose(const mat33_t A)
{
  mat33_t At;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      At.m[i][j] = A.m[j][i];
  return At;
}

// compute the product A * b, where A is a 3x3 matrix and b is a 3x1 column vector
static vec3_t mat33_vec3_multiply(const mat33_t A, const vec3_t b)
{
  vec3_t result;
  for (int i = 0; i < 3; i++)
  {
    result.v[i] = A.m[i][0] * b.v[0] + A.m[i][1] * b.v[1] + A.m[i][2] * b.v[2];
  }
  return result;
}

// compute the 4x1 column vector product A * b, where A is a 4x12 matrix and b is a 12x1 column vector
static vec4_t mat412_vec12_multiply(const float A[4][12], const vec12_t b)
{
  vec4_t result;
  for (int i = 0; i < 4; i++)
  {
    float sum = 0.0f;
    for (int j = 0; j < 12; j++)
    {
      sum += A[i][j] * b.v[j];
    }
    result.v[i] = sum;
  }
  return result;
}

// compute the product A * B, where A, B are 3x3 matrices
static mat33_t mat33_mat33_multiply(const mat33_t A, const mat33_t B)
{
  mat33_t result;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      float sum = 0.0f;
      for (int k = 0; k < 3; k++)
      {
        sum += A.m[i][k] * B.m[k][j];
      }
      result.m[i][j] = sum;
    }
  return result;
}

// compute the right minus operation of the SO3 group
static vec3_t SO3_minus_right(const mat33_t R1, const mat33_t R2)
{
  vec3_t result = {{0.0f}};
  mat33_t R_rel = mat33_mat33_multiply(mat33_transpose(R2), R1); // this matrix goes inside the SO3 Log
  float tr = R_rel.m[0][0] + R_rel.m[1][1] + R_rel.m[2][2];
  float cos_theta = clamp_value((tr - 1.0f) / 2.0f, -1.0f, 1.0f); // clamp for numerical safety

  // do the SO3 Log operation
  float theta = acosf(cos_theta);

  // case 1: theta close to zero
  if (fabs(theta) < tol)
  {
    return result;
  }

  // case 2: theta close to pi (180 degrees)
  if (fabs((float)M_PI - theta) < tol)
  {
    float r00 = R_rel.m[0][0], r01 = R_rel.m[0][1], r02 = R_rel.m[0][2];
    float r10 = R_rel.m[1][0], r11 = R_rel.m[1][1], r12 = R_rel.m[1][2];
    float r20 = R_rel.m[2][0], r21 = R_rel.m[2][1], r22 = R_rel.m[2][2];

    float multiplier;
    if (!(fabs(r22 + 1.0f) < tol))
    {
      multiplier = theta / sqrtf(2.0f * (1.0f + r22));
      result.v[0] = multiplier * r02;
      result.v[1] = multiplier * r12;
      result.v[2] = multiplier * (1.0f + r22);
    }
    else if (!(fabs(r11 + 1.0f) < tol))
    {
      multiplier = theta / sqrtf(2.0f * (1.0f + r11));
      result.v[0] = multiplier * r01;
      result.v[1] = multiplier * (1.0f + r11);
      result.v[2] = multiplier * r21;
    }
    else if (!(fabs(r00 + 1.0f) < tol))
    {
      // fallback to first row
      multiplier = theta / sqrtf(2.0f * (1.0f + r00));
      result.v[0] = multiplier * (1.0f + r00);
      result.v[1] = multiplier * r10;
      result.v[2] = multiplier * r20;
    }
    return result;
  }

  // general case
  float scale = theta / (2.0f * sinf(theta));
  float tau_hat[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      tau_hat[i][j] = scale * (R_rel.m[i][j] - R_rel.m[j][i]);

  result.v[0] = tau_hat[2][1];
  result.v[1] = tau_hat[0][2];
  result.v[2] = tau_hat[1][0];

  return result;
}

// compute the state error (state_current - state_reference)
static vec12_t compute_state_error(const cf_state_t state_cur, const cf_state_t state_ref)
{
  // extract state components
  vec3_t rw_1 = state_cur.rw;
  mat33_t Rwb_1 = state_cur.Rwb;
  vec3_t vb_1 = state_cur.vb;
  vec3_t ob_1 = state_cur.ob;

  vec3_t rw_2 = state_ref.rw;
  mat33_t Rwb_2 = state_ref.Rwb;
  vec3_t vb_2 = state_ref.vb;
  vec3_t ob_2 = state_ref.ob;

  // compute errors
  vec3_t rw_error;
  for (int i = 0; i < 3; i++)
  {
    rw_error.v[i] = rw_1.v[i] - rw_2.v[i];
  }

  vec3_t Rwb_error;
  Rwb_error = SO3_minus_right(Rwb_1, Rwb_2);

  vec3_t vb_error;
  for (int i = 0; i < 3; i++)
  {
    vb_error.v[i] = vb_1.v[i] - vb_2.v[i];
  }

  vec3_t ob_error;
  for (int i = 0; i < 3; i++)
  {
    ob_error.v[i] = ob_1.v[i] - ob_2.v[i];
  }

  // concatenate into the error vector
  vec12_t state_error;
  for (int i = 0; i < 3; i++)
  {
    state_error.v[i] = rw_error.v[i];
    state_error.v[i + 3] = Rwb_error.v[i];
    state_error.v[i + 6] = vb_error.v[i];
    state_error.v[i + 9] = ob_error.v[i];
  }

  return state_error;
}

void hover_control_init(void)
{
  // compute the motors angular speeds (in rad/sec) needed, in order for the crazyflie to hover
  for (int i = 0; i < 4; i++)
  {
    hover_speeds.v[i] = sqrtf(m_cf * g / 4.0f / kf) + hover_adjust; // approximately 1900.0f rad/sec
    hover_thrusts.v[i] = m_cf * g / 4.0f + hover_adjust;            // approximately 0.081 N
  }
  return;
}

void Kinf_LQR_init(const unsigned int choice)
{
  if (choice == 0)
  {
    memcpy(Kinf, Kinf_default, sizeof(Kinf));
  }
  else if (choice == 1)
  {
    static const float Kinf_alt[4][12] = {
        {-4.19094067e+01f, 4.16796142e+01f, 4.97156665e+01f, -2.55682211e+02f, -2.59870002e+02f, -4.86336484e+01f, -6.32252245e+01f, 6.27184737e+01f, 1.10171880e+02f, -5.24415893e+01f, -5.43865958e+01f, -6.18570533e+01f},
        {4.12709645e+01f, 4.06818935e+01f, 4.97156665e+01f, -2.44279904e+02f, 2.54297833e+02f, 4.88008138e+01f, 6.21515659e+01f, 6.08947563e+01f, 1.10171880e+02f, -4.83435542e+01f, 5.29181909e+01f, 6.18886113e+01f},
        {4.09944100e+01f, -4.13254791e+01f, 4.97156665e+01f, 2.49880849e+02f, 2.45170915e+02f, -4.92034052e+01f, 6.13231763e+01f, -6.19759793e+01f, 1.10171880e+02f, 4.98171358e+01f, 4.77869881e+01f, -6.19525544e+01f},
        {-4.03559678e+01f, -4.10360285e+01f, 4.97156665e+01f, 2.50081266e+02f, -2.39598746e+02f, 4.90362398e+01f, -6.02495177e+01f, -6.16372506e+01f, 1.10171880e+02f, 5.09680078e+01f, -4.63185832e+01f, 6.19209964e+01f}};
    memcpy(Kinf, Kinf_alt, sizeof(Kinf));
  }
  else
  {
    memcpy(Kinf, Kinf_default, sizeof(Kinf));
  }
  return;
}

void appMain()
{
  DEBUG_PRINT("Waiting for activation of the controller and the estimator ...\n");

  while (1)
  {
    vTaskDelay(M2T(2000));

    // DEBUG_PRINT("My LQR Controller and My EKF Estimator!\n");
  }
}

void controllerOutOfTreeInit(void)
{
  // controllerPidInit();
  hover_control_init();
  Kinf_LQR_init(Kinf_choice);
  isInit = true;
  return;
}

bool controllerOutOfTreeTest(void)
{
  // return controllerPidTest();
  return isInit;
}

void controllerOutOfTree(control_t *control, const setpoint_t *setpoint,
                         const sensorData_t *sensors,
                         const state_t *state,
                         const stabilizerStep_t tick)
{
  // controllerPid(control, setpoint, sensors, state, tick);

  // This loop runs at approximately RATE_MAIN_LOOP = 1000 Hz rate.

  if (RATE_DO_EXECUTE(update_rate, tick)) // RATE_HL_COMMANDER is RATE_100_HZ = 100 Hz
  {
    // get the current state of the crazyflie
    vec3_t vw_cur = {{state->velocity.x, state->velocity.y, state->velocity.z}};
    quat_t qwb_cur = {{state->attitudeQuaternion.w, state->attitudeQuaternion.x, state->attitudeQuaternion.y, state->attitudeQuaternion.z}};
    mat33_t Rwb_cur = Rq_mat(qwb_cur);
    // vec3_t rpyb_cur = {{radians(state->attitude.roll), -radians(state->attitude.pitch), radians(state->attitude.yaw)}};
    vec3_t vb_cur = mat33_vec3_multiply(mat33_transpose(Rwb_cur), vw_cur);
    cf_state_t state_cur = {
        .rw = {{state->position.x, state->position.y, state->position.z}},
        .Rwb = Rwb_cur,
        .vb = vb_cur,
        .ob = {{radians(sensors->gyro.x), radians(sensors->gyro.y), radians(sensors->gyro.z)}},
    };

    // get the reference state of the crazyflie (in which state it should end up)
    vec3_t rpyb_ref = {{0., 0., radians(setpoint->attitude.yaw)}};
    // vec3_t rpyb_ref = {{radians(setpoint->attitude.roll), -radians(setpoint->attitude.pitch), radians(setpoint->attitude.yaw)}};
    quat_t qwb_ref = qrpy_quat(rpyb_ref);
    mat33_t Rwb_ref = Rq_mat(qwb_ref);
    cf_state_t state_ref = {
        .rw = {{setpoint->position.x, setpoint->position.y, setpoint->position.z}},
        .Rwb = Rwb_ref,
        .vb = {{0.0, 0.0, 0.0}},
        .ob = {{0.0, 0.0, 0.0}},
    };

    // compute the control vector signal
    vec12_t state_error = compute_state_error(state_cur, state_ref);
    vec4_t control_feedback = mat412_vec12_multiply(Kinf, state_error);
    for (int i = 0; i < 4; i++)
    {
      if (control_inputs == CONTROL_SPEEDS) // using the angular speeds of the rotors as control inputs
      {
        control_speeds.v[i] = hover_speeds.v[i] - control_feedback.v[i];
        control_speeds.v[i] = clamp_value(control_speeds.v[i], 0.0f, max_control_speed);
        control_thrusts.v[i] = kf * powf(control_speeds.v[i], 2.0f);
      }
      else if (control_inputs == CONTROL_THRUSTS) // using the thrusts provided by the rotors as control inputs
      {
        control_thrusts.v[i] = hover_thrusts.v[i] - control_feedback.v[i];
        control_thrusts.v[i] = clamp_value(control_thrusts.v[i], 0.0f, max_control_thrust);
      }
    }
    // if (RATE_DO_EXECUTE(1, debug_print_counter))
    //   DEBUG_PRINT("%lu: [%.3f, %.3f, %.3f, %.3f], [%.3f, %.3f, %.3f, %.3f]\n",
    //               tick,
    //               (double)hover_speeds.v1, (double)hover_speeds.v2, (double)hover_speeds.v3, (double)hover_speeds.v4,
    //               (double)control_feedback.v1, (double)control_feedback.v2, (double)control_feedback.v3, (double)control_feedback.v4);
  }

  // everything mentioned below is for a single motor
  // rotor_speed (rad/sec) vs. PWM:         omegar = sqrtf(8e-4f * PWM^2 + 53.33f * PWM)
  // PWM vs. rotor_speed (rad/sec):         PWM = -33333.0f + sqrtf(1250.0f * omegar^2 + 1111111111.0f)
  // thrust (N) vs. rotor_speed (rad/sec):  Fi = kf * omegar^2 = 2.25*e-8f * omegar^2
  // thrust (N) vs. PWM:                    Fi = 1.8e-11f * PWM^2 + 1.2e-6f * PWM
  // PWM vs. normal. thrust (N):            PWM = -33333.0f + sqrtf(8663836225.0f * norm_Fi + 1111111111.0f)
  // normal. PWM vs. normal. thrust (N):    norm_PWM = -0.50863f + sqrtf(0.25871f + 2.01727f * norm_Fi)

  if (setpoint->mode.z == modeDisable)
  {
    if (control_mode == NORM_FORCES_CONTROL) // using the normalized forces control mode
    {
      for (int i = 0; i < 4; i++)
      {
        control->normalizedForces[i] = 0.0f;
      }
      control->controlMode = controlModeForce;
    }
    else if (control_mode == FORCE_TORQUE_CONTROL) // using the forces-torques control mode
    {
      control->thrustSi = 0.0f;
      control->torqueX = 0.0f;
      control->torqueY = 0.0f;
      control->torqueZ = 0.0f;
      control->controlMode = controlModeForceTorque;
    }
  }
  else
  {
    if (control_mode == NORM_FORCES_CONTROL) // using the normalized forces control mode
    {
      for (int i = 0; i < 4; i++)
      {
        control_norm_thrusts.v[i] = clamp_value(control_thrusts.v[i] / max_control_thrust, 0.0f, 1.0f);
        control->normalizedForces[i] = control_norm_thrusts.v[i];
      }
      control->controlMode = controlModeForce;
    }
    else if (control_mode == FORCE_TORQUE_CONTROL) // using the forces-torques control mode
    {
      control_thrust_total = 0.0f;
      for (int i = 0; i < 4; i++)
      {
        control_thrust_total += control_thrusts.v[i];
      }
      const float cos_comp = l_cf * kf * cosf(body_yaw0);
      const float sin_comp = l_cf * kf * sinf(body_yaw0);
      control_body_torques.x = cos_comp * (powf(control_speeds.v1, 2.0f) - powf(control_speeds.v3, 2.0f)) - sin_comp * (powf(control_speeds.v4, 2.0f) - powf(control_speeds.v2, 2.0f));
      control_body_torques.y = sin_comp * (powf(control_speeds.v1, 2.0f) - powf(control_speeds.v3, 2.0f)) + cos_comp * (powf(control_speeds.v4, 2.0f) - powf(control_speeds.v2, 2.0f));
      control_body_torques.z = kt * (powf(control_speeds.v2, 2.0f) + powf(control_speeds.v4, 2.0f) - powf(control_speeds.v1, 2.0f) - powf(control_speeds.v3, 2.0f));
      control->thrustSi = control_thrust_total;
      control->torqueX = control_body_torques.x;
      control->torqueY = control_body_torques.y;
      control->torqueZ = control_body_torques.z;
      control->controlMode = controlModeForceTorque;
    }
  }

  // print some data for debugging
  if (RATE_DO_EXECUTE(1, debug_print_counter))
  {
    debug_print_counter = 0;
    DEBUG_PRINT("My LQR Controller is running!\n");
    // DEBUG_PRINT("Current state [%lu]: rw(m) = [%.3f, %.3f, %.3f],\t\t qwb = [%.3f, %.3f, %.3f, %.3f],\t\t rpyb(deg) = [%.3f, %.3f, %.3f],\n\t\t vb(m/s) = [%.3f, %.3f, %.3f],\t\t ob(deg/s) = [%.3f, %.3f, %.3f]\n",
    //             tick,
    //             (double)state_cur.rw.x, (double)state_cur.rw.y, (double)state_cur.rw.z,
    //             (double)state_cur.qwb.w, (double)state_cur.qwb.x, (double)state_cur.qwb.y, (double)state_cur.qwb.z,
    //             (double)degrees(rpyb_cur.x), (double)degrees(rpyb_cur.y), (double)degrees(rpyb_cur.z),
    //             (double)state_cur.vb.x, (double)state_cur.vb.y, (double)state_cur.vb.z,
    //             (double)degrees(state_cur.ob.x), (double)degrees(state_cur.ob.y), (double)degrees(state_cur.ob.z));
    // DEBUG_PRINT("Reference state [%lu]: rw(m) = [%.3f, %.3f, %.3f],\t\t qwb = [%.3f, %.3f, %.3f, %.3f],\t\t rpyb(deg) = [%.3f, %.3f, %.3f]\n",
    //             tick,
    //             (double)state_ref.rw.x, (double)state_ref.rw.y, (double)state_ref.rw.z,
    //             (double)state_ref.qwb.w, (double)state_ref.qwb.x, (double)state_ref.qwb.y, (double)state_ref.qwb.z,
    //             (double)degrees(rpyb_ref.x), (double)degrees(rpyb_ref.y), (double)degrees(rpyb_ref.z));
    // DEBUG_PRINT("State error [%lu]: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n",
    //             tick,
    //             (double)state_error[0], (double)state_error[1], (double)state_error[2],
    //             (double)state_error[3], (double)state_error[4], (double)state_error[5],
    //             (double)state_error[6], (double)state_error[7], (double)state_error[8],
    //             (double)state_error[9], (double)state_error[10], (double)state_error[11]);
    // DEBUG_PRINT("Motors speeds (rad/sec) [%lu]: m1 = %.0f,\t\t m2 = %.0f,\t\t m3 = %.0f,\t\t m4 = %.0f\n",
    //             tick,
    //             (double)control_speeds.v1, (double)control_speeds.v2, (double)control_speeds.v3, (double)control_speeds.v4);
    // DEBUG_PRINT("Motors thrusts (N) [%lu]: m1 = %.3f (%.3f),\t\t m2 = %.3f (%.3f),\t\t m3 = %.3f (%.3f),\t\t m4 = %.3f (%.3f)\n",
    //             tick,
    //             (double)control_thrusts.v1, (double)control_norm_thrusts.v1, (double)control_thrusts.v2, (double)control_norm_thrusts.v2,
    //             (double)control_thrusts.v3, (double)control_norm_thrusts.v3, (double)control_thrusts.v4, (double)control_norm_thrusts.v4);
    // DEBUG_PRINT("\n");
  }
  debug_print_counter += 1;
}

// /**
//  * Parameters for the LQR controller
//  */
// PARAM_GROUP_START(lqr_controller_params)
// /**
//  * @ brief Update rate (Hz)
//  */
// PARAM_ADD(PARAM_UINT32, update_rate, &update_rate)
// PARAM_GROUP_STOP(lqr_controller_params)

// /**
//  * Logging variables for the LQR controller
//  */
// LOG_GROUP_START(lqr_controller_logs)
// /**
//  * @brief Control normalized thrust for motor 4
//  */
// LOG_ADD(LOG_UINT32, update_rate, &update_rate)
// LOG_GROUP_STOP(lqr_controller_logs)
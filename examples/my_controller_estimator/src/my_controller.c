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

static bool isInit = false;
static const double tol = 1e-10f;

typedef struct mat_3_3_s
{
  float m[3][3];
} mat_3_3_t; // 3x3 matrix

typedef union vec_3_u
{
  float v[3];
  struct
  {
    float x, y, z;
  };
} vec_3_t; // 3x1 column vector

typedef struct vec_4_u
{
  float v[4];
  struct
  {
    float v1, v2, v3, v4;
  };
} vec_4_t; // 4x1 column vector

typedef struct vec_12_u
{
  float v[12];
  struct
  {
    float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;
  };
} vec_12_t; // 12x1 column vector

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
  vec_3_t rw; // position w.r.t. world frame
  quat_t qwb; // body orientation in quaternion form w.r.t. world frame
  vec_3_t vb; // linear velocity w.r.t body frame
  vec_3_t ob; // angular velocity w.r.t body frame
} cf_state_t; // crazyflie's state structure

// // define some logging variables
// static float posw_x;
// static float posw_y;
// static float posw_z;
// static float qwb_w;
// static float qwb_x;
// static float qwb_y;
// static float qwb_z;
// static float roll;
// static float pitch;
// static float yaw;
// static float velb_x;
// static float velb_y;
// static float velb_z;
// static float omegab_x;
// static float omegab_y;
// static float omegab_z;
// static float posw_x_ref;
// static float posw_y_ref;
// static float posw_z_ref;
// static float yaw_ref;
// static float posw_x_error;
// static float posw_y_error;
// static float posw_z_error;
// static float orient_x_error;
// static float orient_y_error;
// static float orient_z_error;
// static float velb_x_error;
// static float velb_y_error;
// static float velb_z_error;
// static float omegab_x_error;
// static float omegab_y_error;
// static float omegab_z_error;
// static float ctrl_speed_m1;
// static float ctrl_speed_m2;
// static float ctrl_speed_m3;
// static float ctrl_speed_m4;
// static float ctrl_thrust_m1;
// static float ctrl_thrust_m2;
// static float ctrl_thrust_m3;
// static float ctrl_thrust_m4;
// static float ctrl_norm_thrust_m1;
// static float ctrl_norm_thrust_m2;
// static float ctrl_norm_thrust_m3;
// static float ctrl_norm_thrust_m4;
static unsigned int debug_print_counter = 0; // a counter for debug printing rate

// define some crazyflie model parameters
// Quadrotor system:
// state:      x = [rw, qwb, vb, omegab] \in R^(13x1)
//             where   rw \in R^(3x1) is the position in the world frame
//                     qwb \in R^(4x1) is the orientation quaternion of the body frame w.r.t. the world frame
//                     vb \in R^(3x1) is the linear velocity in the body frame
//                     omegab \in R^(3x1) is the angular velocity in the body frame
// control:    u = [u1, u2, u3, u4] \in R^(4x1)
//             where u_i is the angular speed of the i-th motor in rad/s
// rotors:     Fi = Kf * ui^2
//             Ti = Kt * ui^2
//             where   Fi is the thrust force produced by the ith rotor
//                     Ti is the torque produced by the ith rotor
//             the motors are numbered in a clockwise manner, with motor 4 being in xy direction
//             the rotors 1, 3 rotate counter-clockwise, and the rotors 2, 4 rotate clockwise
//             the motor arms form right angles (90 degrees) with each other
// note:   the state x comes with the quaternion qwb of the rotation matrix Rwb,
//         but we use the rotation matrix Rwb directly for the dynamics and the jacobians calculations
//         Rwb \in R^(3x3) is the rotation matrix of the body frame w.r.t. the world frame
//         Rwb \in R^(3x3) has dimension 3
static const float g = 9.81f;                              // gravity's acceleration (in m/sec^2)
static const float m_cf = 0.033f;                          // mass (in kg)
static const float l = 0.046f;                             // arm length (in m)
static const float body_yaw0 = -3.0f / 4.0f * (float)M_PI; // assuming body_yaw0 is for the motor 1 at positive y direction, motor 2 at positive x direction and clockwise motor numbers
// static const mat_3_3_t CRAZYFLIE_INERTIA =
//     {{{16.6e-6f, 0.83e-6f, 0.72e-6f},
//       {0.83e-6f, 16.6e-6f, 1.8e-6f},
//       {0.72e-6f, 1.8e-6f, 29.3e-6f}}};
static const float kf = 2.25e-08f; // the coefficient parameter of the square model: thrust (N) vs. rotor_speed (rad/sec), for a single motor
static const float kt = 1.34e-10f; // the coefficient parameter of the square model: torque (N*m) vs. rotor_speed (rad/sec), for a single motor, kt = 0.00596 * kf

// define the LQR controller variables
static float Kinf[4][12] = {{0.0f}}; // initialize the LQR controller's Kinf matrix
static unsigned int Kinf_choice = 0; // parameter for the choice of the LQR controller's Kinf matrix
// the default Kinf constant LQR controller's matrix
static const float Kinf_default[4][12] = {
    {-4.34612169e+02f, 2.25644573e+02f, 1.07444819e+03f, -9.41655365e+02f, -1.20221510e+03f, -8.49498402e+02f, -3.56292177e+02f, 2.43187738e+02f, 7.49740794e+02f, -1.74737382e+02f, -1.98468534e+02f, -6.25310833e+02f},
    {4.04150916e+02f, 1.71361278e+02f, 1.07444819e+03f, -7.09858035e+02f, 1.11832284e+03f, 8.47710105e+02f, 3.31365009e+02f, 1.84191075e+02f, 7.49740794e+02f, -1.30593138e+02f, 1.84689421e+02f, 6.23855890e+02f},
    {2.77671758e+02f, -1.89076156e+02f, 1.07444819e+03f, 7.83778559e+02f, 7.52894884e+02f, -8.41753850e+02f, 2.25767115e+02f, -2.03282872e+02f, 7.49740794e+02f, 1.44306492e+02f, 1.21498330e+02f, -6.19160783e+02f},
    {-2.47210505e+02f, -2.07929696e+02f, 1.07444819e+03f, 8.67734841e+02f, -6.69002624e+02f, 8.43542147e+02f, -2.00839947e+02f, -2.24095941e+02f, 7.49740794e+02f, 1.61024028e+02f, -1.07719217e+02f, 6.20615726e+02f}};
// // very good Kinf!
// static const float Kinf_default[4][12] = {
//     {-8.46581487e+02f, 7.36733795e+02f, 4.72268962e+03f, -1.29534439e+03f, -1.49422694e+03f, -1.08286382e+03f, -5.19059177e+02f, 4.50885761e+02f, 1.08974537e+03f, -1.84844968e+02f, -2.14105733e+02f, -7.91565474e+02f},
//     {7.87892956e+02f, 5.28998202e+02f, 4.72268962e+03f, -9.22677153e+02f, 1.39111275e+03f, 1.07797145e+03f, 4.83152862e+02f, 3.22528056e+02f, 1.08974537e+03f, -1.30529303e+02f, 1.99404240e+02f, 7.87900017e+02f},
//     {4.72748481e+02f, -5.87769873e+02f, 4.72268962e+03f, 1.02593178e+03f, 8.18283109e+02f, -1.06465301e+03f, 2.87199386e+02f, -3.58484283e+02f, 1.08974537e+03f, 1.45250021e+02f, 1.14787654e+02f, -7.77962027e+02f},
//     {-4.14059950e+02f, -6.77962125e+02f, 4.72268962e+03f, 1.19208976e+03f, -7.15168913e+02f, 1.06954537e+03f, -2.51293070e+02f, -4.14929533e+02f, 1.08974537e+03f, 1.70124250e+02f, -1.00086161e+02f, 7.81627484e+02f}};
// static const float Kinf_default[4][12] = {
//     {-2.88838980e+02f, 1.72822432e+02f, 1.48239954e+03f, -9.05656797e+02f, -1.15819524e+03f, -7.87679074e+02f, -2.99871730e+02f, 2.17298829e+02f, 1.18444332e+03f, -1.86751780e+02f, -2.22895120e+02f, -9.57572347e+02f},
//     {2.69415982e+02f, 1.16048132e+02f, 1.48239954e+03f, -6.05576128e+02f, 1.08051434e+03f, 7.82918086e+02f, 2.79726551e+02f, 1.45706509e+02f, 1.18444332e+03f, -1.24307441e+02f, 2.07986837e+02f, 9.51750150e+02f},
//     {1.37083025e+02f, -1.29827257e+02f, 1.48239954e+03f, 6.77742354e+02f, 5.43364223e+02f, -7.70465812e+02f, 1.41701229e+02f, -1.63028394e+02f, 1.18444332e+03f, 1.39178931e+02f, 1.03287704e+02f, -9.36534774e+02f},
//     {-1.17660027e+02f, -1.59043307e+02f, 1.48239954e+03f, 8.33490571e+02f, -4.65683322e+02f, 7.75226801e+02f, -1.21556050e+02f, -1.99976943e+02f, 1.18444332e+03f, 1.71880290e+02f, -8.83794208e+01f, 9.42356971e+02f}};
// static const float Kinf_default[4][12] = {
//     {-1.92056058e+02f, 1.66843336e+02f, 1.06262736e+03f, -8.90113532e+02f, -1.02741107e+03f, -5.94822233e+02f, -2.50677705e+02f, 2.17581060e+02f, 9.44979595e+02f, -1.80833322e+02f, -2.09335144e+02f, -7.82149502e+02f},
//     {1.78782352e+02f, 1.19372305e+02f, 1.06262736e+03f, -6.33254219e+02f, 9.56631055e+02f, 5.92214129e+02f, 2.33370104e+02f, 1.55394019e+02f, 9.44979595e+02f, -1.27864513e+02f, 1.94964326e+02f, 7.78669504e+02f},
//     {1.06288602e+02f, -1.32665067e+02f, 1.06262736e+03f, 7.04133110e+02f, 5.60787294e+02f, -5.85078143e+02f, 1.38124477e+02f, -1.72726252e+02f, 9.44979595e+02f, 1.42254844e+02f, 1.12558588e+02f, -7.69171561e+02f},
//     {-9.30148967e+01f, -1.53550575e+02f, 1.06262736e+03f, 8.19234641e+02f, -4.90007277e+02f, 5.87686247e+02f, -1.20816876e+02f, -2.00248827e+02f, 9.44979595e+02f, 1.66442991e+02f, -9.81877706e+01f, 7.72651558e+02f}};
// static const float Kinf_default[4][12] = {
//     {-1.60612643e+02f, 1.48522601e+02f, 4.89398684e+02f, -7.43444853e+02f, -8.07643316e+02f, -2.95438739e+02f, -1.85683580e+02f, 1.71411262e+02f, 4.10460699e+02f, -1.61161356e+02f, -1.75953946e+02f, -4.29522740e+02f},
//     {1.49709730e+02f, 1.20792548e+02f, 4.89398684e+02f, -5.99722005e+02f, 7.52957260e+02f, 2.95984404e+02f, 1.73090277e+02f, 1.38981015e+02f, 4.10460699e+02f, -1.28838993e+02f, 1.64078225e+02f, 4.30205660e+02f},
//     {1.18211221e+02f, -1.31719136e+02f, 4.89398684e+02f, 6.54522773e+02f, 5.83423640e+02f, -2.96707597e+02f, 1.35709874e+02f, -1.51601300e+02f, 4.10460699e+02f, 1.40738802e+02f, 1.24494585e+02f, -4.31003456e+02f},
//     {-1.07308307e+02f, -1.37596013e+02f, 4.89398684e+02f, 6.88644086e+02f, -5.28737584e+02f, 2.96161931e+02f, -1.23116571e+02f, -1.58790976e+02f, 4.10460699e+02f, 1.49261547e+02f, -1.12618864e+02f, 4.30320536e+02f}};
// static const float Kinf_default[4][12] = {
//     {-3.54916121e+02f, 3.28975469e+02f, 1.07924184e+03f, -9.90321347e+02f, -1.07468498e+03f, -4.17564094e+02f, -3.22610413e+02f, 2.98392582e+02f, 6.69950553e+02f, -1.63787964e+02f, -1.78905688e+02f, -4.34054503e+02f},
//     {3.30720564e+02f, 2.68621595e+02f, 1.07924184e+03f, -8.00249593e+02f, 1.00165218e+03f, 4.18147131e+02f, 3.00641002e+02f, 2.42727533e+02f, 6.69950553e+02f, -1.30801132e+02f, 1.66799409e+02f, 4.34511257e+02f},
//     {2.63784855e+02f, -2.92868422e+02f, 1.07924184e+03f, 8.73430098e+02f, 7.79943711e+02f, -4.18736914e+02f, 2.37711881e+02f, -2.64742687e+02f, 6.69950553e+02f, 1.42930817e+02f, 1.26362912e+02f, -4.34781319e+02f},
//     {-2.39589298e+02f, -3.04728642e+02f, 1.07924184e+03f, 9.17140842e+02f, -7.06910907e+02f, 4.18153877e+02f, -2.15742469e+02f, -2.76377428e+02f, 6.69950553e+02f, 1.51658279e+02f, -1.14256633e+02f, 4.34324564e+02f}};
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

// static vec_4_t hover_speeds = {{1900.0f, 1900.0f, 1900.0f, 1900.0f}}; // the angular speeds (in rad/sec) of the 4 rotors, for the crazyflie to hover
static vec_4_t hover_speeds;                                              // the angular speeds (in rad/sec) of the 4 rotors, for the crazyflie to hover
static float hover_adjust = 0.0f;                                         // adjust hover speeds for hover calibration
static vec_4_t control_speeds = {{0.0f}};                                 // the controlled angular speeds (in rad/sec) of the 4 rotors
static float max_control_speed = 25000.0f * (2.0f * (float)M_PI / 60.0f); // the maximum angular speed (in rad/sec) of a rotor, approximately 2618.0f rad/sec
static vec_4_t control_thrusts = {{0.0f}};                                // the controlled thrusts (in N) of the 4 rotors
static unsigned int update_rate = RATE_HL_COMMANDER;                      // RATE_HL_COMMANDER;                      // the update rate of the control loop (100 Hz default rate)
static bool do_norm_forces_control = true;                                // if true, then controlModeForce, else controlModeForceTorque

// for the forces-torques control mode (controlModeForceTorque), do_norm_forces_control = false
static float control_thrust_total = 0.0f;       // the total thrust (in N) generated by the 4 rotors
static vec_3_t control_body_torques = {{0.0f}}; // the body torques for each axis (x, y, z)

// for the normalized forces control mode (controlModeForce), do_norm_forces_control = true
static const float max_thrust = 0.156f;         // the maximum thrust (in N) generated by only 1 motor
static vec_4_t control_norm_thrusts = {{0.0f}}; // the controlled normalized thrusts, in [0, 1], of the 4 rotors

// functions definitions
// float capAngle(float);
float clamp_to_unit_interval(float);
quat_t qrpy_quat(vec_3_t);
mat_3_3_t Rq_mat(quat_t);
static mat_3_3_t mat33_transpose(const mat_3_3_t);
static vec_3_t mat33_vec3_multiply(const mat_3_3_t, const vec_3_t);
static vec_4_t mat412_vec12_multiply(const float[4][12], const float[12]);
static vec_3_t SO3_minus_right(const mat_3_3_t, const mat_3_3_t);
static vec_12_t compute_state_error(const cf_state_t, const cf_state_t);
void hover_control_init(void);
void Kinf_LQR_init(const unsigned int);

// float capAngle(float angle)
// {
//   float result = angle;
//   while (result > 180.0f)
//   {
//     result -= 360.0f;
//   }
//   while (result < -180.0f)
//   {
//     result += 360.0f;
//   }
//   return result;
// }

// clamp a float value betweeen 0 and 1
float clamp_to_unit_interval(float value)
{
  if (value < 0.0f)
    return 0.0f;
  if (value > 1.0f)
    return 1.0f;
  return value;
}

// convert roll, pitch, yaw angles to the corresponding quaternion
// there is also the function "struct quat rpy2quat(struct vec rpy)" of "math3d.h"
quat_t qrpy_quat(vec_3_t rpy) // rpy carries the roll, pitch, and yaw angles in radians
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
mat_3_3_t Rq_mat(quat_t q) // q is the quaternion
{
  mat_3_3_t R;
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
static mat_3_3_t mat33_transpose(const mat_3_3_t A)
{
  mat_3_3_t At;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      At.m[i][j] = A.m[j][i];
  return At;
}

// compute the product A * b, where A is a 3x3 matrix and b is a 3x1 column vector
static vec_3_t mat33_vec3_multiply(const mat_3_3_t A, const vec_3_t b)
{
  vec_3_t result;
  for (int i = 0; i < 3; i++)
  {
    result.v[i] = A.m[i][0] * b.v[0] + A.m[i][1] * b.v[1] + A.m[i][2] * b.v[2];
  }
  return result;
}

// compute the 4x1 column vector product A * b, where A is a 4x12 matrix and b is a 12x1 column vector
static vec_4_t mat412_vec12_multiply(const float A[4][12], const float b[12])
{
  vec_4_t result;
  for (int i = 0; i < 4; i++)
  {
    float sum = 0.0f;
    for (int j = 0; j < 12; j++)
    {
      sum += A[i][j] * b[j];
    }
    result.v[i] = sum;
  }
  return result;
}

// // compute the product A * b, where A is a MxN matrix and b is a Nx1 column vector
// void mat_vec_multiply(const float **A, int M, int N, const float *b, float *result)
// {
//   for (int i = 0; i < M; i++)
//   {
//     float sum = 0.0f;
//     for (int j = 0; j < N; j++)
//     {
//       sum += A[i][j] * b[j];
//     }
//     result[i] = sum;
//   }
// }

// compute the product A * B, where A, B are 3x3 matrices
static mat_3_3_t mat33_mat33_multiply(const mat_3_3_t A, const mat_3_3_t B)
{
  mat_3_3_t result;
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
static vec_3_t SO3_minus_right(const mat_3_3_t R1, const mat_3_3_t R2)
{
  mat_3_3_t R_rel = mat33_mat33_multiply(mat33_transpose(R2), R1); // this matrix goes inside the SO3 Log
  vec_3_t result;
  float tr = R_rel.m[0][0] + R_rel.m[1][1] + R_rel.m[2][2];
  float cos_theta = (tr - 1.0f) / 2.0f;

  // clamp for numerical safety
  if (cos_theta > 1.0f)
    cos_theta = 1.0f;
  if (cos_theta < -1.0f)
    cos_theta = -1.0f;

  // do the SO3 Log operation
  float theta = acosf(cos_theta);

  // case 1: theta close to zero
  if (fabs(theta) < tol)
  {
    result.v[0] = 0.0f;
    result.v[1] = 0.0f;
    result.v[2] = 0.0f;
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
static vec_12_t compute_state_error(const cf_state_t state_cur, const cf_state_t state_ref)
{
  // extract state components
  vec_3_t rw_1 = state_cur.rw;
  quat_t qwb_1 = state_cur.qwb;
  mat_3_3_t Rwb_1 = Rq_mat(qwb_1);
  vec_3_t vb_1 = state_cur.vb;
  vec_3_t ob_1 = state_cur.ob;

  vec_3_t rw_2 = state_ref.rw;
  quat_t qwb_2 = state_ref.qwb;
  mat_3_3_t Rwb_2 = Rq_mat(qwb_2);
  vec_3_t vb_2 = state_ref.vb;
  vec_3_t ob_2 = state_ref.ob;

  // compute errors
  float rw_error[3];
  for (int i = 0; i < 3; i++)
  {
    rw_error[i] = rw_1.v[i] - rw_2.v[i];
  }

  vec_3_t Rwb_error;
  Rwb_error = SO3_minus_right(Rwb_1, Rwb_2);

  float vb_error[3];
  for (int i = 0; i < 3; i++)
  {
    vb_error[i] = vb_1.v[i] - vb_2.v[i];
  }

  float ob_error[3];
  for (int i = 0; i < 3; i++)
  {
    ob_error[i] = ob_1.v[i] - ob_2.v[i];
  }

  // concatenate into the error vector
  vec_12_t state_error;
  for (int i = 0; i < 3; i++)
  {
    state_error.v[i] = rw_error[i];
    state_error.v[i + 3] = Rwb_error.v[i];
    state_error.v[i + 6] = vb_error[i];
    state_error.v[i + 9] = ob_error[i];
  }

  return state_error;
}

void hover_control_init(void)
{
  // compute the motors angular speeds (in rad/sec) needed, in order for the crazyflie to hover
  for (int i = 0; i < 4; i++)
  {
    hover_speeds.v[i] = sqrtf(m_cf * g / 4.0f / kf) + hover_adjust; // approximately 1900.0f rad/sec
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
        {-1.59870249e+02f, 1.47942014e+02f, 4.85003521e+02f, -8.50979523e+02f, -9.24225961e+02f, -4.17274650e+02f, -2.36761509e+02f, 2.18804865e+02f, 5.76010628e+02f, -1.62351522e+02f, -1.77346369e+02f, -4.33753143e+02f},
        {1.49003459e+02f, 1.20466369e+02f, 4.85003521e+02f, -6.86742336e+02f, 8.61576234e+02f, 4.17984850e+02f, 2.20679491e+02f, 1.77747390e+02f, 5.76010628e+02f, -1.29658748e+02f, 1.65369703e+02f, 4.34337720e+02f},
        {1.18016121e+02f, -1.31356572e+02f, 4.85003521e+02f, 7.49521929e+02f, 6.68400251e+02f, -4.18866731e+02f, 1.73835419e+02f, -1.93863700e+02f, 5.76010628e+02f, 1.41659288e+02f, 1.25213524e+02f, -4.34901394e+02f},
        {-1.07149330e+02f, -1.37051810e+02f, 4.85003521e+02f, 7.88199929e+02f, -6.05750524e+02f, 4.18156531e+02f, -1.57753401e+02f, -2.02688554e+02f, 5.76010628e+02f, 1.50350982e+02f, -1.13236858e+02f, 4.34316817e+02f}};
    memcpy(Kinf, Kinf_alt, sizeof(Kinf));
  }
  else if (choice == 2)
  {
    static const float Kinf_alt[4][12] = {
        {-4.19094067e+01f, 4.16796142e+01f, 4.97156665e+01f, -2.55682211e+02f, -2.59870002e+02f, -4.86336484e+01f, -6.32252245e+01f, 6.27184737e+01f, 1.10171880e+02f, -5.24415893e+01f, -5.43865958e+01f, -6.18570533e+01f},
        {4.12709645e+01f, 4.06818935e+01f, 4.97156665e+01f, -2.44279904e+02f, 2.54297833e+02f, 4.88008138e+01f, 6.21515659e+01f, 6.08947563e+01f, 1.10171880e+02f, -4.83435542e+01f, 5.29181909e+01f, 6.18886113e+01f},
        {4.09944100e+01f, -4.13254791e+01f, 4.97156665e+01f, 2.49880849e+02f, 2.45170915e+02f, -4.92034052e+01f, 6.13231763e+01f, -6.19759793e+01f, 1.10171880e+02f, 4.98171358e+01f, 4.77869881e+01f, -6.19525544e+01f},
        {-4.03559678e+01f, -4.10360285e+01f, 4.97156665e+01f, 2.50081266e+02f, -2.39598746e+02f, 4.90362398e+01f, -6.02495177e+01f, -6.16372506e+01f, 1.10171880e+02f, 5.09680078e+01f, -4.63185832e+01f, 6.19209964e+01f}};
    memcpy(Kinf, Kinf_alt, sizeof(Kinf));
  }
  else if (choice == 3)
  {
    static const float Kinf_alt[4][12] = {
        {-9.91747496e+01f, 9.65265963e+01f, 1.56197707e+02f, -5.63413513e+02f, -5.83567245e+02f, -1.46499122e+02f, -1.47555659e+02f, 1.43317411e+02f, 2.34439227e+02f, -1.09300429e+02f, -1.14438427e+02f, -1.61232488e+02f},
        {9.51179962e+01f, 8.86224205e+01f, 1.56197707e+02f, -5.10118829e+02f, 5.58871360e+02f, 1.47929272e+02f, 1.41461742e+02f, 1.31110850e+02f, 2.34439227e+02f, -9.71517412e+01f, 1.09453672e+02f, 1.62535571e+02f},
        {8.94092949e+01f, -9.26930136e+01f, 1.56197707e+02f, 5.34890127e+02f, 5.11231853e+02f, -1.51221601e+02f, 1.32057053e+02f, -1.37224940e+02f, 2.34439227e+02f, 1.02150182e+02f, 9.63730239e+01f, -1.65509841e+02f},
        {-8.53525415e+01f, -9.24560032e+01f, 1.56197707e+02f, 5.38642215e+02f, -4.86535968e+02f, 1.49791450e+02f, -1.25963137e+02f, -1.37203320e+02f, 2.34439227e+02f, 1.04301989e+02f, -9.13882692e+01f, 1.64206758e+02f}};
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
  // This loop runs at approximately RATE_MAIN_LOOP = 1000 Hz rate.

  if (RATE_DO_EXECUTE(update_rate, tick)) // RATE_HL_COMMANDER is RATE_100_HZ = 100 Hz
  {
    // get the current state of the crazyflie
    vec_3_t vw_cur = {{state->velocity.x, state->velocity.y, state->velocity.z}};
    quat_t qwb_cur = {{state->attitudeQuaternion.w, state->attitudeQuaternion.x, state->attitudeQuaternion.y, state->attitudeQuaternion.z}};
    mat_3_3_t Rwb_cur = Rq_mat(qwb_cur);
    // vec_3_t rpyb_cur = {{radians(state->attitude.roll), -radians(state->attitude.pitch), radians(state->attitude.yaw)}};
    vec_3_t vb_cur = mat33_vec3_multiply(mat33_transpose(Rwb_cur), vw_cur);
    cf_state_t state_cur = {
        .rw = {{state->position.x, state->position.y, state->position.z}},
        .qwb = {{qwb_cur.w, qwb_cur.x, qwb_cur.y, qwb_cur.z}},
        .vb = {{vb_cur.x, vb_cur.y, vb_cur.z}},
        .ob = {{radians(sensors->gyro.x), radians(sensors->gyro.y), radians(sensors->gyro.z)}},
    };

    // get the reference state of the crazyflie (in which state it should end up)
    vec_3_t rpyb_ref = {{0., 0., radians(setpoint->attitude.yaw)}};
    // vec_3_t rpyb_ref = {{radians(setpoint->attitude.roll), -radians(setpoint->attitude.pitch), radians(setpoint->attitude.yaw)}};
    quat_t qwb_ref = qrpy_quat(rpyb_ref);
    cf_state_t state_ref = {
        .rw = {{setpoint->position.x, setpoint->position.y, setpoint->position.z}},
        .qwb = {{qwb_ref.w, qwb_ref.x, qwb_ref.y, qwb_ref.z}},
        .vb = {{0.0, 0.0, 0.0}},
        .ob = {{0.0, 0.0, 0.0}},
    };

    // compute the control vector signal
    vec_12_t state_error = compute_state_error(state_cur, state_ref);
    vec_4_t control_feedback = mat412_vec12_multiply(Kinf, state_error.v);
    for (int i = 0; i < 4; i++)
    {
      control_speeds.v[i] = hover_speeds.v[i] - control_feedback.v[i];
      if (control_speeds.v[i] < 0.0f)
        control_speeds.v[i] = 0.0f;
      if (control_speeds.v[i] > max_control_speed)
        control_speeds.v[i] = max_control_speed;
      control_thrusts.v[i] = kf * powf(control_speeds.v[i], 2.0f);
    }
    // if (RATE_DO_EXECUTE(1, debug_print_counter))
    //   DEBUG_PRINT("%lu: [%.3f, %.3f, %.3f, %.3f], [%.3f, %.3f, %.3f, %.3f]\n",
    //               tick,
    //               (double)hover_speeds.v1, (double)hover_speeds.v2, (double)hover_speeds.v3, (double)hover_speeds.v4,
    //               (double)control_feedback.v1, (double)control_feedback.v2, (double)control_feedback.v3, (double)control_feedback.v4);
  }

  // everything mentioned below is for a single motor
  // rotor_speed (rad/sec) vs. PWM:         omegar = sqrt(8e-4f * PWM^2 + 53.33f * PWM)
  // PWM vs. rotor_speed (rad/sec):         PWM = -33333.0f + sqrtf(1250.0f * omegar^2 + 1111111111.0f)
  // thrust (N) vs. rotor_speed (rad/sec):  Fi = kf * ui^2 = 2.25*e-8f * ui^2
  // thrust (N) vs. PWM:                    Fi = 1.8e-11f * PWM^2 + 1.2e-6f * PWM
  // PWM vs. normal. thrust (N):            PWM = -33333.0f + sqrtf(8663836225.0f * norm_Fi + 1111111111.0f)
  // normal. PWM vs. normal. thrust (N):    norm_PWM = -0.50863f + sqrtf(0.25871f + 2.01727f * norm_Fi)

  if (setpoint->mode.z == modeDisable)
  {
    if (do_norm_forces_control) // using the normalized forces control mode
    {
      for (int i = 0; i < 4; i++)
      {
        control->normalizedForces[i] = 0.0f;
      }
      control->controlMode = controlModeForce;
    }
    else // using the forces-torques control mode
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
    if (do_norm_forces_control) // using the normalized forces control mode
    {
      for (int i = 0; i < 4; i++)
      {
        control_norm_thrusts.v[i] = clamp_to_unit_interval(control_thrusts.v[i] / max_thrust);
        control->normalizedForces[i] = control_norm_thrusts.v[i];
      }
      control->controlMode = controlModeForce;
    }
    else // using the forces-torques control mode
    {
      control_thrust_total = 0.0f;
      for (int i = 0; i < 4; i++)
      {
        control_thrust_total += control_thrusts.v[i];
      }
      const float cos_comp = l * kf * cosf(body_yaw0);
      const float sin_comp = l * kf * sinf(body_yaw0);
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
    // DEBUG_PRINT("Kinf [%lu]: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n"
    //             "             [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n"
    //             "             [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n"
    //             "             [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]",
    //             tick,
    //             (double)Kinf[0][0], (double)Kinf[0][1], (double)Kinf[0][2], (double)Kinf[0][3],
    //             (double)Kinf[0][4], (double)Kinf[0][5], (double)Kinf[0][6], (double)Kinf[0][7],
    //             (double)Kinf[0][8], (double)Kinf[0][9], (double)Kinf[0][10], (double)Kinf[0][11],
    //             (double)Kinf[1][0], (double)Kinf[1][1], (double)Kinf[1][2], (double)Kinf[1][3],
    //             (double)Kinf[1][4], (double)Kinf[1][5], (double)Kinf[1][6], (double)Kinf[1][7],
    //             (double)Kinf[1][8], (double)Kinf[1][9], (double)Kinf[1][10], (double)Kinf[1][11],
    //             (double)Kinf[2][0], (double)Kinf[2][1], (double)Kinf[2][2], (double)Kinf[2][3],
    //             (double)Kinf[2][4], (double)Kinf[2][5], (double)Kinf[2][6], (double)Kinf[2][7],
    //             (double)Kinf[2][8], (double)Kinf[2][9], (double)Kinf[2][10], (double)Kinf[2][11],
    //             (double)Kinf[3][0], (double)Kinf[3][1], (double)Kinf[3][2], (double)Kinf[3][3],
    //             (double)Kinf[3][4], (double)Kinf[3][5], (double)Kinf[3][6], (double)Kinf[3][7],
    //             (double)Kinf[3][8], (double)Kinf[3][9], (double)Kinf[3][10], (double)Kinf[3][11]);
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

  // // store data to the logging variables
  // posw_x = state_cur.rw.x;
  // posw_y = state_cur.rw.y;
  // posw_z = state_cur.rw.z;
  // qwb_w = state_cur.qwb.w;
  // qwb_x = state_cur.qwb.x;
  // qwb_y = state_cur.qwb.y;
  // qwb_z = state_cur.qwb.z;
  // roll = rpyb_cur.x;
  // pitch = rpyb_cur.y;
  // yaw = rpyb_cur.z;
  // velb_x = state_cur.vb.x;
  // velb_y = state_cur.vb.y;
  // velb_z = state_cur.vb.z;
  // omegab_x = state_cur.ob.x;
  // omegab_y = state_cur.ob.y;
  // omegab_z = state_cur.ob.z;
  // posw_x_ref = state_ref.rw.x;
  // posw_y_ref = state_ref.rw.y;
  // posw_z_ref = state_ref.rw.z;
  // yaw_ref = rpyb_ref.z;
  // posw_x_error = state_error[0];
  // posw_y_error = state_error[1];
  // posw_z_error = state_error[2];
  // orient_x_error = state_error[3];
  // orient_y_error = state_error[4];
  // orient_z_error = state_error[5];
  // velb_x_error = state_error[6];
  // velb_y_error = state_error[7];
  // velb_z_error = state_error[8];
  // omegab_x_error = state_error[9];
  // omegab_y_error = state_error[10];
  // omegab_z_error = state_error[11];
  // ctrl_speed_m1 = control_speeds.v1;
  // ctrl_speed_m2 = control_speeds.v2;
  // ctrl_speed_m3 = control_speeds.v3;
  // ctrl_speed_m4 = control_speeds.v4;
  // ctrl_thrust_m1 = control_thrusts.v1;
  // ctrl_thrust_m2 = control_thrusts.v2;
  // ctrl_thrust_m3 = control_thrusts.v3;
  // ctrl_thrust_m4 = control_thrusts.v4;
  // ctrl_norm_thrust_m1 = control_norm_thrusts.v1;
  // ctrl_norm_thrust_m2 = control_norm_thrusts.v2;
  // ctrl_norm_thrust_m3 = control_norm_thrusts.v3;
  // ctrl_norm_thrust_m4 = control_norm_thrusts.v4;

  // controllerPid(control, setpoint, sensors, state, tick);
}

PARAM_GROUP_START(lqr_controller_params)
// /**
//  * @ brief Update rate (Hz)
//  */
// PARAM_ADD(PARAM_UINT32, update_rate, &update_rate)
// /**
//  * @ brief Choice of LQR controller's K infinity matrix
//  */
// PARAM_ADD(PARAM_UINT32, Kinf_choice, &Kinf_choice)
// /**
//  * @ brief Additional hover speed (rad/sec) for each motor, for adjustment purposes
//  */
// PARAM_ADD(PARAM_FLOAT, hover_adjust, &hover_adjust)
PARAM_GROUP_STOP(lqr_controller_params)

/**
 * Logging variables for the LQR controller
 */
LOG_GROUP_START(lqr_controller_logs)
// /**
//  * @brief Position x (m) - world frame
//  */
// LOG_ADD(LOG_FLOAT, posw_x, &posw_x)
// /**
//  * @brief Position y (m) - world frame
//  */
// LOG_ADD(LOG_FLOAT, posw_y, &posw_y)
// /**
//  * @brief Position z (m) - world frame
//  */
// LOG_ADD(LOG_FLOAT, posw_z, &posw_z)
// /**
//  * @brief Body Quaternion w - world frame
//  */
// LOG_ADD(LOG_FLOAT, qwb_w, &qwb_w)
// /**
//  * @brief Body Quaternion x - world frame
//  */
// LOG_ADD(LOG_FLOAT, qwb_x, &qwb_x)
// /**
//  * @brief Body Quaternion y - world frame
//  */
// LOG_ADD(LOG_FLOAT, qwb_y, &qwb_y)
// /**
//  * @brief Body Quaternion z - world frame
//  */
// LOG_ADD(LOG_FLOAT, qwb_z, &qwb_z)
// /**
//  * @brief Roll angle (deg)
//  */
// LOG_ADD(LOG_FLOAT, roll, &roll)
// /**
//  * @brief Pitch angle (deg)
//  */
// LOG_ADD(LOG_FLOAT, pitch, &pitch)
// /**
//  * @brief Yaw angle (deg)
//  */
// LOG_ADD(LOG_FLOAT, yaw, &yaw)
// /**
//  * @brief Linear velocity x (m/s) - body frame
//  */
// LOG_ADD(LOG_FLOAT, velb_x, &velb_x)
// /**
//  * @brief Linear velocity y (m/s) - body frame
//  */
// LOG_ADD(LOG_FLOAT, velb_y, &velb_y)
// /**
//  * @brief Linear velocity z (m/s) - body frame
//  */
// LOG_ADD(LOG_FLOAT, velb_z, &velb_z)
// /**
//  * @brief Angular velocity x (deg/s) - body frame
//  */
// LOG_ADD(LOG_FLOAT, omegab_x, &omegab_x)
// /**
//  * @brief Angular velocity y (deg/s) - body frame
//  */
// LOG_ADD(LOG_FLOAT, omegab_y, &omegab_y)
// /**
//  * @brief Angular velocity z (deg/s) - body frame
//  */
// LOG_ADD(LOG_FLOAT, omegab_z, &omegab_z)
// /**
//  * @brief Reference position x (m) - world frame
//  */
// LOG_ADD(LOG_FLOAT, posw_x_ref, &posw_x_ref)
// /**
//  * @brief Reference position y (m) - world frame
//  */
// LOG_ADD(LOG_FLOAT, posw_y_ref, &posw_y_ref)
// /**
//  * @brief Reference position z (m) - world frame
//  */
// LOG_ADD(LOG_FLOAT, posw_z_ref, &posw_z_ref)

// /**
//  * @brief Position x error
//  */
// LOG_ADD(LOG_FLOAT, posw_x_error, &posw_x_error)
// /**
//  * @brief Position y error
//  */
// LOG_ADD(LOG_FLOAT, posw_y_error, &posw_y_error)
// /**
//  * @brief Position z error
//  */
// LOG_ADD(LOG_FLOAT, posw_z_error, &posw_z_error)
// /**
//  * @brief Orientation x error
//  */
// LOG_ADD(LOG_FLOAT, orient_x_error, &orient_x_error)
// /**
//  * @brief Orientation y error
//  */
// LOG_ADD(LOG_FLOAT, orient_y_error, &orient_y_error)
// /**
//  * @brief Orientation z error
//  */
// LOG_ADD(LOG_FLOAT, orient_z_error, &orient_z_error)
// /**
//  * @brief Linear velocity x error
//  */
// LOG_ADD(LOG_FLOAT, velb_x_error, &velb_x_error)
// /**
//  * @brief Linear velocity y error
//  */
// LOG_ADD(LOG_FLOAT, velb_y_error, &velb_y_error)
// /**
//  * @brief Linear velocity z error
//  */
// LOG_ADD(LOG_FLOAT, velb_z_error, &velb_z_error)
// /**
//  * @brief Angular velocity x error
//  */
// LOG_ADD(LOG_FLOAT, omegab_x_error, &omegab_x_error)
// /**
//  * @brief Angular velocity y error
//  */
// LOG_ADD(LOG_FLOAT, omegab_y_error, &omegab_y_error)
// /**
//  * @brief Angular velocity z error
//  */
// LOG_ADD(LOG_FLOAT, omegab_z_error, &omegab_z_error)

// /**
//  * @brief Control speed (rad/sec) for motor 1
//  */
// LOG_ADD(LOG_FLOAT, ctrl_speed_m1, &ctrl_speed_m1)
// /**
//  * @brief Control speed (rad/sec) for motor 2
//  */
// LOG_ADD(LOG_FLOAT, ctrl_speed_m2, &ctrl_speed_m2)
// /**
//  * @brief Control speed (rad/sec) for motor 3
//  */
// LOG_ADD(LOG_FLOAT, ctrl_speed_m3, &ctrl_speed_m3)
// /**
//  * @brief Control speed (rad/sec) for motor 4
//  */
// LOG_ADD(LOG_FLOAT, ctrl_speed_m4, &ctrl_speed_m4)
// /**
//  * @brief Control thrust (N) for motor 1
//  */
// LOG_ADD(LOG_FLOAT, ctrl_thrust_m1, &ctrl_thrust_m1)
// /**
//  * @brief Control thrust (N) for motor 2
//  */
// LOG_ADD(LOG_FLOAT, ctrl_thrust_m2, &ctrl_thrust_m2)
// /**
//  * @brief Control thrust (N) for motor 3
//  */
// LOG_ADD(LOG_FLOAT, ctrl_thrust_m3, &ctrl_thrust_m3)
// /**
//  * @brief Control thrust (N) for motor 4
//  */
// LOG_ADD(LOG_FLOAT, ctrl_thrust_m4, &ctrl_thrust_m4)
// /**
//  * @brief Control normalized thrust for motor 1
//  */
// LOG_ADD(LOG_FLOAT, ctrl_norm_thrust_m1, &ctrl_norm_thrust_m1)
// /**
//  * @brief Control normalized thrust for motor 2
//  */
// LOG_ADD(LOG_FLOAT, ctrl_norm_thrust_m2, &ctrl_norm_thrust_m2)
// /**
//  * @brief Control normalized thrust for motor 3
//  */
// LOG_ADD(LOG_FLOAT, ctrl_norm_thrust_m3, &ctrl_norm_thrust_m3)
// /**
//  * @brief Control normalized thrust for motor 4
//  */
// LOG_ADD(LOG_FLOAT, ctrl_norm_thrust_m4, &ctrl_norm_thrust_m4)
LOG_GROUP_STOP(lqr_controller_logs)
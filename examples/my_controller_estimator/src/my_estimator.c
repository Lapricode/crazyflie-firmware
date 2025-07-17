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
 *
 * hello_world.c - App layer application of a simple hello world debug print every
 *   2 seconds.
 */

#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "stm32fxxx.h"
#include "FreeRTOS.h"
#include "queue.h"
#include "static_mem.h"
#include "task.h"

#include "estimator.h"
#include "estimator_kalman.h"

#include "log.h"
#include "param.h"
#include "math.h"
#include "math3d.h"
#include "debug.h"
#include "physicalConstants.h"
#include "power_distribution.h"
#include "platform_defaults.h"
#include "stabilizer_types.h"
#include "eventtrigger.h"

// Measurement models
#include "mm_distance.h"
#include "mm_absolute_height.h"
#include "mm_position.h"
#include "mm_pose.h"
#include "mm_tdoa.h"
#include "mm_flow.h"
#include "mm_tof.h"
#include "mm_yaw_error.h"
#include "mm_sweep_angles.h"
#include "mm_tdoa_robust.h"
#include "mm_distance_robust.h"

#define DEBUG_MODULE "MY_ESTIMATOR"
#define SIZE3 3
#define SIZE4 4
#define STATE_DIM 12
#define MEASURE_DIM 9
// #define MEASUREMENTS_QUEUE_SIZE (20)
// static xQueueHandle measurementsQueue;
// STATIC_MEM_QUEUE_ALLOC(measurementsQueue, MEASUREMENTS_QUEUE_SIZE, sizeof(measurement_t));

// // events
// EVENTTRIGGER(estPosition, uint8, source)
// EVENTTRIGGER(estSweepAngle, uint8, sensorId, uint8, baseStationId, uint8, sweepId, float, t, float, sweepAngle)
// EVENTTRIGGER(estGyroscope)
// EVENTTRIGGER(estAcceleration)

static bool isInit = false;
static const double tol = 1e-10f;

typedef struct mat33_s
{
  float m[3][3];
} mat33_t; // 3x3 matrix

typedef union vec3_u
{
  float v[SIZE3];
  struct
  {
    float x, y, z;
  };
} vec3_t; // 3x1 column vector

typedef struct vec4_u
{
  float v[SIZE4];
  struct
  {
    float v1, v2, v3, v4;
  };
} vec4_t; // 4x1 column vector

typedef struct vecS_u
{
  float v[STATE_DIM];
  struct
  {
    float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13;
  };
} vecS_t; // STATE_DIMx1 column vector

typedef struct vecM_u
{
  float v[STATE_DIM];
  struct
  {
    float v1, v2, v3, v4, v5, v6, v7, v8, v9;
  };
} vecM_t; // MEASURE_DIMx1 column vector

typedef union quat_u
{
  float q[SIZE4];
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

typedef struct matSS_s
{
  float m[STATE_DIM][STATE_DIM];
} matSS_t; // STATE_DIMxSTATE_DIM matrix

typedef struct matMM_s
{
  float m[MEASURE_DIM][MEASURE_DIM];
} matMM_t; // MEASURE_DIMxMEASURE_DIM matrix

typedef struct matSM_s
{
  float m[STATE_DIM][MEASURE_DIM];
} matSM_t; // STATE_DIMxMEASURE_DIM matrix

typedef struct matMS_s
{
  float m[MEASURE_DIM][STATE_DIM];
} matMS_t; // MEASURE_DIMxSTATE_DIM matrix

// define some logging variables
static unsigned int debug_print_counter = 0; // a counter for debug printing rate

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
static const mat33_t I_cf =
    {{{16.6e-6f, 0.83e-6f, 0.72e-6f},
      {0.83e-6f, 16.6e-6f, 1.8e-6f},
      {0.72e-6f, 1.8e-6f, 29.3e-6f}}}; // the crazyflie's moment of inertia
static const mat33_t I_inv_cf =
    {{{1.0f / 16.6e-6f, 0.0f, 0.0f},
      {0.0f, 1.0f / 16.6e-6f, 0.0f},
      {0.0f, 0.0f, 1.0f / 29.3e-6f}}}; // the crazyflie's inverted moment of inertia
static const float kf = 2.25e-08f;     // the coefficient parameter of the square model: thrust (N) vs. rotor_speed (rad/sec), for a single motor
static const float kt = 1.34e-10f;     // the coefficient parameter of the square model: torque (N*m) vs. rotor_speed (rad/sec), for a single motor, kt = 0.00596 * kf

// define variables and parameters for the Extended Kalman Filter (EKF)
static const uint32_t predict_rate = RATE_100_HZ;
static const float prediction_update_interval_sec = 1.0f / (float)predict_rate;
static const float dt = prediction_update_interval_sec;
static vec3_t vb_prev = {{0.0f}};
static matSS_t Q = {
    // process noise covariance
    {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
    }};
static matMM_t R = {
    // observation noise covariance
    {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
    }};
static matSS_t P = {{{0.0f}}};              // covariance estimate
static const float max_covariance = 100.0f; // maximum allowed covariance
static const float min_covariance = 1e-6f;  // minimum allowed covariance

// functions definitions
// static void kalman_predict(kalmanCoreData_t *this);
// static void kalman_update(kalmanCoreData_t *this);
static float clamp_value(float, float, float);
static mat33_t Rq_mat(quat_t);
static mat33_t vec3_hat(const vec3_t);
static vec3_t vec3_vec3_cross(const vec3_t, const vec3_t);
static vec3_t vec3_vec3_add(const vec3_t, const vec3_t);
static vec3_t vec3_vec3_sub(const vec3_t, const vec3_t);
static vec3_t vec3_scale(const vec3_t, const float);
static mat33_t mat33_mat33_add(const mat33_t, const mat33_t);
static mat33_t mat33_mat33_sub(const mat33_t, const mat33_t);
static vec3_t mat33_vec3_multiply(const mat33_t, const vec3_t);
static vecS_t matSM_vecM_multiply(const matSM_t, const vecM_t);
static mat33_t mat33_mat33_multiply(const mat33_t, const mat33_t);
static matSS_t matSS_matSS_multiply(const matSS_t, const matSS_t);
static matSM_t matSS_matSM_multiply(const matSS_t, const matSM_t);
static matMM_t matMS_matSM_multiply(const matMS_t, const matSM_t);
static matSS_t matSM_matMS_multiply(const matSM_t, const matMS_t);
static matSM_t matSM_matMM_multiply(const matSM_t, const matMM_t);
static mat33_t mat33_transpose(const mat33_t);
static matSS_t matSS_transpose(const matSS_t);
static matSM_t matMS_transpose(const matMS_t);
static matMM_t matMM_invert(const matMM_t);
static mat33_t SO3_plus_right(const mat33_t, const vec3_t);
static cf_state_t state_plus_right(const cf_state_t, const vecS_t);

static cf_state_t transition_model(const cf_state_t, const vec4_t);
static matSS_t transition_jacobian(const cf_state_t, const vec4_t);
static vecM_t observation_model(const cf_state_t);
static matMS_t observation_jacobian(const cf_state_t);
static cf_state_t predict_state(const cf_state_t, const vec4_t);
static matSS_t predict_covariance(const matSS_t, const matSS_t);
static matSM_t update_kalman_gain(const matSS_t, const matMS_t);
static cf_state_t update_state(const cf_state_t, const matSM_t, const vecM_t);
static matSS_t update_covariance(const matSS_t, const matSM_t, const matMS_t);

// clamp a float value betweeen two limit values
float clamp_value(float value, float min_lim, float max_lim)
{
  if (value < min_lim)
    return min_lim;
  if (value > max_lim)
    return max_lim;
  return value;
}

// convert a quaternion to the corresponding rotation matrix
// there is also the function "struct mat33 quat2rotmat(struct quat q)"" of "math3d.h"
mat33_t Rq_mat(quat_t q) // q is the quaternion
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

// compute the skew symmetric hat matrix that corresponds to the given vector
static mat33_t vec3_hat(const vec3_t v)
{
  mat33_t v_hat = {{{0.0f}}};
  v_hat.m[0][1] = -v.z;
  v_hat.m[0][2] = v.y;
  v_hat.m[1][0] = v.z;
  v_hat.m[1][2] = -v.x;
  v_hat.m[2][0] = -v.y;
  v_hat.m[2][1] = v.x;
  return v_hat;
}

// compute the cross producct of two 3d vectors
static vec3_t vec3_vec3_cross(const vec3_t a, const vec3_t b)
{
  vec3_t result;
  result.x = a.y * b.z - a.z * b.y;
  result.y = a.z * b.x - a.x * b.z;
  result.z = a.x * b.y - a.y * b.x;
  return result;
}

// compute the addition of two 3d vectors
static vec3_t vec3_vec3_add(const vec3_t a, const vec3_t b)
{
  vec3_t result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  result.z = a.z + b.z;
  return result;
}

// compute the subtraction of two 3d vectors
static vec3_t vec3_vec3_sub(const vec3_t a, const vec3_t b)
{
  vec3_t result;
  result.x = a.x - b.x;
  result.y = a.y - b.y;
  result.z = a.z - b.z;
  return result;
}

// compute the scaling of a 3d vector with a float value
static vec3_t vec3_scale(const vec3_t v, const float s)
{
  vec3_t result;
  result.x = v.x * s;
  result.y = v.y * s;
  result.z = v.z * s;
  return result;
}

// compute the addition of two 3x3 matrices
static mat33_t mat33_mat33_add(const mat33_t A, const mat33_t B)
{
  mat33_t C;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      C.m[i][j] = A.m[i][j] + B.m[i][j];
  return C;
}

// compute the subtraction of two 3x3 matrices
static mat33_t mat33_mat33_sub(const mat33_t A, const mat33_t B)
{
  mat33_t C;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      C.m[i][j] = A.m[i][j] - B.m[i][j];
  return C;
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

// matSM * vecM -> vecS
static vecS_t matSM_vecM_multiply(const matSM_t A, const vecM_t b)
{
  vecS_t result;
  for (int i = 0; i < STATE_DIM; i++)
  {
    result.v[i] = A.m[i][0] * b.v[0] + A.m[i][1] * b.v[1] + A.m[i][2] * b.v[2] +
                  A.m[i][3] * b.v[3] + A.m[i][4] * b.v[4] + A.m[i][5] * b.v[5] +
                  A.m[i][6] * b.v[6] + A.m[i][7] * b.v[7] + A.m[i][8] * b.v[8];
  }
  return result;
}

// compute the product A * B, where A is a 3x3 matrix and B is a 3x3 matrix
static mat33_t mat33_mat33_multiply(const mat33_t A, const mat33_t B)
{
  mat33_t C;
  for (int i = 0; i < SIZE3; i++)
  {
    for (int j = 0; j < SIZE3; j++)
    {
      C.m[i][j] = 0.0;
      for (int k = 0; k < SIZE3; k++)
        C.m[i][j] += A.m[i][k] * B.m[k][j];
    }
  }
  return C;
}

// matSS * matSS -> matSS
static matSS_t matSS_matSS_multiply(const matSS_t A, const matSS_t B)
{
  matSS_t C;
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      C.m[i][j] = 0.0;
      for (int k = 0; k < STATE_DIM; k++)
        C.m[i][j] += A.m[i][k] * B.m[k][j];
    }
  }
  return C;
}

// matSS * matSM -> matSM
static matSM_t matSS_matSM_multiply(const matSS_t A, const matSM_t B)
{
  matSM_t C;
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      C.m[i][j] = 0.0f;
      for (int k = 0; k < STATE_DIM; k++)
      {
        C.m[i][j] += A.m[i][k] * B.m[k][j];
      }
    }
  }
  return C;
}

// matMS * matSM -> matMM
static matMM_t matMS_matSM_multiply(const matMS_t A, const matSM_t B)
{
  matMM_t C;
  for (int i = 0; i < MEASURE_DIM; i++)
  {
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      C.m[i][j] = 0.0f;
      for (int k = 0; k < STATE_DIM; k++)
      {
        C.m[i][j] += A.m[i][k] * B.m[k][j];
      }
    }
  }
  return C;
}

// matSM * matMS -> matSS
static matSS_t matSM_matMS_multiply(const matSM_t A, const matMS_t B)
{
  matSS_t C;
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      C.m[i][j] = 0.0f;
      for (int k = 0; k < MEASURE_DIM; k++)
      {
        C.m[i][j] += A.m[i][k] * B.m[k][j];
      }
    }
  }
  return C;
}

// matSM * matMM -> matSM
static matSM_t matSM_matMM_multiply(const matSM_t A, const matMM_t B)
{
  matSM_t C;
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      C.m[i][j] = 0.0f;
      for (int k = 0; k < MEASURE_DIM; k++)
      {
        C.m[i][j] += A.m[i][k] * B.m[k][j];
      }
    }
  }
  return C;
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

// compute the STATE_DIMxSTATE_DIM transposed matrix A^T, where A is a STATE_DIMxSTATE_DIM matrix
static matSS_t matSS_transpose(const matSS_t A)
{
  matSS_t At;
  for (int i = 0; i < STATE_DIM; i++)
    for (int j = 0; j < STATE_DIM; j++)
      At.m[i][j] = A.m[j][i];
  return At;
}

// compute the STATE_DIMxMEASURE_DIM transposed matrix A^T, where A is a MEASURE_DIMxSTATE_DIM matrix
static matSM_t matMS_transpose(const matMS_t A)
{
  matSM_t At;
  for (int i = 0; i < STATE_DIM; i++)
    for (int j = 0; j < MEASURE_DIM; j++)
      At.m[i][j] = A.m[j][i];
  return At;
}

// compute the inverse A^(-1) of a STATE_DIMxSTATE_DIM matrix A using Gauss-Jordan elimination
static matMM_t matMM_invert(const matMM_t A)
{
  matMM_t A_copy = A;

  // initialize inverse matrix as the identity matrix
  matMM_t A_inv;
  for (int i = 0; i < MEASURE_DIM; i++)
  {
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      A_inv.m[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int i = 0; i < MEASURE_DIM; i++)
  {
    // check for zero diagonal element
    if (fabs(A_copy.m[i][i]) < tol)
    {
      // find a row to swap
      int swapped = 0;
      for (int k = i + 1; k < MEASURE_DIM; k++)
      {
        if (fabs(A_copy.m[k][i]) > tol)
        {
          // swap rows in both A and the inverse
          for (int j = 0; j < MEASURE_DIM; j++)
          {
            float tmp = A_copy.m[i][j];
            A_copy.m[i][j] = A_copy.m[k][j];
            A_copy.m[k][j] = tmp;

            tmp = A_inv.m[i][j];
            A_inv.m[i][j] = A_inv.m[k][j];
            A_inv.m[k][j] = tmp;
          }
          swapped = 1;
          break;
        }
      }
      if (!swapped) // singular matrix
      {
        return A_inv;
      }
    }

    // normalize the pivot row
    float temp;
    temp = A_copy.m[i][i];
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      A_copy.m[i][j] /= temp;
      A_inv.m[i][j] /= temp;
    }

    // eliminate other rows
    for (int k = 0; k < MEASURE_DIM; k++)
    {
      if (k == i)
        continue;
      temp = A_copy.m[k][i];
      for (int j = 0; j < MEASURE_DIM; j++)
      {
        A_copy.m[k][j] -= A_copy.m[i][j] * temp;
        A_inv.m[k][j] -= A_inv.m[i][j] * temp;
      }
    }
  }

  return A_inv;
}

// right plus operation for the orientation of the crazyflie
static mat33_t SO3_plus_right(const mat33_t R, const vec3_t tau)
{
  // compute Exp(tau) for SO(3)
  float theta = sqrtf(tau.x * tau.x + tau.y * tau.y + tau.z * tau.z);
  vec3_t u = {{0.0, 0.0, 1.0}};
  if ((double)theta >= tol)
  {
    u.x = tau.x / theta;
    u.y = tau.y / theta;
    u.z = tau.z / theta;
  }

  // compute R_exp = I + sin(theta)*u_hat + (1 - cos(theta))*u_hat^2
  mat33_t u_hat = vec3_hat(u);
  mat33_t u_hat_squared = mat33_mat33_multiply(u_hat, u_hat);
  mat33_t R_exp;
  for (int i = 0; i < SIZE3; ++i)
    for (int j = 0; j < SIZE3; ++j)
      R_exp.m[i][j] = (i == j ? 1.0f : 0.0f) + sinf(theta) * u_hat.m[i][j] + (1.0f - cosf(theta)) * u_hat_squared.m[i][j];

  mat33_t R_plus = mat33_mat33_multiply(R, R_exp);

  return R_plus;
}

// right plus operation for the state of the crazyflie (state + vector)
static cf_state_t state_plus_right(const cf_state_t s, const vecS_t ds)
{
  vec3_t rw = s.rw;
  mat33_t Rwb = s.Rwb;
  vec3_t vb = s.vb;
  vec3_t ob = s.ob;

  vec3_t drw = {{ds.v[0], ds.v[1], ds.v[2]}};
  vec3_t dRwb = {{ds.v[SIZE3], ds.v[SIZE3 + 1], ds.v[SIZE3 + 2]}};
  vec3_t dvb = {{ds.v[2 * SIZE3], ds.v[2 * SIZE3 + 1], ds.v[2 * SIZE3 + 2]}};
  vec3_t dob = {{ds.v[3 * SIZE3], ds.v[3 * SIZE3 + 1], ds.v[3 * SIZE3 + 2]}};

  vec3_t rw_sum = vec3_vec3_add(rw, drw);
  mat33_t Rwb_sum = SO3_plus_right(Rwb, dRwb);
  vec3_t vb_sum = vec3_vec3_add(vb, dvb);
  vec3_t ob_sum = vec3_vec3_add(ob, dob);

  cf_state_t state_plus = {
      .rw = rw_sum,
      .Rwb = Rwb_sum,
      .vb = vb_sum,
      .ob = ob_sum,
  };

  return state_plus;
}

// the nonlinear transition model, the crazyflie dynamics
static cf_state_t transition_model(const cf_state_t x, const vec4_t u)
{
  // vec3_t rw = x.rw;
  mat33_t Rwb = x.Rwb;
  vec3_t vb = x.vb;
  vec3_t ob = x.ob;

  // compute rw_dot
  vec3_t rw_dot = mat33_vec3_multiply(Rwb, vb);

  // compute vb_dot
  vec3_t gw = {{0.0f, 0.0f, -m_cf * g}};
  vec3_t gb = mat33_vec3_multiply(mat33_transpose(Rwb), gw);
  vec3_t thrust_b = {{0.0f, 0.0f, u.v[0] + u.v[1] + u.v[2] + u.v[3]}};
  vec3_t Fb = vec3_vec3_add(gb, thrust_b);
  vec3_t ob_cross_vb = vec3_vec3_cross(ob, vb);
  vec3_t vb_dot = vec3_vec3_sub(vec3_scale(Fb, 1.0f / m_cf), ob_cross_vb);

  // compute ob_dot
  float T13 = l_cf * (u.v[0] - u.v[2]);
  float T42 = l_cf * (u.v[3] - u.v[1]);
  vec3_t Tb = {{
      T13 * cosf(body_yaw0) - T42 * sinf(body_yaw0),
      T13 * sinf(body_yaw0) + T42 * cosf(body_yaw0),
      (kt / kf) * (-u.v[0] - u.v[2] + u.v[1] + u.v[3]),
  }};
  vec3_t Iob = mat33_vec3_multiply(I_cf, ob);
  vec3_t ob_cross_Iob = vec3_vec3_cross(ob, Iob);
  vec3_t torque_term = vec3_vec3_sub(Tb, ob_cross_Iob);
  vec3_t ob_dot = mat33_vec3_multiply(I_inv_cf, torque_term);

  // build dx_times_dt vector of size STATE_DIM
  vecS_t dx_times_dt;
  vec3_t ob_dt = vec3_scale(ob, dt);
  for (int i = 0; i < SIZE3; i++)
  {
    dx_times_dt.v[i] = rw_dot.v[i] * dt;
    dx_times_dt.v[SIZE3 + i] = ob_dt.v[i] * dt;
    dx_times_dt.v[2 * SIZE3 + i] = vb_dot.v[i] * dt;
    dx_times_dt.v[3 * SIZE3 + +i] = ob_dot.v[i] * dt;
  }

  // compute the next state estimate
  cf_state_t x_next = state_plus_right(x, dx_times_dt);
  return x_next;
}

static matSS_t transition_jacobian(const cf_state_t x, const vec4_t u)
{
  // vec3_t rw = x.rw;
  mat33_t Rwb = x.Rwb;
  vec3_t vb = x.vb;
  vec3_t ob = x.ob;

  matSS_t F = {{{0.0f}}};

  mat33_t vb_hat = vec3_hat(vb);
  mat33_t ob_hat = vec3_hat(ob);
  mat33_t drwdot_over_dRwb = mat33_mat33_multiply(Rwb, vb_hat);
  vec3_t gw = {{0.0f, 0.0f, -g}};
  mat33_t dvbdot_over_dRwb = mat33_mat33_multiply(mat33_transpose(Rwb), mat33_mat33_multiply(vec3_hat(gw), Rwb));
  mat33_t Iob_hat = vec3_hat(mat33_vec3_multiply(I_cf, ob));
  mat33_t ob_hat_times_I = mat33_mat33_multiply(ob_hat, I_cf);
  mat33_t dobdot_over_ob = mat33_mat33_multiply(I_inv_cf, mat33_mat33_sub(Iob_hat, ob_hat_times_I));

  // F[0:3][3:6]
  for (int i = 0; i < SIZE3; i++)
    for (int j = SIZE3; j < 2 * SIZE3; j++)
      F.m[i][j] = drwdot_over_dRwb.m[i][j];

  // F[0:3][6:9]
  for (int i = 0; i < SIZE3; i++)
    for (int j = 2 * SIZE3; j < 3 * SIZE3; j++)
      F.m[i][j] = Rwb.m[i][j];

  // # A[3:6, 3:6] = SO3.jacobian_rotation_action_1(Rwb, omegab)
  // # A[3:6, 9:12] = SO3.jacobian_rotation_action_2(Rwb, omegab)

  // // F[3:6][3:6]
  // for (int i = SIZE3; i < 2 * SIZE3; i++)
  //   for (int j = SIZE3; j < 2 * SIZE3; j++)
  //     F.m[i][j] = dvbdot_over_dRwb.m[i][j];

  // // F[3:6][9:12]
  // for (int i = SIZE3; i < 2 * SIZE3; i++)
  //   for (int j = 3 * SIZE3; j < 4 * SIZE3; j++)
  //     F.m[i][j] = dvbdot_over_dRwb.m[i][j];

  // F[6:9][3:6]
  for (int i = 2 * SIZE3; i < 3 * SIZE3; i++)
    for (int j = SIZE3; j < 2 * SIZE3; j++)
      F.m[i][j] = dvbdot_over_dRwb.m[i][j];

  // F[6:9][6:9]
  for (int i = 2 * SIZE3; i < 3 * SIZE3; i++)
    for (int j = 2 * SIZE3; j < 3 * SIZE3; j++)
      F.m[i][j] = -ob_hat.m[i][j];

  // F[6:9][9:12]
  for (int i = 2 * SIZE3; i < 3 * SIZE3; i++)
    for (int j = 3 * SIZE3; j < 4 * SIZE3; j++)
      F.m[i][j] = vb_hat.m[i][j];

  // F[9:12][9:12]
  for (int i = 3 * SIZE3; i < 4 * SIZE3; i++)
    for (int j = 3 * SIZE3; j < 4 * SIZE3; j++)
      F.m[i][j] = dobdot_over_ob.m[i][j];

  return F;
}

// the observation model, the crazyflie measurements
static vecM_t observation_model(const cf_state_t x)
{
  vec3_t rw = x.rw;
  // mat33_t Rwb = x.Rwb;
  vec3_t vb = x.vb;
  vec3_t ob = x.ob;

  vec3_t ab = {{0.0f}};
  for (int i = 0; i < SIZE3; i++)
  {
    ab.v[i] = (vb.v[i] - vb_prev.v[i]) / dt;
  }
  vb_prev.x = vb.x;
  vb_prev.y = vb.y;
  vb_prev.z = vb.z;

  vecM_t h = {{0.0f}};
  for (int i = 0; i < SIZE3; i++)
  {
    h.v[i] = rw.v[i];
    h.v[i + 3] = ob.v[i];
    h.v[i + 6] = ab.v[i];
  }
  return h;
}

static matMS_t observation_jacobian(const cf_state_t x)
{
  // vec3_t rw = x.rw;
  // mat33_t Rwb = x.Rwb;
  // vec3_t vb = x.vb;
  // vec3_t ob = x.ob;

  matMS_t H = {{{0.0f}}};
  for (int i = 0; i < SIZE3; i++)
  {
    H.m[i][i] = 1.0f;
    H.m[SIZE3 + i][STATE_DIM - SIZE3 + i] = 1.0f;
    H.m[MEASURE_DIM - SIZE3 + i][2 * SIZE3 + i] = 1.0f / dt;
  }
  return H;
}

static cf_state_t predict_state(const cf_state_t x_kpr_kpr, const vec4_t u_pr)
{
  cf_state_t x_k_kpr = transition_model(x_kpr_kpr, u_pr);
  return x_k_kpr;
}

static matSS_t predict_covariance(const matSS_t P_kpr_kpr, const matSS_t F_k)
{
  matSS_t F_k_tr = matSS_transpose(F_k);
  matSS_t temp = matSS_matSS_multiply(matSS_matSS_multiply(F_k, P_kpr_kpr), F_k_tr);
  matSS_t P_k_kpr = {{{0.0f}}};
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      P_k_kpr.m[i][j] = temp.m[i][j] + Q.m[i][j];
      P_k_kpr.m[i][j] = clamp_value(P_k_kpr.m[i][j], min_covariance, max_covariance);
    }
  }
  return P_k_kpr;
}

static matSM_t update_kalman_gain(matSS_t P_k_kpr, matMS_t H_k)
{
  matSM_t H_k_tr = matMS_transpose(H_k);
  matSM_t temp1 = matSS_matSM_multiply(P_k_kpr, H_k_tr);
  matMM_t temp2 = matMS_matSM_multiply(H_k, temp1);
  matMM_t temp3 = {{{0.0f}}};
  for (int i = 0; i < MEASURE_DIM; i++)
  {
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      temp3.m[i][j] = temp2.m[i][j] + R.m[i][j];
    }
  }
  matSM_t K_k = matSM_matMM_multiply(temp1, matMM_invert(temp3));
  return K_k;
}

static cf_state_t update_state(const cf_state_t x_k_kpr, const matSM_t K_k, const vecM_t z_k)
{
  vecM_t temp1;
  vecM_t h_k_kpr = observation_model(x_k_kpr);
  for (int i = 0; i < MEASURE_DIM; i++)
  {
    temp1.v[i] = z_k.v[i] - h_k_kpr.v[i];
  }
  vecS_t temp2 = matSM_vecM_multiply(K_k, temp1);
  cf_state_t x_k_k = state_plus_right(x_k_kpr, temp2);
  return x_k_k;
}

static matSS_t update_covariance(const matSS_t P_k_kpr, const matSM_t K_k, const matMS_t H_k)
{
  matSS_t temp1 = matSM_matMS_multiply(K_k, H_k);
  matSS_t temp2;
  for (int i = 0; i < MEASURE_DIM; i++)
  {
    for (int j = 0; j < MEASURE_DIM; j++)
    {
      if (i == j)
        temp2.m[i][j] = 1.0f - temp1.m[i][j];
      else
        temp2.m[i][j] = -temp1.m[i][j];
    }
  }
  matSS_t P_k_k = matSS_matSS_multiply(temp2, P_k_kpr);
  for (int i = 0; i < STATE_DIM; i++)
    for (int j = 0; j < STATE_DIM; j++)
      P_k_k.m[i][j] = clamp_value(P_k_k.m[i][j], min_covariance, max_covariance);
  return P_k_k;
}

void estimatorOutOfTreeInit(void)
{
  // estimatorKalmanInit();
  // initialize the covariances
  isInit = true;
  return;
}

bool estimatorOutOfTreeTest(void)
{
  // return estimatorKalmanTest();
  return isInit;
}

void estimatorOutOfTree(state_t *state, const stabilizerStep_t tick)
{
  // estimatorKalman(state, tick);

  // get the current state of the crazyflie
  vec3_t vw_cur = {{state->velocity.x, state->velocity.y, state->velocity.z}};
  quat_t qwb_cur = {{state->attitudeQuaternion.w, state->attitudeQuaternion.x, state->attitudeQuaternion.y, state->attitudeQuaternion.z}};
  mat33_t Rwb_cur = Rq_mat(qwb_cur);
  vec3_t vb_cur = mat33_vec3_multiply(mat33_transpose(Rwb_cur), vw_cur);
  cf_state_t x_kpr_kpr = {
      .rw = {{state->position.x, state->position.y, state->position.z}},
      .Rwb = Rwb_cur,
      .vb = vb_cur,
      .ob = {{0.0f}}, // need to fix this
                      // .ob = {{radians(sensors->gyro.x), radians(sensors->gyro.y), radians(sensors->gyro.z)}},
  };

  // get the current control input
  vec4_t u_kpr = {{0.0f}};

  // get the current measurements
  vecM_t z_k = {{0.0f}};

  // predict estimate
  cf_state_t x_k_kpr = predict_state(x_kpr_kpr, u_kpr);
  matSS_t F_k = transition_jacobian(x_kpr_kpr, u_kpr);
  P = predict_covariance(P, F_k);

  // update estimate
  matMS_t H_k = observation_jacobian(x_k_kpr);
  matSM_t K_k = update_kalman_gain(P, H_k);
  cf_state_t x_k_k = update_state(x_k_k, K_k, z_k);
  P = update_covariance(P, K_k, H_k);

  // bool quadIsFlying = supervisorIsFlying();

  // print some data for debugging
  if (RATE_DO_EXECUTE(1, debug_print_counter))
  {
    debug_print_counter = 0;
    DEBUG_PRINT("My EKF Estimator is running!\n");
  }
  debug_print_counter += 1;

  return;
}
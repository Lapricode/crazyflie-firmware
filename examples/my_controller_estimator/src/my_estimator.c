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

#include "app.h"
#include "FreeRTOS.h"
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

static bool isInit = false;
static const double tol = 1e-10f;

typedef union vec_3_u
{
  float v[SIZE3];
  struct
  {
    float x, y, z;
  };
} vec_3_t; // 3x1 column vector

typedef struct vec_4_u
{
  float v[SIZE4];
  struct
  {
    float v1, v2, v3, v4;
  };
} vec_4_t; // 4x1 column vector

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
  vec_3_t rw; // position w.r.t. world frame
  quat_t qwb; // body orientation in quaternion form w.r.t. world frame
  vec_3_t vb; // linear velocity w.r.t body frame
  vec_3_t ob; // angular velocity w.r.t body frame
} cf_state_t; // crazyflie's state structure

typedef struct matSS_s
{
  float m[STATE_DIM][STATE_DIM];
} matSS_t; // STATE_DIMxSTATE_DIM matrix

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

// // define some crazyflie model parameters
// // Quadrotor system:
// // state:      x = [rw, qwb, vb, omegab] \in R^(13x1)
// //             where   rw \in R^(3x1) is the position in the world frame
// //                     qwb \in R^(4x1) is the orientation quaternion of the body frame w.r.t. the world frame
// //                     vb \in R^(3x1) is the linear velocity in the body frame
// //                     omegab \in R^(3x1) is the angular velocity in the body frame
// // control:    u = [u1, u2, u3, u4] \in R^(4x1)
// //             where u_i is the angular speed of the i-th motor in rad/s
// // rotors:     Fi = Kf * ui^2
// //             Ti = Kt * ui^2
// //             where   Fi is the thrust force produced by the ith rotor
// //                     Ti is the torque produced by the ith rotor
// //             the motors are numbered in a clockwise manner, with motor 4 being in xy direction
// //             the rotors 1, 3 rotate counter-clockwise, and the rotors 2, 4 rotate clockwise
// //             the motor arms form right angles (90 degrees) with each other
// // note:   the state x comes with the quaternion qwb of the rotation matrix Rwb,
// //         but we use the rotation matrix Rwb directly for the dynamics and the jacobians calculations
// //         Rwb \in R^(3x3) is the rotation matrix of the body frame w.r.t. the world frame
// //         Rwb \in R^(3x3) has dimension 3
// static const float g = 9.81f;                              // gravity's acceleration (in m/sec^2)
// static const float m_cf = 0.033f;                          // mass (in kg)
// static const float l = 0.046f;                             // arm length (in m)
// static const float body_yaw0 = -3.0f / 4.0f * (float)M_PI; // assuming body_yaw0 is for the motor 1 at positive y direction, motor 2 at positive x direction and clockwise motor numbers
// // static const mat_3_3_t CRAZYFLIE_INERTIA =
// //     {{{16.6e-6f, 0.83e-6f, 0.72e-6f},
// //       {0.83e-6f, 16.6e-6f, 1.8e-6f},
// //       {0.72e-6f, 1.8e-6f, 29.3e-6f}}};
// static const float kf = 2.25e-08f; // the coefficient parameter of the square model: thrust (N) vs. rotor_speed (rad/sec), for a single motor
// static const float kt = 1.34e-10f; // the coefficient parameter of the square model: torque (N*m) vs. rotor_speed (rad/sec), for a single motor, kt = 0.00596 * kf

// define parameters for the Extended Kalman Filter (EKF)
static const uint32_t predict_rate = RATE_100_HZ;
static const float prediction_update_interval_ms = 1000.0f / (float)predict_rate;
static float Q[STATE_DIM][STATE_DIM] = {
    // process noise covariance
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
};
static float R[MEASURE_DIM][MEASURE_DIM] = {
    // observation noise covariance
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
};
// static float P[STATE_DIM][STATE_DIM];       // covariance estimate
// static float S[STATE_DIM][STATE_DIM];       // innovation covariance
static const float max_covariance = 100.0f; // maximum allowed covariance
static const float min_covariance = 1e-6f;  // minimum allowed covariance

// static void kalman_predict(kalmanCoreData_t *this);
// static void kalman_update(kalmanCoreData_t *this);
static float clamp_value(float, float, float);
static vec_3_t mat33_vec3_multiply(const mat_3_3_t, const vec_3_t);
static void matSS_matSS_multiply(const float[STATE_DIM][STATE_DIM], const float[STATE_DIM][STATE_DIM], float[STATE_DIM][STATE_DIM]);
static void matSS_matSM_multiply();
static void matMS_matSM_multiply();
static void matSM_matMS_multiply();
static void matSS_transpose(const float[STATE_DIM][STATE_DIM], float[STATE_DIM][STATE_DIM]);
static void matMS_transpose(const float[MEASURE_DIM][STATE_DIM], float[STATE_DIM][MEASURE_DIM]);
static void matS_invert(const float[STATE_DIM][STATE_DIM], float[STATE_DIM][STATE_DIM]);
static void SO3_plus_right(const float[SIZE3][SIZE3], const float[SIZE3], float[SIZE3][SIZE3]);
static cf_state_t predict_state(const cf_state_t, const vec_4_t);
static matSS_t predict_covariance(const);
static matSM_t update_kalman_gain();
static cf_state_t update_state(const cf_state_t, const float[STATE_DIM][MEASURE_DIM], const float[MEASURE_DIM]);
static matSS_t update_covariance(void);

// clamp a float value betweeen two limit values
float clamp_value(float value, float min_lim, float max_lim)
{
  if (value < min_lim)
    return min_lim;
  if (value > max_lim)
    return max_lim;
  return value;
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

// compute the STATE_DIMxSTATE_DIM matrix multiplication C = A * B, where A, B are square matrices of size STATE_DIMxSTATE_DIM
static void matSS_matSS_multiply(const float A[STATE_DIM][STATE_DIM], const float B[STATE_DIM][STATE_DIM], float C[STATE_DIM][STATE_DIM])
{
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < STATE_DIM; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
  }
  return;
}

// compute the STATE_DIMxSTATE_DIM matrix multiplication C = A * B, where A, B are square matrices of size STATE_DIMxSTATE_DIM
static void matSS_matS_multiply(const float A[STATE_DIM][STATE_DIM], const float B[STATE_DIM][STATE_DIM], float C[STATE_DIM][STATE_DIM])
{
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < STATE_DIM; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
  }
  return;
}

// compute the STATE_DIMxSTATE_DIM matrix multiplication C = A * B, where A, B are square matrices of size STATE_DIMxSTATE_DIM
static void matSM_matMS_multiply(const float A[STATE_DIM][STATE_DIM], const float B[STATE_DIM][STATE_DIM], float C[STATE_DIM][STATE_DIM])
{
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < STATE_DIM; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
  }
  return;
}

// compute the STATE_DIMxSTATE_DIM matrix multiplication C = A * B, where A, B are square matrices of size STATE_DIMxSTATE_DIM
static void matMS_matSM_multiply(const float A[STATE_DIM][STATE_DIM], const float B[STATE_DIM][STATE_DIM], float C[STATE_DIM][STATE_DIM])
{
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < STATE_DIM; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
  }
  return;
}

// compute the STATE_DIMxSTATE_DIM transposed matrix A^T, where A is a STATE_DIMxSTATE_DIM matrix
static void matSS_transpose(const float[STATE_DIM][STATE_DIM], float[STATE_DIM][STATE_DIM])
{
  float At[STATE_DIM][STATE_DIM];
  for (int i = 0; i < STATE_DIM, i++)
    for (int j = 0; j < STATE_DIM, j++)
      At[i][j] = A[j][i];
  return;
}

// compute the STATE_DIMxMEASURE_DIM transposed matrix A^T, where A is a MEASURE_DIMxSTATE_DIM matrix
static void matMS_transpose(const float[MEASURE_DIM][STATE_DIM], float[STATE_DIM][MEASURE_DIM])
{
  float At[STATE_DIM][MEASURE_DIM];
  for (int i = 0; i < STATE_DIM, i++)
    for (int j = 0; j < MEASURE_DIM, j++)
      At[i][j] = A[j][i];
  return;
}

// compute the inverse A^(-1) of a STATE_DIMxSTATE_DIM matrix A using Gauss-Jordan elimination
static void matS_invert(const float S[STATE_DIM][STATE_DIM], float S_inv[STATE_DIM][STATE_DIM])
{
  float S_copy[STATE_DIM][STATE_DIM];
  memcpy(S_copy, S, sizeof(float) * STATE_DIM * STATE_DIM);

  // initialize inverse matrix as the identity matrix
  for (int i = 0; i < STATE_DIM; i++)
  {
    for (int j = 0; j < STATE_DIM; j++)
    {
      S_inv[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int i = 0; i < STATE_DIM; i++)
  {
    // check for zero diagonal element
    if (fabs(S_copy[i][i]) < tol)
    {
      // find a row to swap
      int swapped = 0;
      for (int k = i + 1; k < STATE_DIM; k++)
      {
        if (fabs(S_copy[k][i]) > tol)
        {
          // swap rows in both S and the inverse
          for (int j = 0; j < STATE_DIM; j++)
          {
            float tmp = S_copy[i][j];
            S_copy[i][j] = S_copy[k][j];
            S_copy[k][j] = tmp;

            tmp = S_inv[i][j];
            S_inv[i][j] = S_inv[k][j];
            S_inv[k][j] = tmp;
          }
          swapped = 1;
          break;
        }
      }
      if (!swapped) // singular matrix
      {
        return;
      }
    }

    // normalize the pivot row
    float temp;
    temp = S_copy[i][i];
    for (int j = 0; j < STATE_DIM; j++)
    {
      S_copy[i][j] /= temp;
      S_inv[i][j] /= temp;
    }

    // eliminate other rows
    for (int k = 0; k < STATE_DIM; k++)
    {
      if (k == i)
        continue;
      temp = S_copy[k][i];
      for (int j = 0; j < STATE_DIM; j++)
      {
        S_copy[k][j] -= S_copy[i][j] * temp;
        S_inv[k][j] -= S_inv[i][j] * temp;
      }
    }
  }

  return;
}

// right plus operation for the orientation of the quadrotor
static void SO3_plus_right(const float R[SIZE3][SIZE3], const float tau[SIZE3], float R_plus[SIZE3][SIZE3])
{
  float R_exp[SIZE3][SIZE3];

  // compute Exp(tau) for SO(3)
  float theta = sqrtf(tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
  float u[3] = {0.0, 0.0, 1.0};
  if ((double)theta >= tol)
  {
    u[0] = tau[0] / theta;
    u[1] = tau[1] / theta;
    u[2] = tau[2] / theta;
  }

  // compute R_exp = I + sin(theta)*u_hat + (1 - cos(theta))*u_hat^2
  float u_hat[3][3];
  u_hat[0][0] = 0.0;
  u_hat[0][1] = -u[2];
  u_hat[0][2] = u[1];
  u_hat[1][0] = u[2];
  u_hat[1][1] = 0.0;
  u_hat[1][2] = -u[0];
  u_hat[2][0] = -u[1];
  u_hat[2][1] = u[0];
  u_hat[2][2] = 0.0;
  float u_hat_squared[3][3];
  mat_3_3_mult(u_hat, u_hat, u_hat_squared);
  for (int i = 0; i < SIZE3; ++i)
    for (int j = 0; j < SIZE3; ++j)
      R_exp[i][j] = (i == j ? 1.0f : 0.0f) + sinf(theta) * u_hat[i][j] + (1.0f - cosf(theta)) * u_hat_squared[i][j];

  mat_3_3_mult(R, R_exp, R_plus);

  return;
}

static cf_state_t predict_state(const cf_state_t x_kpr_kpr, const vec_4_t u)
{
  cf_state_t x_k_kpr;
  return x_k_kpr;
}

static matSS_t predict_covariance()
{
  return;
}

static matSM_t update_kalman_gain()
{
  return;
}

static cf_state_t update_state(const cf_state_t x_k_kpr, const float K[STATE_DIM][MEASURE_DIM], const float zk[MEASURE_DIM])
{
  cf_state_t x_k_k;
  return x_k_k;
}

static matSS_t update_covariance(void)
{
  return;
}

void estimatorOutOfTreeInit(void)
{
  estimatorKalmanInit();
  // initialize the covariances
  isInit = true;
  return;
}

bool estimatorOutOfTreeTest(void)
{
  return estimatorKalmanTest();
  return isInit;
}

void estimatorOutOfTree(state_t *state, const stabilizerStep_t tick)
{
  estimatorKalman(state, tick);

  // print some data for debugging
  if (RATE_DO_EXECUTE(1, debug_print_counter))
  {
    debug_print_counter = 0;
    DEBUG_PRINT("My EKF Estimator is running!\n");
  }
  debug_print_counter += 1;

  return;
}
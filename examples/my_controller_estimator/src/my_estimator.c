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
#define SIZE9 9
#define SIZE12 12

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
static float Q[SIZE12][SIZE12] = {
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
static float R[SIZE9][SIZE9] = {
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
static float P[SIZE12][SIZE12];             // covariance estimate
static float S[SIZE12][SIZE12];             // innovation covariance
static const float max_covariance = 100.0f; // maximum allowed covariance
static const float min_covariance = 1e-6f;  // minimum allowed covariance

// static void kalman_predict(kalmanCoreData_t *this);
// static void kalman_update(kalmanCoreData_t *this);

// compute the 3x3 matrix multiplication C = A * B, where A, B are square matrices of size 3x3
void mat_3_3_mult(const float A[SIZE3][SIZE3], const float B[SIZE3][SIZE3], float C[SIZE3][SIZE3])
{
  for (int i = 0; i < SIZE3; i++)
  {
    for (int j = 0; j < SIZE3; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
  }
}

// compute the 12x12 matrix multiplication C = A * B, where A, B are square matrices of size 12x12
void mat_12_12_mult(const float A[SIZE12][SIZE12], const float B[SIZE12][SIZE12], float C[SIZE12][SIZE12])
{
  for (int i = 0; i < SIZE12; i++)
  {
    for (int j = 0; j < SIZE12; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < SIZE12; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
  }
}

// compute the inverse S^(-1) of a 12x12 matrix S using Gauss-Jordan elimination
void invert_matrix(float S[SIZE12][SIZE12], float S_inv[SIZE12][SIZE12])
{
  float S_copy[SIZE12][SIZE12];
  memcpy(S_copy, S, sizeof(float) * SIZE12 * SIZE12);

  // initialize inverse matrix as the identity matrix
  for (int i = 0; i < SIZE12; i++)
  {
    for (int j = 0; j < SIZE12; j++)
    {
      S_inv[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int i = 0; i < SIZE12; i++)
  {
    // check for zero diagonal element
    if (fabs(S_copy[i][i]) < EPSILON)
    {
      // find a row to swap
      int swapped = 0;
      for (int k = i + 1; k < SIZE12; k++)
      {
        if (fabs(S_copy[k][i]) > EPSILON)
        {
          // swap rows in both S and the inverse
          for (int j = 0; j < SIZE12; j++)
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
        return 0;
      }
    }

    // normalize the pivot row
    float temp;
    temp = S_copy[i][i];
    for (int j = 0; j < SIZE12; j++)
    {
      S_copy[i][j] /= temp;
      S_inv[i][j] /= temp;
    }

    // eliminate other rows
    for (int k = 0; k < SIZE12; k++)
    {
      if (k == i)
        continue;
      temp = S_copy[k][i];
      for (int j = 0; j < SIZE12; j++)
      {
        S_copy[k][j] -= S_copy[i][j] * temp;
        S_inv[k][j] -= S_inv[i][j] * temp;
      }
    }
  }

  return;
}

// right plus operation for the orientation of the quadrotor
void right_plus(const float R[SIZE3][3], const float tau[SIZE3], float R_plus[SIZE3][SIZE3])
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
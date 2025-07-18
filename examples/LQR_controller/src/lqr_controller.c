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

static const float kf = 2.25e-08f; // the coefficient parameter of the square model: thrust (N) vs. rotor_speed (rad/sec), for a single motor
static const float kt = 1.34e-10f; // the coefficient parameter of the square model: torque (N*m) vs. rotor_speed (rad/sec), for a single motor, kt = 0.00596 * kf

// define the LQR controller variables
static float Kinf[4][12] = {{0.0f}}; // initialize the LQR controller's Kinf matrix

// our custom Kinf constant LQR controller's matrix
static const float custom_Kinf[4][12] = {
    {-4.34612169e+02f, 2.25644573e+02f, 1.07444819e+03f, -9.41655365e+02f, -1.20221510e+03f, -8.49498402e+02f, -3.56292177e+02f, 2.43187738e+02f, 7.49740794e+02f, -1.74737382e+02f, -1.98468534e+02f, -6.25310833e+02f},
    {4.04150916e+02f, 1.71361278e+02f, 1.07444819e+03f, -7.09858035e+02f, 1.11832284e+03f, 8.47710105e+02f, 3.31365009e+02f, 1.84191075e+02f, 7.49740794e+02f, -1.30593138e+02f, 1.84689421e+02f, 6.23855890e+02f},
    {2.77671758e+02f, -1.89076156e+02f, 1.07444819e+03f, 7.83778559e+02f, 7.52894884e+02f, -8.41753850e+02f, 2.25767115e+02f, -2.03282872e+02f, 7.49740794e+02f, 1.44306492e+02f, 1.21498330e+02f, -6.19160783e+02f},
    {-2.47210505e+02f, -2.07929696e+02f, 1.07444819e+03f, 8.67734841e+02f, -6.69002624e+02f, 8.43542147e+02f, -2.00839947e+02f, -2.24095941e+02f, 7.49740794e+02f, 1.61024028e+02f, -1.07719217e+02f, 6.20615726e+02f}};

static float state_error[12] = {0.0f}; // the state error
// static vec_4_t hover_speeds = {{1900.0f, 1900.0f, 1900.0f, 1900.0f}}; // the angular speeds (in rad/sec) of the 4 rotors, for the crazyflie to hover
static vec_4_t hover_speeds;                                              // the angular speeds (in rad/sec) of the 4 rotors, for the crazyflie to hover
static float hover_adjust = 0.0f;                                         // adjust hover speeds for hover calibration
static vec_4_t control_speeds = {{0.0f, 0.0f, 0.0f, 0.0f}};               // the controlled angular speeds (in rad/sec) of the 4 rotors
static float max_control_speed = 25000.0f * (2.0f * (float)M_PI / 60.0f); // the maximum angular speed (in rad/sec) of a rotor, approximately 2618.0f rad/sec
static vec_4_t control_thrusts = {{0.0f, 0.0f, 0.0f, 0.0f}};              // the controlled thrusts (in N) of the 4 rotors
static unsigned int update_rate = RATE_HL_COMMANDER;                      // RATE_HL_COMMANDER;                      // the update rate of the control loop (100 Hz default rate)
static bool do_norm_forces_control = true;                                // if true, then controlModeForce, else controlModeForceTorque

// for the forces-torques control mode (controlModeForceTorque), do_norm_forces_control = false
static float control_thrust_total = 0.0f;                   // the total thrust (in N) generated by the 4 rotors
static vec_3_t control_body_torques = {{0.0f, 0.0f, 0.0f}}; // the body torques for each axis (x, y, z)

// for the normalized forces control mode (controlModeForce), do_norm_forces_control = true
static const float max_thrust = 0.156f;                           // the maximum thrust (in N) generated by only 1 motor
static vec_4_t control_norm_thrusts = {{0.0f, 0.0f, 0.0f, 0.0f}}; // the controlled normalized thrusts, in [0, 1], of the 4 rotors

// functions definitions
// float capAngle(float);
float clamp_to_unit_interval(float);
quat_t qrpy_quat(vec_3_t);
mat_3_3_t Rq_mat(quat_t);
static mat_3_3_t mat33_transpose(mat_3_3_t);
static vec_3_t mat33_vec3_multiply(mat_3_3_t, vec_3_t);
static vec_3_t SO3_minus_right(mat_3_3_t, mat_3_3_t);
static void compute_state_error(const cf_state_t, const cf_state_t, float *);
void hover_control_init(void);
void Kinf_LQR_init();

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
static mat_3_3_t mat33_transpose(mat_3_3_t A)
{
  mat_3_3_t At;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      At.m[i][j] = A.m[j][i];
  return At;
}

// compute the product A * b, where A is a 3x3 matrix and b is a 3x1 column vector
static vec_3_t mat33_vec3_multiply(mat_3_3_t A, vec_3_t b)
{
  vec_3_t result;
  for (int i = 0; i < 3; i++)
  {
    result.v[i] = A.m[i][0] * b.v[0] + A.m[i][1] * b.v[1] + A.m[i][2] * b.v[2];
  }
  return result;
}

// compute the 4x1 column vector product A * b, where A is a 4x12 matrix and b is a 12x1 column vector
static void mat412_vec12_multiply(const float A[4][12], const float b[12], vec_4_t *result)
{
  for (int i = 0; i < 4; i++)
  {
    float sum = 0.0f;
    for (int j = 0; j < 12; j++)
    {
      sum += A[i][j] * b[j];
    }
    result->v[i] = sum;
  }
  return;
}

// compute the product A * B, where A, B are 3x3 matrices
static mat_3_3_t mat33_mat33_multiply(mat_3_3_t A, mat_3_3_t B)
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
static vec_3_t SO3_minus_right(mat_3_3_t R1, mat_3_3_t R2)
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
static void compute_state_error(const cf_state_t state_cur, const cf_state_t state_ref, float *error)
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

  // concatenate into error vector
  for (int i = 0; i < 3; i++)
  {
    error[i] = rw_error[i];
    error[i + 3] = Rwb_error.v[i];
    error[i + 6] = vb_error[i];
    error[i + 9] = ob_error[i];
  }

  return;
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

void Kinf_LQR_init()
{
  memcpy(Kinf, custom_Kinf, sizeof(Kinf));
  return;
}

void appMain()
{
  DEBUG_PRINT("Waiting for activation of the controller and the estimator ...\n");

  while (1)
  {
    vTaskDelay(M2T(2000));
  }
}

void controllerOutOfTreeInit(void)
{
  hover_control_init();
  Kinf_LQR_init();
  isInit = true;
  return;
}

bool controllerOutOfTreeTest(void)
{
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
    compute_state_error(state_cur, state_ref, state_error);
    vec_4_t control_feedback;
    mat412_vec12_multiply(Kinf, state_error, &control_feedback);
    for (int i = 0; i < 4; i++)
    {
      control_speeds.v[i] = hover_speeds.v[i] - control_feedback.v[i];
      if (control_speeds.v[i] < 0.0f)
        control_speeds.v[i] = 0.0f;
      if (control_speeds.v[i] > max_control_speed)
        control_speeds.v[i] = max_control_speed;
      control_thrusts.v[i] = kf * powf(control_speeds.v[i], 2.0f);
    }
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
  }
  debug_print_counter += 1;
}

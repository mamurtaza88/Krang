/*
 * Copyright (c) 2014-2016, Humanoid Lab, Georgia Tech Research Corporation
 * Copyright (c) 2014-2017, Graphics Lab, Georgia Tech Research Corporation
 * Copyright (c) 2016-2017, Personal Robotics Lab, Carnegie Mellon University
 * All rights reserved.
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include "Controller.hpp"
#define MAXBUFSIZE  ((int) 1e6)

Eigen::MatrixXd readMatrix(const char *filename)
{
  int cols = 0, rows = 0;
  double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
  std::ifstream infile;
  infile.open(filename);
  while (! infile.eof())
  {
    std::string line;
    getline(infile, line);

    int temp_cols = 0;
    std::stringstream stream(line);
    while(! stream.eof())
      stream >> buff[cols*rows+temp_cols++];

      if (temp_cols == 0)
        continue;

      if (cols == 0)
        cols = temp_cols;

        rows++;
  }

  infile.close();

  rows--;

    // Populate matrix with numbers.
  Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
     for (int j = 0; j < cols; j++)
        result(i,j) = buff[ cols*i+j ];

  return result;
};


//==============================================================================
Controller::Controller(dart::dynamics::SkeletonPtr _robot,
                       dart::dynamics::BodyNode* _LeftendEffector,
                       dart::dynamics::BodyNode* _RightendEffector)
    : mRobot(_robot),
      mLeftEndEffector(_LeftendEffector),
      mRightEndEffector(_RightendEffector) {
  assert(_robot != nullptr);
  assert(_LeftendEffector != nullptr);
  assert(_RightendEffector != nullptr);

  int dof = mRobot->getNumDofs();
  // std::cout << "[controller] DoF: " << dof << std::endl;

  mForces.setZero(19);

  mSteps = 0;
  // *************** Read Initial Pose for Pose Regulation and Generate
  // reference zCOM
  Eigen::Matrix<double, 25, 1> qInit;// Eigen::Matrix<double, 25, 1> qInit;
  Eigen::Matrix3d Rot0;
  qInit = mRobot->getPositions();
  mBaseTf = mRobot->getBodyNode(0)->getTransform().matrix();
  double psiInit = atan2(mBaseTf(0, 0), -mBaseTf(1, 0));
  Rot0 << cos(psiInit), sin(psiInit), 0, -sin(psiInit), cos(psiInit), 0, 0, 0,
      1;
  double qBody1Init =
      atan2(mBaseTf(0, 1) * cos(psiInit) + mBaseTf(1, 1) * sin(psiInit),
            mBaseTf(2, 1));
  mqBodyInit(0) = qBody1Init;
  mqBodyInit.tail(17) = qInit.tail(17);
  dart::dynamics::BodyNode* LWheel = mRobot->getBodyNode("LWheel");
  dart::dynamics::BodyNode* RWheel = mRobot->getBodyNode("RWheel");
  Eigen::Vector3d bodyCOM = Rot0 * (mRobot->getCOM() - qInit.segment(3, 3));
  bodyCOM(1) = 0;
  mInitCOMDistance = bodyCOM.norm();

  // ************** Remove position limits
  for (int i = 6; i < dof - 1; ++i)
    _robot->getJoint(i)->setPositionLimitEnforced(false);
  // std::cout << "Position Limit Enforced set to false" << std::endl;

  // ************** Set joint damping
  for (int i = 6; i < dof - 1; ++i)
    _robot->getJoint(i)->setDampingCoefficient(0, 0.5);
  // std::cout << "Damping coefficients set" << std::endl;

  mdqFilt = new filter(25, 100);

  // ************** Wheel Radius and Distance between wheels
  mR = 0.25, mL = 0.68;

  // *********************************** Tunable Parameters
  config4cpp::Configuration* cfg = config4cpp::Configuration::create();
  const char* scope = "";
  const char* configFile = "/home/krang/SethResearch/13EE2-3DUnification-UnlockedJoints-Murtaza/examples/3dofddp/controlParams.cfg";
  const char* str;
  std::istringstream stream;
  double newDouble;
  Eigen::Matrix<double, 18, 1> tauLim;

  mKpEE.setZero();
  mKvEE.setZero();
  mWEER.setZero();
  mWEEL.setZero();
  mWBal.setZero();
  mWMatPose.setZero();
  mWMatSpeedReg.setZero();
  mWMatReg.setZero();
  try {
    cfg->parse(configFile);

    // Waist Locked?
    mWaistLocked = cfg->lookupBoolean(scope, "waistLocked");

    // -- COM Angle Based Control or not
    mCOMAngleControl = cfg->lookupBoolean(scope, "COMAngleControl");
    mMaintainInitCOMDistance =
        cfg->lookupBoolean(scope, "maintainInitCOMDistance");

    // -- Torque Limits
    str = cfg->lookupString(scope, "tauLim");
    stream.str(str);
    for (int i = 0; i < 18; i++) stream >> tauLim(i);
    stream.clear();

    // -- Gains
    mKpEE(0, 0) = cfg->lookupFloat(scope, "KpEE");
    mKpEE(1, 1) = mKpEE(0, 0);
    mKpEE(2, 2) = mKpEE(0, 0);
    mKvEE(0, 0) = cfg->lookupFloat(scope, "KvEE");
    mKvEE(1, 1) = mKvEE(0, 0);
    mKvEE(2, 2) = mKvEE(0, 0);
    mKpOr(0, 0) = cfg->lookupFloat(scope, "KpOr");
    mKpOr(1, 1) = mKpOr(0, 0);
    mKpOr(2, 2) = mKpOr(0, 0);
    mKvOr(0, 0) = cfg->lookupFloat(scope, "KvOr");
    mKvOr(1, 1) = mKvOr(0, 0);
    mKvOr(2, 2) = mKvOr(0, 0);
    mKpCOM = cfg->lookupFloat(scope, "KpCOM");
    mKvCOM = cfg->lookupFloat(scope, "KvCOM");
    mKvSpeedReg = cfg->lookupFloat(scope, "KvSpeedReg");
    mKpPose = cfg->lookupFloat(scope, "KpPose");
    mKvPose = cfg->lookupFloat(scope, "KvPose");

    // -- Weights
    // Right Arm
    if (mWaistLocked)
      str = cfg->lookupString(scope, "wEERWaistLocked");
    else
      str = cfg->lookupString(scope, "wEER");
    stream.str(str);
    for (int i = 0; i < 3; i++) stream >> mWEER(i, i);
    stream.clear();
    mWOrR = cfg->lookupFloat(scope, "wOrR");

   // // /* ********************* RES-CLF Joint Reading Parameter ************** */
    mJoint1_Pose = cfg->lookupFloat(scope, "Joint1");
    mJoint2_Pose = cfg->lookupFloat(scope, "Joint2");
    mJoint3_Pose = cfg->lookupFloat(scope, "Joint3");
    mJoint4_Pose = cfg->lookupFloat(scope, "Joint4");
    mJoint5_Pose = cfg->lookupFloat(scope, "Joint5");
    mGamma = cfg->lookupFloat(scope, "Gamma");
    mRelaxation = cfg->lookupFloat(scope, "Relaxation");
    //mJoint_SpeedReg = cfg->lookupFloat(scope, "Joint_SpeedReg");
    std::cout << "Gamma = " << mGamma << std::endl;
    std::cout << "mRelaxation = " << mRelaxation << std::endl;
    std::cout << "Controller Constructor Check 6" << std::endl;
    //mSpeedReg = cfg->lookupBoolean(scope, "SpeedReg");
    mFullSpeedReg = cfg->lookupBoolean(scope, "FullSpeedReg");

    std::cout << "mJoint1_Pose = " << mJoint1_Pose << std::endl;
    std::cout << "mJoint2_Pose = " << mJoint2_Pose << std::endl;
    std::cout << "mJoint3_Pose = " << mJoint3_Pose << std::endl;
    std::cout << "mJoint4_Pose = " << mJoint4_Pose << std::endl;
    std::cout << "mJoint5_Pose = " << mJoint5_Pose << std::endl;
    //std::cout << "mJoint_SpeedReg = " << mJoint_SpeedReg << std::endl;

    if(mFullSpeedReg) std::cout << "Doing Full Speed Regulation" << std::endl;
    
    // Left Arm
    if (mWaistLocked)
      str = cfg->lookupString(scope, "wEELWaistLocked");
    else
      str = cfg->lookupString(scope, "wEEL");
    stream.str(str);
    for (int i = 0; i < 3; i++) stream >> mWEEL(i, i);
    stream.clear();
    mWOrL = cfg->lookupFloat(scope, "wOrL");

    // Balance
    str = cfg->lookupString(scope, "wBal");
    stream.str(str);
    for (int i = 0; i < 3; i++) stream >> mWBal(i, i);
    stream.clear();
    // Regulation
    const char* s[] = {"wRegBase", "wRegWaist", "wRegTorso", "wRegKinect",
                       "wRegArm1", "wRegArm2",  "wRegArm3",  "wRegArm4",
                       "wRegArm5", "wRegArm6",  "wRegArm7"};
    for (int i = 0; i < 11; i++) {
      str = cfg->lookupString(scope, s[i]);
      stream.str(str);
      stream >> mWMatPose(i, i);
      stream >> mWMatSpeedReg(i, i);
      stream >> mWMatReg(i, i);
      if (i > 3) {
        mWMatPose(i + 7, i + 7) = mWMatPose(i, i);
        mWMatSpeedReg(i + 7, i + 7) = mWMatSpeedReg(i, i);
        mWMatReg(i + 7, i + 7) = mWMatReg(i, i);
      }
      stream.clear();
    }

  } catch (const config4cpp::ConfigurationException& ex) {
    std::cerr << ex.c_str() << std::endl;
    cfg->destroy();
  }
  std::cout << "COMAngleControl: " << (mCOMAngleControl ? "true" : "false")
            << std::endl;
  std::cout << "maintainInitCOMDistance: "
            << (mMaintainInitCOMDistance ? "true" : "false") << std::endl;
  std::cout << "tauLim: " << mTauLim.transpose() << std::endl;
  std::cout << "KpEE: " << mKpEE(0, 0) << ", " << mKpEE(1, 1) << ", "
            << mKpEE(2, 2) << std::endl;
  std::cout << "KvEE: " << mKvEE(0, 0) << ", " << mKvEE(1, 1) << ", "
            << mKvEE(2, 2) << std::endl;
  std::cout << "KpCOM: " << mKpCOM << std::endl;
  std::cout << "KvCOM: " << mKvCOM << std::endl;
  std::cout << "KvSpeedReg: " << mKvSpeedReg << std::endl;
  std::cout << "KpPose: " << mKpPose << std::endl;
  std::cout << "KvPose: " << mKvPose << std::endl;
  std::cout << "wEER: " << mWEER.diagonal().transpose() << std::endl;
  std::cout << "wEEL: " << mWEEL.diagonal().transpose() << std::endl;
  // std::cout << "wBal: " << mWBal(0, 0) << ", " << mWBal(1, 1) << ", " <<
  // mWBal(2, 2) << std::endl;
  std::cout << "wBal: " << mWBal.diagonal().transpose() << std::endl;
  std::cout << "wMatPoseReg: ";
  for (int i = 0; i < 18; i++) std::cout << mWMatPose(i, i) << ", ";
  std::cout << std::endl;
  std::cout << "wMatSpeedReg: ";
  for (int i = 0; i < 18; i++) std::cout << mWMatSpeedReg(i, i) << ", ";
  std::cout << std::endl;
  std::cout << "wMatReg: ";
  for (int i = 0; i < 18; i++) std::cout << mWMatReg(i, i) << ", ";
  std::cout << std::endl;
  std::cout << "waistLocked: " << (mWaistLocked ? "true" : "false")
            << std::endl;
  cfg->destroy();

  // PBal and bBal size based on mCOMAngleControl
  if (mCOMAngleControl) {
    mPBal = Eigen::MatrixXd::Zero(1, 18);
    mbBal = Eigen::VectorXd::Zero(1);
  } else {
    mPBal = Eigen::MatrixXd::Zero(3, 18);
    mbBal = Eigen::VectorXd::Zero(3);
  }

  // *********************************** Transform Jacobians
  mJtf.topRightCorner(8, 17) = Eigen::Matrix<double, 8, 17>::Zero();
  mJtf.bottomLeftCorner(17, 3) = Eigen::Matrix<double, 17, 3>::Zero();
  mJtf.bottomRightCorner(17, 17) = Eigen::Matrix<double, 17, 17>::Identity();
  mdJtf.setZero();

  // ******************************** zero Cols
  mZeroCol.setZero();
  mZero7Col.setZero();

//**********************************Set RES-CLF Parameters **************

  // mF.setZero();
  // mG.setZero();
  if(mFullSpeedReg){
    mP = readMatrix("/home/krang/SethResearch/13EE2-3DUnification-UnlockedJoints-Murtaza/examples/3dofddp/P_space_SpeedReg_Full.txt");
    mF = readMatrix("/home/krang/SethResearch/13EE2-3DUnification-UnlockedJoints-Murtaza/examples/3dofddp/F_space_SpeedReg_Full.txt");
    mG = readMatrix("/home/krang/SethResearch/13EE2-3DUnification-UnlockedJoints-Murtaza/examples/3dofddp/G_space_SpeedReg_Full.txt"); 
  }
  size_t n = mF.rows();
  size_t n1 = mG.cols();
  mEta = Eigen::VectorXd::Zero(n);
  // **************************** if waist locked, dimesion of decision variable
  // in QP should be reduced by one
  if (mWaistLocked)
    mOptDim = 17;
  else
    mOptDim = 18;
  mddqBodyRef = Eigen::VectorXd::Zero(mOptDim);
  mMM = Eigen::MatrixXd::Zero(mOptDim, mOptDim);
  mhh = Eigen::VectorXd::Zero(mOptDim);
  mTauLim = Eigen::VectorXd::Zero(mOptDim);
  mTauLim << tauLim(0), tauLim.tail(mOptDim - 1);

  mddqBodyRef1 = Eigen::VectorXd::Zero(n1+1);
}

//==============================================================================
Controller::~Controller() {}

//==============================================================================
struct OptParams {
  Eigen::MatrixXd P;
  Eigen::VectorXd b;
};

//==============================================================================
void printMatrix(Eigen::MatrixXd A) {
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      std::cout << A(i, j) << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//==============================================================================
void constraintFunc(unsigned m, double* result, unsigned n, const double* x,
                    double* grad, void* f_data) {
  OptParams* constParams = reinterpret_cast<OptParams*>(f_data);
  // std::cout << "done reading optParams " << std::endl;

  if (grad != NULL) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        grad[i * n + j] = constParams->P(i, j);
      }
    }
  }
  // std::cout << "done with gradient" << std::endl;

  Eigen::MatrixXd X = Eigen::VectorXd::Zero(n);
  for (size_t i = 0; i < n; i++) X(i) = x[i];
  // std::cout << "done reading x" << std::endl;

  Eigen::VectorXd mResult;
  mResult = constParams->P * X - constParams->b;
  for (size_t i = 0; i < m; i++) {
    result[i] = mResult(i);
  }
  // std::cout << "done calculating the result"
}

//==============================================================================
double optFunc(const std::vector<double>& x, std::vector<double>& grad,
               void* my_func_data) {
  OptParams* optParams = reinterpret_cast<OptParams*>(my_func_data);
  // std::cout << "done reading optParams " << std::endl;
  // Eigen::Matrix<double, 18, 1> X(x.data());
  size_t n = x.size();
  Eigen::VectorXd X = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; i++) X(i) = x[i];
  // std::cout << "done reading x" << std::endl;

  if (!grad.empty()) {
    Eigen::MatrixXd mGrad =
        optParams->P.transpose() * (optParams->P * X - optParams->b);
    // std::cout << "done calculating gradient" << std::endl;
    Eigen::VectorXd::Map(&grad[0], mGrad.size()) = mGrad;
    // std::cout << "done changing gradient cast" << std::endl;
  }
  // std::cout << "about to return something" << std::endl;
  return (0.5 * pow((optParams->P * X - optParams->b).norm(), 2));
}

// // /* ***************************************Optimization Function for RESCLF ************* */
struct OptParams_RESCLF {
  Eigen::MatrixXd L_G;
  Eigen::VectorXd L_F;
  Eigen::VectorXd V_x;
  double gamma;
  double relaxation;
};

double optFunc_RESCLF(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data) {
  OptParams_RESCLF* optParams_RESCLF = reinterpret_cast<OptParams_RESCLF*>(my_func_data);
  // std::cout << "done reading optParams " << std::endl;
  // Eigen::Matrix<double, 18, 1> X(x.data());
  size_t n = x.size();
  Eigen::VectorXd X = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; i++) X(i) = x[i];
  // std::cout << "done reading x" << std::endl;

  if (!grad.empty()) {
    Eigen::MatrixXd mGrad = 2 * X;
    mGrad(n-1,0) = 2*optParams_RESCLF->relaxation*X(n-1);
    // std::cout << "done calculating gradient in optFunc_RESCLF" << std::endl;
    Eigen::VectorXd::Map(&grad[0], mGrad.size()) = mGrad;
    // std::cout << "done changing gradient cast in optFunc_RESCLF" << std::endl;
  }
  // std::cout << "about to return something from optFunc_RESCLF" << std::endl;
  double output = 0;
    for(int i = 0; i < n-1;i++){
      output = output + pow(X(i),2);
    }
    output = output + optParams_RESCLF->relaxation*pow(X(n-1),2);
    // std::cout << "Returning output" << std::endl;
  return output;
}

double constraintFunc_RESCLF1(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  OptParams_RESCLF* constParams = reinterpret_cast<OptParams_RESCLF*>(my_func_data);
  
  // std::cout << "size of constParams->L_G = " << constParams->L_G.rows() << "*" << constParams->L_G.cols() << std::endl;
  // double gamma = optParams->gamma;
  size_t n = x.size();
  // Eigen::Matrix<double, 8, 1> X(x.data());
  Eigen::VectorXd X = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; i++) X(i) = x[i];

  // std::cout << "About to compute Gradient" << std::endl;
  if (!grad.empty()) {
    Eigen::MatrixXd mGrad = Eigen::VectorXd::Zero(n);
    
    for(int i = 0;i < n-1 ;i++){
      mGrad(i,0) = constParams->L_G(0,i);
    }
    
    mGrad(n-1,0) = -1;
    Eigen::VectorXd::Map(&grad[0], mGrad.size()) = mGrad;
    
  }
   // std::cout << "Gradient Done" << std::endl;  
  
  Eigen::Matrix<double,1 ,1> mResult;
  mResult = constParams->L_G * X.head(n-1) + constParams->L_F + constParams->gamma * constParams->V_x - X.segment<1>(n-1);
  // mResult = mResult.col(0) - X(0);
  
  double result;
  result = mResult(0,0);
  return result;

}

void constraintFunc_RESCLF2(unsigned m, double* result, unsigned n, const double* x,
                    double* grad, void* f_data) {
  OptParams* constParams = reinterpret_cast<OptParams*>(f_data);
  // std::cout << "done reading optParams in constraintFunc_RESCLF2" << std::endl;

  // std::cout << "value of m = " << m << "\nvalue of n = " << n << std::endl;
  if (grad != NULL) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n-1; j++) {
        grad[i * n + j] = constParams->P(i, j);
      }
    }
  }
  // std::cout << "done with gradient" << std::endl;

  Eigen::VectorXd X = Eigen::VectorXd::Zero(n);
  for (size_t i = 0; i < n; i++) X(i) = x[i];
  // std::cout << "done reading x" << std::endl;
  
  Eigen::VectorXd mResult;
  mResult = constParams->P * X.head(n-1) - constParams->b;//X.block<17-7,1>(0,0) - constParams->b;
  for (size_t i = 0; i < m; i++) {
    result[i] = mResult(i);
  }
  // std::cout << "done calculating the result" << std::endl;
}



//==============================================================================
void Controller::updatePositions() {
  mBaseTf = mRobot->getBodyNode(0)->getTransform().matrix();
  mq = mRobot->getPositions();
  mxyz0 = mq.segment(3, 3);  // position of frame 0 in the world frame
                             // represented in the world frame
  mpsi = atan2(mBaseTf(0, 0), -mBaseTf(1, 0));
  mqBody1 = atan2(mBaseTf(0, 1) * cos(mpsi) + mBaseTf(1, 1) * sin(mpsi),
                  mBaseTf(2, 1));
  mqBody(0) = mqBody1;
  mqBody.tail(17) = mq.tail(17);
  mRot0 << cos(mpsi), sin(mpsi), 0, -sin(mpsi), cos(mpsi), 0, 0, 0, 1;
}

//==============================================================================
void Controller::updateSpeeds() {
  mdqFilt->AddSample(mRobot->getVelocities());
  mdq = mdqFilt->average;
  mdxyz0 = mBaseTf.matrix().block<3, 3>(0, 0) *
           mdq.segment(3, 3);  // velocity of frame 0 in the world frame
                               // represented in the world frame
  mdx = mdq(4) * sin(mqBody1) - mdq(5) * cos(mqBody1);
  mdqBody1 = -mdq(0);
  mdpsi = (mBaseTf.block<3, 3>(0, 0) * mdq.head(3))(2);
  mdqBody(0) = mdqBody1;
  mdqBody.tail(17) = mdq.tail(17);
  mdqMin(0) = mdx;
  mdqMin(1) = mdpsi;
  mdqMin.tail(18) = mdqBody;
  mdRot0 << (-sin(mpsi) * mdpsi), (cos(mpsi) * mdpsi), 0, (-cos(mpsi) * mdpsi),
      (-sin(mpsi) * mdpsi), 0, 0, 0, 0;
}

//==============================================================================
void Controller::updateTransformJacobian() {
  // ********************************* Transform Jacobian
  // Coordinate Transformation to minimum set of coordinates
  // dq0 = -dq_1
  // dq1 = dpsi*cos(q_1)
  // dq2 = dpsi*sin(q_1)
  // dq3 = 0
  // dq4 = dx*sin(q_1)
  // dq5 = -dx*cos(q_1)
  // dq6 = dx/R - (L/(2*R))*dpsi - dq_1
  // dq7 = dx/R + (L/(2*R))*dpsi - dq_1
  // dq8 = dq_2
  // dq9 = dq_3
  // [dq0 dq1 dq2 dq3 dq4 dq5 dq6 dq7]' = J*[dx dpsi dq_1]';
  // where

  mJtf.topLeftCorner(8, 3) << 0, 0, -1, 0, cos(mqBody1), 0, 0, sin(mqBody1), 0,
      0, 0, 0, sin(mqBody1), 0, 0, -cos(mqBody1), 0, 0, 1 / mR, -mL / (2 * mR),
      -1, 1 / mR, mL / (2 * mR), -1;

  mdJtf.topLeftCorner(8, 3) << 0, 0, 0, 0, -sin(mqBody1) * mdqBody1, 0, 0,
      cos(mqBody1) * mdqBody1, 0, 0, 0, 0, cos(mqBody1) * mdqBody1, 0, 0,
      sin(mqBody1) * mdqBody1, 0, 0, 0, 0, 0, 0, 0, 0;

  if (mSteps < 0) {
    std::cout << "Jtf: " << std::endl;
    for (int i = 0; i < mJtf.rows(); i++) {
      for (int j = 0; j < mJtf.cols(); j++) std::cout << mJtf(i, j) << ", ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "dJtf: " << std::endl;
    for (int i = 0; i < mdJtf.rows(); i++) {
      for (int j = 0; j < mdJtf.cols(); j++) std::cout << mdJtf(i, j) << ", ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

//==============================================================================
void Controller::setLeftArmOptParams(
    const Eigen::Vector3d& _LeftTargetPosition) {
  static Eigen::Vector3d xEELref, xEEL, dxEEL, ddxEELref;
  static Eigen::Matrix<double, 3, 15> JEEL_small, dJEEL_small;
  static Eigen::Matrix<double, 3, 25> JEEL_full, dJEEL_full;
  static Eigen::Matrix<double, 3, 18> JEEL, dJEEL;

  xEELref = _LeftTargetPosition;
  if (mSteps == 1) {
    std::cout << "xEELref: " << xEELref(0) << ", " << xEELref(1) << ", "
              << xEELref(2) << std::endl;
  }

  // x, dx, ddxref
  xEEL = mRot0 * (mLeftEndEffector->getTransform().translation() - mxyz0);
  dxEEL = mRot0 * (mLeftEndEffector->getLinearVelocity() - mdxyz0) +
          mdRot0 * (mLeftEndEffector->getTransform().translation() - mxyz0);
  ddxEELref = -mKpEE * (xEEL - xEELref) - mKvEE * dxEEL;

  // Jacobian
  JEEL_small = mLeftEndEffector->getLinearJacobian();
  JEEL_full << JEEL_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      JEEL_small.block<3, 2>(0, 6), mZeroCol, JEEL_small.block<3, 7>(0, 8),
      mZero7Col;
  JEEL = (mRot0 * JEEL_full * mJtf).topRightCorner(3, 18);

  // Jacobian Derivative
  dJEEL_small = mLeftEndEffector->getLinearJacobianDeriv();
  dJEEL_full << dJEEL_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      dJEEL_small.block<3, 2>(0, 6), mZeroCol, dJEEL_small.block<3, 7>(0, 8),
      mZero7Col;
  dJEEL = (mdRot0 * JEEL_full * mJtf + mRot0 * dJEEL_full * mJtf +
           mRot0 * JEEL_full * mdJtf)
              .topRightCorner(3, 18);

  // P and b
  mPEEL << mWEEL * JEEL;
  mbEEL = -mWEEL * (dJEEL * mdqBody - ddxEELref);


//  /* Parameters for RESCLF */
  mEEL << xEEL - xEELref,dxEEL;
  mJEEL = JEEL;
  mdJEEL = dJEEL;

  mdxEEL = dxEEL;
  mddxEEL = ddxEELref;
}

//==============================================================================
void Controller::setRightArmOptParams(
    const Eigen::Vector3d& _RightTargetPosition) {
  static Eigen::Vector3d xEERref, xEER, dxEER, ddxEERref;
  static Eigen::Matrix<double, 3, 15> JEER_small, dJEER_small;
  static Eigen::Matrix<double, 3, 25> JEER_full, dJEER_full;
  static Eigen::Matrix<double, 3, 18> JEER, dJEER;

  // x, dx, ddxref
  xEERref = _RightTargetPosition;
  if (mSteps == 1) {
    std::cout << "xEErefR: " << xEERref(0) << ", " << xEERref(1) << ", "
              << xEERref(2) << std::endl;
  }
  xEER = mRot0 * (mRightEndEffector->getTransform().translation() - mxyz0);
  dxEER = mRot0 * (mRightEndEffector->getLinearVelocity() - mdxyz0) +
          mdRot0 * (mRightEndEffector->getTransform().translation() - mxyz0);
  ddxEERref = -mKpEE * (xEER - xEERref) - mKvEE * dxEER;

  // Jacobian
  JEER_small = mRightEndEffector->getLinearJacobian();
  JEER_full << JEER_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      JEER_small.block<3, 2>(0, 6), mZeroCol, mZero7Col,
      JEER_small.block<3, 7>(0, 8);
  JEER = (mRot0 * JEER_full * mJtf).topRightCorner(3, 18);

  // Jacobian Derivative
  dJEER_small = mRightEndEffector->getLinearJacobianDeriv();
  dJEER_full << dJEER_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      dJEER_small.block<3, 2>(0, 6), mZeroCol, mZero7Col,
      dJEER_small.block<3, 7>(0, 8);
  dJEER = (mdRot0 * JEER_full * mJtf + mRot0 * dJEER_full * mJtf +
           mRot0 * JEER_full * mdJtf)
              .topRightCorner(3, 18);

  // P and b
  mPEER << mWEER * JEER;
  mbEER = -mWEER * (dJEER * mdqBody - ddxEERref);

  //  /* Parameters for RESCLF */
  mEER << xEER - xEERref,dxEER;
  mJEER = JEER;
  mdJEER = dJEER;

  mdxEER = dxEER;
  mddxEER = ddxEERref;

}

//==============================================================================
void Controller::setLeftOrientationOptParams(
    const Eigen::Vector3d& _LeftTargetRPY) {
  static Eigen::Quaterniond quatRef, quat;
  static double quatRef_w, quat_w;
  static Eigen::Vector3d quatRef_xyz, quat_xyz, quatError_xyz, w, dwref;
  static Eigen::Matrix<double, 3, 15> JwL_small, dJwL_small;
  static Eigen::Matrix<double, 3, 25> JwL_full, dJwL_full;
  static Eigen::Matrix<double, 3, 18> JwL, dJwL;

  // Reference orientation (TargetRPY is assumed to be in Frame 0)
  quatRef = Eigen::Quaterniond(
      Eigen::AngleAxisd(_LeftTargetRPY(0), Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(_LeftTargetRPY(1), Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(_LeftTargetRPY(2), Eigen::Vector3d::UnitZ()));
  quatRef_w = quatRef.w();
  quatRef_xyz << quatRef.x(), quatRef.y(), quatRef.z();
  if (quatRef_w < 0) {
    quatRef_w *= -1.0;
    quatRef_xyz *= -1.0;
  }

  // Current orientation in Frame 0
  Eigen::Vector3d currentRPY = dart::math::matrixToEulerXYZ(
      mRot0 * mLeftEndEffector->getTransform().rotation());
  quat = Eigen::Quaterniond(
      Eigen::AngleAxisd(currentRPY(0), Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(currentRPY(1), Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(currentRPY(2), Eigen::Vector3d::UnitZ()));
  // quat =
  // Eigen::Quaterniond(mRot0*mLeftEndEffector->getTransform().rotation());
  quat_w = quat.w();
  quat_xyz << quat.x(), quat.y(), quat.z();
  if (pow(-quat_w - quatRef_w, 2) + pow((-quat_xyz - quatRef_xyz).norm(), 2) <
      pow(quat_w - quatRef_w, 2) + pow((-quat_xyz - quatRef_xyz).norm(), 2)) {
    quat_w *= -1.0;
    quat_xyz *= -1.0;
  }

  // Orientation error
  quatError_xyz =
      quatRef_w * quat_xyz - quat_w * quatRef_xyz + quatRef_xyz.cross(quat_xyz);

  // Jacobian
  JwL_small = mLeftEndEffector->getAngularJacobian();
  JwL_full << JwL_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      JwL_small.block<3, 2>(0, 6), mZeroCol, JwL_small.block<3, 7>(0, 8),
      mZero7Col;
  JwL = (mRot0 * JwL_full * mJtf).topRightCorner(3, 18);

  // Jacobian Derivative
  dJwL_small = mLeftEndEffector->getAngularJacobianDeriv();
  dJwL_full << dJwL_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      dJwL_small.block<3, 2>(0, 6), mZeroCol, dJwL_small.block<3, 7>(0, 8),
      mZero7Col;
  dJwL = (mdRot0 * JwL_full * mJtf + mRot0 * dJwL_full * mJtf +
          mRot0 * JwL_full * mdJtf)
             .topRightCorner(3, 18);

  // Current angular speed in frame 0 and Reference angular acceleration of the
  // end-effector in frame 0
  w = JwL * mdqBody;
  dwref = -mKpOr * quatError_xyz - mKvOr * w;

  // P and b
  mPOrL = mWOrL * JwL;
  mbOrL = -mWOrL * (dJwL * mdqBody - dwref);

   // Parameter for RESCLF
  mOrL << quatError_xyz, w;
  mJwL = JwL;
  mdJwL = dJwL;
  mOrLRef = dwref;

  mdxOrL = w;
  mddxOrL = dwref;

}

//==============================================================================
void Controller::setRightOrientationOptParams(
    const Eigen::Vector3d& _RightTargetRPY) {
  static Eigen::Quaterniond quatRef, quat;
  static double quatRef_w, quat_w;
  static Eigen::Vector3d quatRef_xyz, quat_xyz, quatError_xyz, w, dwref;
  static Eigen::Matrix<double, 3, 15> JwR_small, dJwR_small;
  static Eigen::Matrix<double, 3, 25> JwR_full, dJwR_full;
  static Eigen::Matrix<double, 3, 18> JwR, dJwR;

  // Reference orientation (TargetRPY is assumed to be in Frame 0)
  quatRef = Eigen::Quaterniond(
      Eigen::AngleAxisd(_RightTargetRPY(0), Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(_RightTargetRPY(1), Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(_RightTargetRPY(2), Eigen::Vector3d::UnitZ()));
  quatRef_w = quatRef.w();
  quatRef_xyz << quatRef.x(), quatRef.y(), quatRef.z();
  if (quatRef_w < 0) {
    quatRef_w *= -1.0;
    quatRef_xyz *= -1.0;
  }

  // Current orientation in Frame 0
  Eigen::Vector3d currentRPY = dart::math::matrixToEulerXYZ(
      mRot0 * mRightEndEffector->getTransform().rotation());
  quat = Eigen::Quaterniond(
      Eigen::AngleAxisd(currentRPY(0), Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(currentRPY(1), Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(currentRPY(2), Eigen::Vector3d::UnitZ()));
  // quat =
  // Eigen::Quaterniond(mRot0*mRightEndEffector->getTransform().rotation());
  quat_w = quat.w();
  quat_xyz << quat.x(), quat.y(), quat.z();
  if (pow(-quat_w - quatRef_w, 2) + pow((-quat_xyz - quatRef_xyz).norm(), 2) <
      pow(quat_w - quatRef_w, 2) + pow((-quat_xyz - quatRef_xyz).norm(), 2)) {
    quat_w *= -1.0;
    quat_xyz *= -1.0;
  }

  // Orientation error
  quatError_xyz =
      quatRef_w * quat_xyz - quat_w * quatRef_xyz + quatRef_xyz.cross(quat_xyz);

  // Jacobian
  JwR_small = mRightEndEffector->getAngularJacobian();
  JwR_full << JwR_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      JwR_small.block<3, 2>(0, 6), mZeroCol, mZero7Col,
      JwR_small.block<3, 7>(0, 8);
  JwR = (mRot0 * JwR_full * mJtf).topRightCorner(3, 18);

  // Jacobian Derivative
  dJwR_small = mRightEndEffector->getAngularJacobianDeriv();
  dJwR_full << dJwR_small.block<3, 6>(0, 0), mZeroCol, mZeroCol,
      dJwR_small.block<3, 2>(0, 6), mZeroCol, mZero7Col,
      dJwR_small.block<3, 7>(0, 8);
  dJwR = (mdRot0 * JwR_full * mJtf + mRot0 * dJwR_full * mJtf +
          mRot0 * JwR_full * mdJtf)
             .topRightCorner(3, 18);

  // Current angular speed in frame 0 and Reference angular acceleration of the
  // end-effector in frame 0
  w = JwR * mdqBody;
  dwref = -mKpOr * quatError_xyz - mKvOr * w;

  // P and b
  mPOrR = mWOrR * JwR;
  mbOrR = -mWOrR * (dJwR * mdqBody - dwref);


// Parameter for RESCLF
  mOrR << quatError_xyz, w;
  mJwR = JwR;
  mdJwR = dJwR;
  mOrRRef = dwref;

  mdxOrR = w;
  mddxOrR = dwref;

}

//==============================================================================
void Controller::setBalanceOptParams(double thref, double dthref,
                                     double ddthref) {
  static Eigen::Vector3d COM, dCOM, COMref, dCOMref, ddCOMref, ddCOMStar, COMActual, dCOMActual, ddCOMStar1;
  static Eigen::Matrix<double, 3, 25> JCOM_full, dJCOM_full;
  static Eigen::Matrix<double, 3, 18> JCOM, dJCOM;
  static Eigen::Matrix<double, 1, 18> Jth, dJth;
  static Eigen::Matrix<double, 1, 3> thVec, dthVec;
  static Eigen::Matrix<double, 2, 18> JCOM_Modified, dJCOM_Modified;
  static double L, th, th_wrong, dth, ddthStar;

  //*********************************** Balance
  // Excluding wheels from COM Calculation
  // Eigen::Vector3d bodyCOM = Rot0*(mRobot->getCOM() - xyz0);
  // Eigen::Vector3d bodyCOMLinearVelocity =
  // Rot0*(mRobot->getCOMLinearVelocity() - dxyz0) + dRot0*(mRobot->getCOM() -
  // xyz0);

  // COM, dCOM, JCOM and dJCOM
  COM = mRot0 * (mRobot->getCOM() - mxyz0);
  dCOM = mRot0 * (mRobot->getCOMLinearVelocity() - mdxyz0) +
         mdRot0 * (mRobot->getCOM() - mxyz0);

  // x, dx, ddxStar
  COM = mRot0 * (mRobot->getCOM() - mxyz0);
  dCOM = mRot0 * (mRobot->getCOMLinearVelocity() - mdxyz0) +
         mdRot0 * (mRobot->getCOM() - mxyz0);
  if (mCOMAngleControl) {
    th = atan2(COM(0), COM(2));
    dth = (cos(th) / COM(2)) * (cos(th) * dCOM(0) - sin(th) * dCOM(2));
    ddthStar = ddthref - mKpCOM * (th - thref) - mKvCOM * (dth - dthref);
  } else {
    if (mMaintainInitCOMDistance)
      L = mInitCOMDistance;
    else
      L = pow(COM(0) * COM(0) + COM(2) * COM(2), 0.5);
    COMref << L * sin(thref), 0, L * cos(thref);
    dCOMref << (L * cos(thref) * dthref), 0.0, (-L * sin(thref) * dthref);
    ddCOMStar << (-L * sin(thref) * dthref * dthref + L * cos(thref) * ddthref),
        0.0, (-L * cos(thref) * dthref * dthref - L * sin(thref) * ddthref);

    th = atan2(COM(0), COM(2));
    dth = (cos(th) / COM(2)) * (cos(th) * dCOM(0) - sin(th) * dCOM(2));

    COMActual << L * sin(th), 0, L * cos(th);
    dCOMActual << (L * cos(th) * dth), 0.0, (-L * sin(th) * dth);

  }

  // Jacobian
  JCOM_full = mRobot->getCOMLinearJacobian();
  JCOM = (mRot0 * JCOM_full * mJtf).topRightCorner(3, 18);
  if (mCOMAngleControl) {
    thVec << cos(th), 0.0, -sin(th);
    Jth = (cos(th) * thVec * JCOM) / COM(2);
  }

  // Jacobian derivative
  dJCOM_full = mRobot->getCOMLinearJacobianDeriv();
  dJCOM = (mdRot0 * JCOM_full * mJtf + mRot0 * dJCOM_full * mJtf +
           mRot0 * JCOM_full * mdJtf)
              .topRightCorner(3, 18);
  if (mCOMAngleControl) {
    dthVec << -sin(th), 0.0, -cos(th);
    dJth = (-sin(th) * thVec * JCOM * dth + cos(th) * dthVec * JCOM * dth +
            cos(th) * thVec * dJCOM - dCOM(2) * Jth) /
           COM(2);
  }

  // P and b
  if (mCOMAngleControl) {
    mPBal << mWBal(0, 0) * Jth;
    mbBal << mWBal(0, 0) * (-dJth * mdqBody + ddthStar);
  } else {
    mPBal << mWBal * JCOM;
    mbBal << mWBal * (-dJCOM * mdqBody + ddCOMStar);
  }
  JCOM_Modified << JCOM.block<1,18>(0,0),JCOM.block<1,18>(2,0);
  dJCOM_Modified << dJCOM.block<1,18>(0,0),dJCOM.block<1,18>(2,0);


  // mCOM_eta << COM-COMref, dCOM - dCOMref;
  mCOM_eta << COMActual-COMref, dCOMActual - dCOMref;
  mJCOM = JCOM;
  mdJCOM = dJCOM;
  mCOMRef = COMref;
  mdCOMRef = dCOMref;
  mddCOMRef = ddCOMStar;


}

//==============================================================================
void Controller::computeDynamics() {
  static Eigen::Matrix<double, 25, 25> M_full;
  static Eigen::Matrix<double, 20, 20> M;
  static Eigen::Matrix<double, 20, 1> h;
  static Eigen::Matrix<double, 19, 1> h_without_psi_equation;
  static double axx, alpha, beta;
  static Eigen::Matrix<double, 18, 1> axq, hh;
  static Eigen::Matrix<double, 18, 19> PP;
  static Eigen::Matrix<double, 18, 18> Aqq, A_qq, B, pre, MM;

  // ***************************** Inertia and Coriolis Matrices
  M_full = mRobot->getMassMatrix();
  M = mJtf.transpose() * M_full * mJtf;
  h = mJtf.transpose() * M_full * mdJtf * mdqMin +
      mJtf.transpose() * mRobot->getCoriolisAndGravityForces();
  h_without_psi_equation(0) = h(0);
  h_without_psi_equation.tail(18) = h.tail(18);

  axx = M(0, 0);
  axq = M.bottomLeftCorner(18, 1);
  Aqq = M.bottomRightCorner(18, 18);
  alpha = axq(0) / (mR * axx);
  beta = 1 / (1 + alpha);
  A_qq = Aqq - (1 / axx) * (axq * axq.transpose());  // AqqSTAR in derivation
  B << axq / (mR * axx), Eigen::Matrix<double, 18, 17>::Zero();
  pre = Eigen::Matrix<double, 18, 18>::Identity() - beta * B;
  PP << -pre * axq / axx, pre;
  MM = pre * A_qq;
  hh = PP * h_without_psi_equation;
  mMM << MM(0, 0), MM.topRightCorner(1, mOptDim - 1),
      MM.bottomLeftCorner(mOptDim - 1, 1),
      MM.bottomRightCorner(mOptDim - 1, mOptDim - 1);
  mhh << hh(0), hh.tail(mOptDim - 1);
  if (mSteps < 0) {
    std::cout << "axx: " << axx << std::endl;
    std::cout << "axq: ";
    for (int i = 0; i < axq.rows(); i++)
      for (int j = 0; j < axq.cols(); j++) std::cout << axq(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "Aqq: ";
    for (int i = 0; i < Aqq.rows(); i++)
      for (int j = 0; j < Aqq.cols(); j++) std::cout << Aqq(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "A_qq: ";
    for (int i = 0; i < A_qq.rows(); i++)
      for (int j = 0; j < A_qq.cols(); j++) std::cout << A_qq(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "B: ";
    for (int i = 0; i < B.rows(); i++)
      for (int j = 0; j < B.cols(); j++) std::cout << B(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "pre: ";
    for (int i = 0; i < pre.rows(); i++)
      for (int j = 0; j < pre.cols(); j++) std::cout << pre(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "PP: ";
    for (int i = 0; i < PP.rows(); i++)
      for (int j = 0; j < PP.cols(); j++) std::cout << PP(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "MM: ";
    for (int i = 0; i < mMM.rows(); i++)
      for (int j = 0; j < mMM.cols(); j++) std::cout << mMM(i, j) << ", ";
    std::cout << std::endl;
    std::cout << "hh: ";
    for (int i = 0; i < mhh.rows(); i++)
      for (int j = 0; j < mhh.cols(); j++) std::cout << mhh(i, j) << ", ";
    std::cout << std::endl;
  }
  std::cout << "Size of M = " << M.rows() << "*" << M.cols() << std::endl;
  std::cout << "Size of h = " << h.rows() << "*" << h.cols() << std::endl;
}

//==============================================================================
void Controller::update(const Eigen::Vector3d& _LeftTargetPosition,
                        const Eigen::Vector3d& _RightTargetPosition,
                        const Eigen::Vector3d& _LeftTargetRPY,
                        const Eigen::Vector3d& _RightTargetRPY, double thref,
                        double dthref, double ddthref, double tau_0) {
  // increase the step counter
  mSteps++;

  // updates mBaseTf, mq, mxyz0, mpsi, mqBody1, mqBody, mRot0
  // Needs mRobot
  updatePositions();
  // updates mdq, mdxyz0, mdx, mdqBody1, mdpsi, mdqBody, mdqMin, dRot0
  // Needs mRobot, mdqFilt, mBaseTf, mqBody1
  updateSpeeds();

  // updates mJtf and mdJtf
  // Needs mqBody1, mdqBody1, mR, mL
  updateTransformJacobian();

  // sets mPEEL and mbEEL
  // Needs mRot0, mLeftEndEffector, mxyz0, mdxyz0, mdRot0, mKpEE, mKvEE, mJtf,
  // mdJtf
  setLeftArmOptParams(_LeftTargetPosition);

  // sets mPEER and mbEER
  // Needs mRot0, mRightEndEffector, mxyz0, mdxyz0, mdRot0, mKpEE, mKvEE, mJtf,
  // mdJtf
  setRightArmOptParams(_RightTargetPosition);

  setLeftOrientationOptParams(_LeftTargetRPY);
  setRightOrientationOptParams(_RightTargetRPY);

  // sets mPBal and mbBal
  // Needs mRot0, mRobot, mxyz0, mdxyz0, mdRot0, mKpCOM, mKvCOM, mJtf, mdJtf
  setBalanceOptParams(thref, dthref, ddthref);

  // set Regulation Opt Params
  mPPose = mWMatPose;
  mbPose << mWMatPose * (-mKpPose * (mqBody - mqBodyInit) - mKvPose * mdqBody);

  mPSpeedReg = mWMatSpeedReg;
  mbSpeedReg << -mWMatSpeedReg * mKvSpeedReg * mdqBody;

  mPReg = mWMatReg;
  mbReg.setZero();

  // set mMM and mhh
  // Needs mRobot, mJtf, mdJtf, mdqMin, mR
  computeDynamics();

  std::cout << "Size of mdJEEL = " << mdJEEL.rows() << "*" << mdJEEL.cols() <<std::endl;
  std::cout << "Size of mdqBody = " << mdqBody.rows() << "*" << mdqBody.cols() <<std::endl;
  std::cout << "Size of mJEEL = " << mJEEL.rows() << "*" << mJEEL.cols() <<std::endl;
  std::cout << "Size of mddqBodyRef = " << mddqBodyRef.rows() << "*" << mddqBodyRef.cols() <<std::endl;
  // mEELRef = mdJEEL*mdqBody + mJEEL*mddqBodyRef;
  // mEERRef = mdJEER*mdqBody + mJEER*mddqBodyRef;
  // mOrLRef = mdJwL*mdqBody + mJwL*mddqBodyRef;
  // mOrRRef = mdJwR*mdqBody + mJwR*mddqBodyRef;
  // mCOMRef = mdJCOM*mdqBody + mJCOM*mddqBodyRef;
   mEELRef = Eigen::VectorXd::Zero(3);
   mEERRef = Eigen::VectorXd::Zero(3);
   mOrLRef = Eigen::VectorXd::Zero(3);
   mOrRRef = Eigen::VectorXd::Zero(3);

  // ***************************** QP
  OptParams optParams;
  Eigen::MatrixXd P(mPEER.rows() + mPOrR.rows() + mPEEL.rows() + mPOrL.rows() +
                        mPBal.rows() + mPPose.rows() + mPSpeedReg.rows() +
                        mPReg.rows(),
                    mOptDim);
  P << mPEER.col(0), mPEER.topRightCorner(mPEER.rows(), mOptDim - 1),
      mPOrR.col(0), mPOrR.topRightCorner(mPOrR.rows(), mOptDim - 1),
      mPEEL.col(0), mPEEL.topRightCorner(mPEEL.rows(), mOptDim - 1),
      mPOrL.col(0), mPOrL.topRightCorner(mPOrL.rows(), mOptDim - 1),
      mPBal.col(0), mPBal.topRightCorner(mPBal.rows(), mOptDim - 1),
      mPPose.col(0), mPPose.topRightCorner(mPPose.rows(), mOptDim - 1),
      mPSpeedReg.col(0),
      mPSpeedReg.topRightCorner(mPSpeedReg.rows(), mOptDim - 1), mPReg.col(0),
      mPReg.topRightCorner(mPReg.rows(), mOptDim - 1);

  Eigen::VectorXd b(mbEER.rows() + mbOrR.rows() + mbEEL.rows() + mbOrL.rows() +
                        mbBal.rows() + mbPose.rows() + mbSpeedReg.rows() +
                        mbReg.rows(),
                    mbEER.cols());
  b << mbEER, mbOrR, mbEEL, mbOrL, mbBal, mbPose, mbSpeedReg, mbReg;
  optParams.P = P;
  optParams.b = b;

  const std::vector<double> inequalityconstraintTol(mOptDim, 1e-3);
  OptParams inequalityconstraintParams[2];
  inequalityconstraintParams[0].P = mMM;
  inequalityconstraintParams[1].P = -mMM;
  inequalityconstraintParams[0].b = -mhh + mTauLim;
  inequalityconstraintParams[1].b = mhh + mTauLim;

  // nlopt::opt opt(nlopt::LN_COBYLA, 30);
  nlopt::opt opt(nlopt::LD_SLSQP, mOptDim);
  double minf;
  opt.set_min_objective(optFunc, &optParams);
  opt.add_inequality_mconstraint(constraintFunc, &inequalityconstraintParams[0],
                                 inequalityconstraintTol);
  opt.add_inequality_mconstraint(constraintFunc, &inequalityconstraintParams[1],
                                 inequalityconstraintTol);
  opt.set_xtol_rel(1e-3);
  if (maxTimeSet) opt.set_maxtime(0.1);
  std::vector<double> ddqBodyRef_vec(mOptDim);
  Eigen::VectorXd::Map(&ddqBodyRef_vec[0], mddqBodyRef.size()) = mddqBodyRef;
  try {
    nlopt::result result = opt.optimize(ddqBodyRef_vec, minf);
  } catch (std::exception& e) {
    // std::cout << "nlopt failed: " << e.what() << std::endl;
  }
  for (int i = 0; i < mOptDim; i++) mddqBodyRef(i) = ddqBodyRef_vec[i];

  Eigen::VectorXd bodyTorques_Munzir = mMM * mddqBodyRef + mhh;



// // /* ************************Computing RES-CLF************* */
  // size_t n = mF.rows();
  // std::cout << "n = mF.rows() = " << n <<  std::endl;

  // std::cout << "Check 1" <<  std::endl;
  // Eigen::VectorXd mEta(n);
  // std::cout << "Size of mEta = " << mEta.rows() << "*" << mEta.cols() <<  std::endl;
  // mEta << mEEL, mEER, mOrL, mOrR, mCOM_eta, mdqBody.block<1,1>(mJoint1_Pose-1,0), mdqBody.block<1,1>(mJoint2_Pose-1,0), mdqBody.block<1,1>(mJoint3_Pose-1,0);

  // std::cout << "Check 2" <<  std::endl;

  // Eigen::MatrixXd MM_inverse = mMM.colPivHouseholderQr().inverse();
  // Eigen::MatrixXd MM_inverse_hh = -mMM.colPivHouseholderQr().inverse()*mhh;

  // Eigen::Matrix<double, 18, 18> A_LgLfy;
  // Eigen::Matrix<double, 18,1> Lf_yx;

  // std::cout << "Check 3" <<  std::endl;

  // A_LgLfy << mJEEL*MM_inverse,
  //            mJEER*MM_inverse,
  //            mJwL*MM_inverse,
  //            mJwR*MM_inverse,
  //            mJCOM*MM_inverse, 
  //            MM_inverse.block<1,18>(mJoint1_Pose-1,0),
  //            MM_inverse.block<1,18>(mJoint2_Pose-1,0),
  //            MM_inverse.block<1,18>(mJoint3_Pose-1,0);

  // std::cout << "Check 4" <<  std::endl;

  // std::cout << "A_LgLfy = " << std::endl;
  // printMatrix(A_LgLfy);
  // std::cout << "mMM = " << std::endl;
  // printMatrix(mMM);
  // std::cout << "MM_inverse = " << std::endl;
  // printMatrix(MM_inverse);
  // std::cout << "mJEEL = " << std::endl;
  // printMatrix(mJEEL);
  // std::cout << "mJEER = " << std::endl;
  // printMatrix(mJEER);
  // std::cout << "mJwL = " << std::endl;
  // printMatrix(mJwL);
  // std::cout << "mJwR = " << std::endl;
  // printMatrix(mJwR);
  // std::cout << "mJCOM = " << std::endl;
  // printMatrix(mJCOM);
  // std::cout << " MM_inverse.block<1,17-7>(mJoint1_Pose-1,0); = " << std::endl;
  // printMatrix( MM_inverse.block<1,18>(mJoint1_Pose-1,0));
  // std::cout << " MM_inverse.block<1,18>(mJoint2_Pose-1,0); = " << std::endl;
  // printMatrix( MM_inverse.block<1,18>(mJoint2_Pose-1,0));
  // std::cout << " MM_inverse.block<1,18>(mJoint3_Pose-1,0); = " << std::endl;
  // printMatrix( MM_inverse.block<1,18>(mJoint3_Pose-1,0));

  // Lf_yx << mdJEEL*mdqBody - mJEEL*MM_inverse*mhh - mEELRef,
  //          mdJEER*mdqBody - mJEER*MM_inverse*mhh - mEERRef,
  //          mdJwL*mdqBody - mJwL*MM_inverse*mhh - mOrLRef,
  //          mdJwR*mdqBody - mJwR*MM_inverse*mhh - mOrRRef,
  //          mdJCOM*mdqBody - mJCOM*MM_inverse*mhh - mCOMRef,  
  //          MM_inverse_hh.block<1,1>(mJoint1_Pose-1,0),
  //          MM_inverse_hh.block<1,1>(mJoint2_Pose-1,0),
  //          MM_inverse_hh.block<1,1>(mJoint3_Pose-1,0);

  // std::cout << "Check 5" <<  std::endl;
  // std::cout << "The determinant of A_LgLfy is " << A_LgLfy.determinant() << std::endl;
  // std::cout << "The determinant of mMM is " << mMM.determinant() << std::endl;

  // Eigen::MatrixXd LfV_x = mEta.transpose()*(mF.transpose()*mP+mP*mF)*mEta;
  // Eigen::MatrixXd LgV_x = 2*mEta.transpose()*mP*mG;
  // Eigen::MatrixXd V_x = mEta.transpose()*mP*mEta;

  // OptParams_RESCLF optParams_RESCLF;
  // OptParams optParams1;

  // double lambda_minQ = 1;  // Provided by the Matlab QQ Matrix
  // double lambda_maxP = 2.7321; /// Provided by the Matlab P Matrix

  // OptParams_RESCLF inequalityconstraintParams_RESCLF;
  // inequalityconstraintParams_RESCLF.L_F = LfV_x;
  // inequalityconstraintParams_RESCLF.L_G = LgV_x;
  // inequalityconstraintParams_RESCLF.V_x = V_x;
  // inequalityconstraintParams_RESCLF.gamma = mGamma;//lambda_minQ/lambda_maxP;
  // inequalityconstraintParams_RESCLF.relaxation = mRelaxation;

  // std::cout << "Check 6" <<  std::endl;
  // // nlopt::opt opt1(nlopt::LN_COBYLA, mOptDim);
  // nlopt::opt opt1(nlopt::LD_SLSQP, 19);
  // // nlopt::opt opt1(nlopt::AUGLAG, 7);

  // double minf1;
  // opt1.set_min_objective(optFunc_RESCLF, &inequalityconstraintParams_RESCLF);
  // opt1.add_inequality_constraint(constraintFunc_RESCLF1, &inequalityconstraintParams_RESCLF,1e-3);
  // OptParams inequalityconstraintParams1[2];

  // const std::vector<double> inequalityconstraintTol1(18, 1e-3);
  
  // bool inververtible_check =  A_LgLfy.colPivHouseholderQr().isInvertible(); 
  // std::cout << "Is The matrix is Invertible?   " << inververtible_check << std::endl;

  // inequalityconstraintParams1[0].P = A_LgLfy.colPivHouseholderQr().inverse();
  // inequalityconstraintParams1[1].P = -A_LgLfy.colPivHouseholderQr().inverse();
  // inequalityconstraintParams1[0].b = A_LgLfy.colPivHouseholderQr().solve(Lf_yx) + mTauLim;
  // inequalityconstraintParams1[1].b = -A_LgLfy.colPivHouseholderQr().solve(Lf_yx) + mTauLim;
  // opt1.add_inequality_mconstraint(constraintFunc_RESCLF2, &inequalityconstraintParams1[0],
  //                                inequalityconstraintTol1);
  // opt1.add_inequality_mconstraint(constraintFunc_RESCLF2, &inequalityconstraintParams1[1],
  //                                inequalityconstraintTol1);


  // std::vector<double> ddqBodyRef_vec1(19);
  // Eigen::VectorXd::Map(&ddqBodyRef_vec1[0], mddqBodyRef.size()) = mddqBodyRef;
  // std::cout << "size of mddqBodyRef = " <<mddqBodyRef.rows() << "*" << mddqBodyRef.cols() <<std::endl;
  // try {
  //   std::cout << "Check 7" <<  std::endl;
  //   opt1.set_xtol_rel(1e-4);
  //   opt1.set_maxtime(1);
  //   nlopt::result result = opt1.optimize(ddqBodyRef_vec1, minf);
  // } 
  // catch (std::exception& e) {
  //   std::cout << "nlopt failed: " << e.what() << std::endl;
  // }
  // std::cout << "Check 8" <<  std::endl;
  // Eigen::Matrix<double, 18, 1> ddq1;
  // for(int i = 0;i < 18;i++) ddq1(i,0) = ddqBodyRef_vec1[i]; 

  // std::cout << "Check 9" <<  std::endl;
  // Eigen::VectorXd bodyTorques1;
  // bodyTorques1 = A_LgLfy.colPivHouseholderQr().solve(-Lf_yx + ddq1);

// // /* ************************Computing RES-CLF-Alternative ************* */    This piece of code works
   Eigen::MatrixXd MM_inverse = mMM.colPivHouseholderQr().inverse();
  // Eigen::MatrixXd MM_inverse_hh = -mMM.colPivHouseholderQr().inverse()*mhh;

  Eigen::Matrix<double,15,18> J_Aug; // = Eigen::MatrixXd::Zero(6,7);
  Eigen::Matrix<double,15,18> dJ_Aug;
  Eigen::Matrix<double,15,1> desired_ddx_aug;

  Eigen::Vector3d ddxEEL = mdJEEL * mdqBody + mJEEL * mddqBodyRef;
  Eigen::Vector3d ddxEER = mdJEER * mdqBody + mJEER * mddqBodyRef;
  Eigen::Vector3d ddxOrL = mdJwL * mdqBody + mJwL * mddqBodyRef;
  Eigen::Vector3d ddxOrR = mdJwR * mdqBody + mJwR * mddqBodyRef;
  Eigen::Vector3d ddCOMRef = mdJCOM * mdqBody + mJCOM * mddqBodyRef;


  J_Aug << mJEEL,mJEER,mJwL,mJwR,mJCOM;
  dJ_Aug << mdJEEL,mdJEER,mdJwL,mdJwR,mdJCOM;

  // desired_ddx_aug << mddxEEL,mddxEER,mddxOrL,mddxOrR,mddCOMRef;
  desired_ddx_aug << ddxEEL,ddxEER,ddxOrL,ddxOrR,ddCOMRef;

  Eigen::MatrixXd pinv_J_Aug = J_Aug.transpose()* (J_Aug * J_Aug.transpose() + 0.000025 * Eigen::MatrixXd::Identity(15,15)).colPivHouseholderQr().inverse();
  Eigen::MatrixXd M2 = mMM*pinv_J_Aug;    //n*6




  Eigen::VectorXd bodyTorques1 = M2*(desired_ddx_aug - dJ_Aug*mdqBody) + mhh;

  std::cout << "bodyTorques from RESCLF = \n" << bodyTorques1 << std::endl;
  std::cout << "bodyTorques from Munzir = \n" << bodyTorques_Munzir << std::endl;


// // /* ************************Computing RES-CLF-Augment ************* */    This piece of code donot work well.
  //  Eigen::MatrixXd MM_inverse = mMM.colPivHouseholderQr().inverse();
  // // Eigen::MatrixXd MM_inverse_hh = -mMM.colPivHouseholderQr().inverse()*mhh;

  // Eigen::MatrixXd MEEL = mJEEL*MM_inverse*mJEEL.transpose();
  // Eigen::MatrixXd MEEL_Inverse = MEEL.inverse();

  // Eigen::MatrixXd MEER = mJEER*MM_inverse*mJEER.transpose();
  // Eigen::MatrixXd MEER_Inverse = MEER.inverse();

  // Eigen::MatrixXd MOrL = mJwL*MM_inverse*mJwL.transpose();
  // Eigen::MatrixXd MOrL_Inverse = MOrL.inverse();

  // Eigen::MatrixXd MOrR = mJwR*MM_inverse*mJwR.transpose();
  // Eigen::MatrixXd MOrR_Inverse = MOrR.inverse();

  // Eigen::Matrix<double, 2, 18> JCOM_Modified;
  // JCOM_Modified.block<1,18>(0,0) = mJCOM.block<1,18>(0,0);
  // JCOM_Modified.block<1,18>(1,0) = mJCOM.block<1,18>(2,0);

  // Eigen::Matrix<double, 2, 18> dJCOM_Modified;
  // dJCOM_Modified.block<1,18>(0,0) = mdJCOM.block<1,18>(0,0);
  // dJCOM_Modified.block<1,18>(1,0) = mdJCOM.block<1,18>(2,0);

  // // Eigen::MatrixXd MCOM = mJCOM*MM_inverse*mJCOM.transpose();
  // Eigen::MatrixXd MCOM = JCOM_Modified*MM_inverse*JCOM_Modified.transpose();
  // Eigen::MatrixXd MCOM_Inverse = MCOM.inverse();

  // Eigen::Matrix<double,14,18> J_Aug; // = Eigen::MatrixXd::Zero(6,7);
  // Eigen::Matrix<double,14,18> dJ_Aug;
  // Eigen::Matrix<double,14,1> desired_ddx_aug;

  // Eigen::Vector3d ddxEEL = mdJEEL * mdqBody + mJEEL * mddqBodyRef;
  // Eigen::Vector3d ddxEER = mdJEER * mdqBody + mJEER * mddqBodyRef;
  // Eigen::Vector3d ddxOrL = mdJwL * mdqBody + mJwL * mddqBodyRef;
  // Eigen::Vector3d ddxOrR = mdJwR * mdqBody + mJwR * mddqBodyRef;
  // // Eigen::Vector3d ddCOMRef = mdJCOM * mdqBody + mJCOM * mddqBodyRef;
  // Eigen::Vector2d ddCOMRef = dJCOM_Modified * mdqBody + JCOM_Modified * mddqBodyRef;

  // // J_Aug << mJEEL,mJEER,mJwL,mJwR,mJCOM;
  // // dJ_Aug << mdJEEL,mdJEER,mdJwL,mdJwR,mdJCOM;
  // J_Aug << mJEEL,mJEER,mJwL,mJwR,JCOM_Modified;
  // dJ_Aug << mdJEEL,mdJEER,mdJwL,mdJwR,dJCOM_Modified;

  // Eigen::Matrix<double, 14,14> M_Aug = Eigen::MatrixXd::Zero(14,14);
  // M_Aug.block<3,3>(0,0) = MEEL_Inverse;
  // M_Aug.block<3,3>(3,3) = MEER_Inverse;
  // M_Aug.block<3,3>(6,6) = MOrL_Inverse;
  // M_Aug.block<3,3>(9,9) = MOrR_Inverse;
  // M_Aug.block<2,2>(12,12) = MCOM_Inverse;

  //   // desired_ddx_aug << mddxEEL,mddxEER,mddxOrL,mddxOrR,mddCOMRef;
  // // desired_ddx_aug << ddxEEL - mdJEEL*mdqBody,ddxEER - mdJEER*mdqBody,ddxOrL - mdJwL*mdqBody,ddxOrR - mdJwR*mdqBody, ddCOMRef - mdJCOM*mdqBody;
  // desired_ddx_aug << ddxEEL - mdJEEL*mdqBody,ddxEER - mdJEER*mdqBody,ddxOrL - mdJwL*mdqBody,ddxOrR - mdJwR*mdqBody, ddCOMRef - dJCOM_Modified*mdqBody;
  // // desired_ddx_aug << ddxEEL - mdJEEL*mdqBody + mJEEL*MM_inverse*mhh,ddxEER - mdJEER*mdqBody + mJEER*MM_inverse*mhh,
  // //                    ddxOrL - mdJwL*mdqBody + mJwL*MM_inverse*mhh, ddxOrR - mdJwR*mdqBody + mJwR*MM_inverse*mhh,
  // //                     ddCOMRef - mdJCOM*mdqBody + mJCOM*MM_inverse*mhh;

  // Eigen::VectorXd bodyTorques1 = J_Aug.transpose()*M_Aug*desired_ddx_aug +  mhh;
  // // Eigen::VectorXd bodyTorques1 = J_Aug.transpose()*M_Aug*desired_ddx_aug;

  // std::cout << "bodyTorques from RESCLF = \n" << bodyTorques1 << std::endl;
  // std::cout << "bodyTorques from Munzir = \n" << bodyTorques_Munzir << std::endl;



// // // /* ************************Computing I/O Feedback Linearization ************* */  //Do not work
//   Eigen::MatrixXd MM_inverse = mMM.colPivHouseholderQr().inverse();
//   Eigen::MatrixXd MM_inverse_hh = -mMM.colPivHouseholderQr().inverse()*mhh;

//   size_t n = mF.rows();
//   // std::cout << "Check 1" <<  std::endl;
//   Eigen::VectorXd mEta(n);
//   if(mSteps == 1){
//     std::cout << "n = mF.rows() = " << n <<  std::endl;
//     std::cout << "Size of mEta = " << mEta.rows() << "*" << mEta.cols() <<  std::endl;
//   }

//   Eigen::Matrix<double, 2, 18> JCOM_Modified;
//   JCOM_Modified.block<1,18>(0,0) = mJCOM.block<1,18>(0,0);
//   JCOM_Modified.block<1,18>(1,0) = mJCOM.block<1,18>(2,0);

//   Eigen::Matrix<double, 2, 18> dJCOM_Modified;
//   dJCOM_Modified.block<1,18>(0,0) = mdJCOM.block<1,18>(0,0);
//   dJCOM_Modified.block<1,18>(1,0) = mdJCOM.block<1,18>(2,0);

//   mEta << mEEL, mEER, mOrL, mOrR, mCOM_eta;

//   // std::cout << "Check 2" <<  std::endl;

//   // Eigen::Vector3d ddxEEL = mdJEEL * mdqBody + mJEEL * mddqBodyRef;
//   // Eigen::Vector3d ddxEER = mdJEER * mdqBody + mJEER * mddqBodyRef;
//   // Eigen::Vector3d ddxOrL = mdJwL * mdqBody + mJwL * mddqBodyRef;
//   // Eigen::Vector3d ddxOrR = mdJwR * mdqBody + mJwR * mddqBodyRef;
//   // Eigen::Vector3d ddCOMRef = mdJCOM * mdqBody + mJCOM * mddqBodyRef;
//   Eigen::Vector2d ddCOMRef = dJCOM_Modified * mdqBody + JCOM_Modified * mddqBodyRef;

//  std::cout << "Check 3" <<  std::endl;

//   Eigen::Matrix<double, 14, 18> A_LgLfy;
//   Eigen::Matrix<double, 14,1> Lf_yx;
//   Eigen::Matrix<double, 14,1> v;

//   A_LgLfy << mJEEL*MM_inverse,
//              mJEER*MM_inverse,
//              mJwL*MM_inverse,
//              mJwR*MM_inverse,
//              JCOM_Modified*MM_inverse;// mJCOM*MM_inverse;

//   std::cout << "Check 4" <<  std::endl;

//   // std::cout << "A_LgLfy = " << std::endl;
//   // printMatrix(A_LgLfy);
//   // std::cout << "mMM = " << std::endl;
//   // printMatrix(mMM);
//   // std::cout << "MM_inverse = " << std::endl;
//   // printMatrix(MM_inverse);
//   // std::cout << "mJEEL = " << std::endl;
//   // printMatrix(mJEEL);
//   // std::cout << "mJEER = " << std::endl;
//   // printMatrix(mJEER);
//   // std::cout << "mJwL = " << std::endl;
//   // printMatrix(mJwL);
//   // std::cout << "mJwR = " << std::endl;
//   // printMatrix(mJwR);
//   // std::cout << "mJCOM = " << std::endl;
//   // printMatrix(mJCOM);

//   Lf_yx << mdJEEL*mdqBody - mJEEL*MM_inverse*mhh,
//            mdJEER*mdqBody - mJEER*MM_inverse*mhh,
//            mdJwL*mdqBody - mJwL*MM_inverse*mhh,
//            mdJwR*mdqBody - mJwR*MM_inverse*mhh,
//            dJCOM_Modified*mdqBody - JCOM_Modified*MM_inverse*mhh;// mdJCOM*mdqBody - mJCOM*MM_inverse*mhh;

//   std::cout << "Check 5" <<  std::endl;
//   // std::cout << "The determinant of A_LgLfy is " << A_LgLfy.determinant() << std::endl;
//   std::cout << "The determinant of mMM is " << mMM.determinant() << std::endl;

//   v << ddxEEL,ddxEER, ddxOrL, ddxOrR,ddCOMRef;

//   Eigen::MatrixXd LfV_x = mEta.transpose()*(mF.transpose()*mP+mP*mF)*mEta;
//   Eigen::MatrixXd LgV_x = 2*mEta.transpose()*mP*mG;
//   Eigen::MatrixXd V_x = mEta.transpose()*mP*mEta;

//   Eigen::Matrix<double,15,18> J_Aug; // = Eigen::MatrixXd::Zero(6,7);
//   Eigen::Matrix<double,15,18> dJ_Aug;
//   Eigen::Matrix<double,15,1> desired_ddx_aug;

//   // Eigen::MatrixXd pinv_A_LgLfy = A_LgLfy.transpose()* (A_LgLfy * A_LgLfy.transpose() + 0.000025 * Eigen::MatrixXd::Identity(15,15)).colPivHouseholderQr().inverse();
//   Eigen::MatrixXd pinv_A_LgLfy = A_LgLfy.transpose()* (A_LgLfy * A_LgLfy.transpose() + 0.0025 * Eigen::MatrixXd::Identity(14,14)).colPivHouseholderQr().inverse();
//   std::cout << "The determinant of (A_LgLfy * A_LgLfy.transpose()) =   " << (A_LgLfy * A_LgLfy.transpose()).determinant() << std::endl;
//   // std::cout << "The determinant of (A_LgLfy * A_LgLfy.transpose() + 0.000025 * Eigen::MatrixXd::Identity(15,15)) =   " << (A_LgLfy * A_LgLfy.transpose() + 0.000025 * Eigen::MatrixXd::Identity(15,15)).determinant() << std::endl;
//    std::cout << "The determinant of (A_LgLfy * A_LgLfy.transpose() + 0.000025 * Eigen::MatrixXd::Identity(15,15)) =   " << (A_LgLfy * A_LgLfy.transpose() + 0.000025 * Eigen::MatrixXd::Identity(14,14)).determinant() << std::endl;
 
 
//   // desired_ddx_aug << mddxEEL,mddxEER,mddxOrL,mddxOrR,mddCOMRef;
//   // desired_ddx_aug << ddxEEL,ddxEER,ddxOrL,ddxOrR,ddCOMRef;
    
//   Eigen::VectorXd bodyTorques1 = pinv_A_LgLfy*(-Lf_yx + v);
//   Eigen::VectorXd bodyTorques2 = pinv_A_LgLfy*(-Lf_yx - v);

//   std::cout << "bodyTorques1(+desired_ddx_aug) from RESCLF = \n" << bodyTorques1 << std::endl;
//   std::cout << "bodyTorques2(-desired_ddx_aug) from RESCLF = \n" << bodyTorques2 << std::endl;
//   std::cout << "bodyTorques from Munzir = \n" << bodyTorques_Munzir << std::endl;

// // /* ************************Computing RES-CLF-Alternative with QP************* */    This piece of code donot works
//    Eigen::MatrixXd MM_inverse = mMM.colPivHouseholderQr().inverse();
//   // Eigen::MatrixXd MM_inverse_hh = -mMM.colPivHouseholderQr().inverse()*mhh;

//   Eigen::Matrix<double,15,18> J_Aug; // = Eigen::MatrixXd::Zero(6,7);
//   Eigen::Matrix<double,15,18> dJ_Aug;
//   Eigen::Matrix<double,15,1> desired_ddx_aug;
//   Eigen::Matrix<double, 3,1> Zero3Row = Eigen::MatrixXd::Zero(3,1);

//   Eigen::Vector3d ddxEEL = mdJEEL * mdqBody + mJEEL * mddqBodyRef;
//   Eigen::Vector3d ddxEER = mdJEER * mdqBody + mJEER * mddqBodyRef;
//   Eigen::Vector3d ddxOrL = mdJwL * mdqBody + mJwL * mddqBodyRef;
//   Eigen::Vector3d ddxOrR = mdJwR * mdqBody + mJwR * mddqBodyRef;
//   Eigen::Vector3d ddCOMRef = mdJCOM * mdqBody + mJCOM * mddqBodyRef;

//   J_Aug << mJEEL,mJEER,mJwL,mJwR,mJCOM;
//   dJ_Aug << mdJEEL,mdJEER,mdJwL,mdJwR,mdJCOM;

//   size_t n = mF.rows();
//   size_t n1 = mG.cols();
//   // std::cout << "Check 1" <<  std::endl;
//   Eigen::VectorXd mEta(n);
//   if(mSteps == 1){
//     std::cout << "n = mF.rows() = " << n <<  std::endl;
//     std::cout << "Size of mEta = " << mEta.rows() << "*" << mEta.cols() <<  std::endl;
//   }

//   mEta << mEEL, mEER, mOrL, mOrR, mCOM_eta;
//   // desired_ddx_aug << mddxEEL,mddxEER,mddxOrL,mddxOrR,mddCOMRef;
//   // desired_ddx_aug << ddxEEL,ddxEER,ddxOrL,ddxOrR,ddCOMRef;
//   // desired_ddx_aug << Zero3Row,Zero3Row,Zero3Row,Zero3Row, ddCOMRef;
//   desired_ddx_aug << Zero3Row,Zero3Row,Zero3Row,Zero3Row,Zero3Row ;

//   const int size_pinvJ = (J_Aug * J_Aug.transpose()).rows();

//   Eigen::MatrixXd pinv_J_Aug = J_Aug.transpose()* (J_Aug * J_Aug.transpose() + 0.000025 * Eigen::MatrixXd::Identity(15,15)).colPivHouseholderQr().inverse();
//   Eigen::MatrixXd M2 = mMM*pinv_J_Aug;    //n*6


//   Eigen::MatrixXd LfV_x = mEta.transpose()*(mF.transpose()*mP+mP*mF)*mEta;
//   Eigen::MatrixXd LgV_x = 2*mEta.transpose()*mP*mG;
//   Eigen::MatrixXd V_x = mEta.transpose()*mP*mEta;


//   OptParams_RESCLF optParams_RESCLF;
//   OptParams optParams1;

//   double lambda_minQ = 1;  // Provided by the Matlab QQ Matrix
//   double lambda_maxP = 2.7321; /// Provided by the Matlab P Matrix

//   OptParams_RESCLF inequalityconstraintParams_RESCLF;
//   inequalityconstraintParams_RESCLF.L_F = LfV_x;
//   inequalityconstraintParams_RESCLF.L_G = LgV_x;
//   inequalityconstraintParams_RESCLF.V_x = V_x;
//   inequalityconstraintParams_RESCLF.gamma = mGamma;//lambda_minQ/lambda_maxP;
//   inequalityconstraintParams_RESCLF.relaxation = mRelaxation;

//   std::cout << "Check 6" <<  std::endl;
//   // nlopt::opt opt1(nlopt::LN_COBYLA, mOptDim);
//   nlopt::opt opt1(nlopt::LD_SLSQP, n1+1);
//   // nlopt::opt opt1(nlopt::AUGLAG, 7);
//   std::cout << "Check 7" <<  std::endl;
//   double minf1;
//   opt1.set_min_objective(optFunc_RESCLF, &inequalityconstraintParams_RESCLF);
//   opt1.add_inequality_constraint(constraintFunc_RESCLF1, &inequalityconstraintParams_RESCLF,1e-3);
//   OptParams inequalityconstraintParams1[2];
// std::cout << "Check 8" <<  std::endl;
//   const std::vector<double> inequalityconstraintTol1(18, 1e-3);
//   //   inequalityconstraintParams[0].P = mMM;
//   // inequalityconstraintParams[1].P = -mMM;
//   // inequalityconstraintParams[0].b = -mhh + mTauLim;
//   // inequalityconstraintParams[1].b = mhh + mTauLim;

//   inequalityconstraintParams1[0].P = M2;
//   inequalityconstraintParams1[1].P = -M2;
//   inequalityconstraintParams1[0].b = -(M2*(desired_ddx_aug - dJ_Aug*mdqBody) + mhh - mTauLim);
//   inequalityconstraintParams1[1].b = -(-M2*(desired_ddx_aug - dJ_Aug*mdqBody) - mhh - mTauLim);
//   std::cout << "Check 9" <<  std::endl;
//   opt1.add_inequality_mconstraint(constraintFunc_RESCLF2, &inequalityconstraintParams1[0],
//                                  inequalityconstraintTol1);
//   opt1.add_inequality_mconstraint(constraintFunc_RESCLF2, &inequalityconstraintParams1[1],
//                                  inequalityconstraintTol1);
//   std::cout << "Check 10" <<  std::endl;
//   std::vector<double> ddqBodyRef_vec1(n1+1);
//   Eigen::VectorXd::Map(&ddqBodyRef_vec1[0], mddqBodyRef1.size()) = mddqBodyRef1;
//   std::cout << "size of mddqBodyRef  = " << mddqBodyRef1.rows() << "*" << mddqBodyRef1.cols() <<std::endl;
//   try {
//     std::cout << "Check 11" <<  std::endl;
//     opt1.set_xtol_rel(1e-4);
//     opt1.set_maxtime(0.1);
//     nlopt::result result = opt1.optimize(ddqBodyRef_vec1, minf);
//   } 
//   catch (std::exception& e) {
//     std::cout << "nlopt failed: " << e.what() << std::endl;
//   }
//   std::cout << "Check 12" <<  std::endl;
//   for(int i = 0;i < n1+1;i++) mddqBodyRef1(i,0) = ddqBodyRef_vec1[i]; 



//   Eigen::VectorXd bodyTorques1 = M2*(desired_ddx_aug + mddqBodyRef1.head(n1) - dJ_Aug*mdqBody) + mhh;

  // std::cout << "bodyTorques from RESCLF = \n" << bodyTorques1 << std::endl;
  // std::cout << "bodyTorques from Munzir = \n" << bodyTorques_Munzir << std::endl;

  // ************************************ Torques
  Eigen::VectorXd bodyTorques = bodyTorques1;
  mForces(0) = -mR / mL * tau_0 - bodyTorques(0) / 2;
  mForces(1) = mR / mL * tau_0 - bodyTorques(0) / 2;
  mForces.tail(mOptDim - 1) = bodyTorques.tail(mOptDim - 1);

  const std::vector<size_t> index{6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
                                  16, 17, 18, 19, 20, 21, 22, 23, 24};
  mRobot->setForces(index, mForces);

  if (mSteps < 0) {
    std::cout << "PEER: " << mPEER.rows() << " x " << mPEER.cols() << std::endl;
    std::cout << "PEEL: " << mPEEL.rows() << " x " << mPEEL.cols() << std::endl;
    std::cout << "PBal: " << mPBal.rows() << " x " << mPBal.cols() << std::endl;
    std::cout << "PPose: " << mPPose.rows() << " x " << mPPose.cols()
              << std::endl;
    std::cout << "PSpeedReg: " << mPSpeedReg.rows() << " x "
              << mPSpeedReg.cols() << std::endl;
    std::cout << "PReg: " << mPReg.rows() << " x " << mPReg.cols() << std::endl;
    // std::cout << "PxdotReg: " << PxdotReg.rows() << " x " << PxdotReg.cols()
    // << std::endl;
    std::cout << "bEER: " << mbEER.rows() << " x " << mbEER.cols() << std::endl;
    std::cout << "bEEL: " << mbEEL.rows() << " x " << mbEEL.cols() << std::endl;
    std::cout << "bBal: " << mbBal.rows() << " x " << mbBal.cols() << std::endl;
    std::cout << "bPose: " << mbPose.rows() << " x " << mbPose.cols()
              << std::endl;
    std::cout << "bSpeedReg: " << mbSpeedReg.rows() << " x "
              << mbSpeedReg.cols() << std::endl;
    std::cout << "bReg: " << mbReg.rows() << " x " << mbReg.cols() << std::endl;
    // std::cout << "bxdotReg: " << bxdotReg.rows() << " x " << bxdotReg.cols()
    // << std::endl;
  }

  if (mSteps < 0) {
    std::cout << "ddqBodyRef: " << std::endl;
    for (int i = 0; i < 18; i++) {
      std::cout << mddqBodyRef(i) << ", ";
    }
    std::cout << std::endl;
    std::cout << "ddqBodyRef_vec: " << std::endl;
    for (int i = 0; i < 18; i++) {
      std::cout << ddqBodyRef_vec[i] << ", ";
    }
    std::cout << std::endl;
  }

  // if(mSteps%(maxTimeSet==1?30:30) == 0)
  if (false) {
    std::cout << "mForces: " << mForces(0);
    for (int i = 1; i < 19; i++) {
      std::cout << ", " << mForces(i);
    }
    std::cout << std::endl;

    // Print the objective function components
    std::cout << "EEL loss: " << pow((mPEEL * mddqBodyRef - mbEEL).norm(), 2)
              << std::endl;
    std::cout << "EER loss: " << pow((mPEER * mddqBodyRef - mbEER).norm(), 2)
              << std::endl;
    std::cout << "OrL loss: " << pow((mPOrL * mddqBodyRef - mbOrL).norm(), 2)
              << std::endl;
    std::cout << "OrR loss: " << pow((mPOrR * mddqBodyRef - mbOrR).norm(), 2)
              << std::endl;
    std::cout << "Bal loss: " << pow((mPBal * mddqBodyRef - mbBal).norm(), 2)
              << std::endl;
    std::cout << "Pose loss: " << pow((mPPose * mddqBodyRef - mbPose).norm(), 2)
              << std::endl;
    std::cout << "Speed Reg loss: "
              << pow((mPSpeedReg * mddqBodyRef - mbSpeedReg).norm(), 2)
              << std::endl;
    std::cout << "Reg loss: " << pow((mPReg * mddqBodyRef - mbReg).norm(), 2)
              << std::endl;
  }
}

//==============================================================================
dart::dynamics::SkeletonPtr Controller::getRobot() const { return mRobot; }

//==============================================================================
dart::dynamics::BodyNode* Controller::getEndEffector(
    const std::string& s) const {
  if (!s.compare("left")) {
    return mLeftEndEffector;
  } else if (!s.compare("right")) {
    return mRightEndEffector;
  }
}

//==============================================================================
void Controller::keyboard(unsigned char /*_key*/, int /*_x*/, int /*_y*/) {}

package Avionics "Library of Avionics models"
  // GUE
  extends Modelica.Icons.Package;
  package Icons "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Icons;
  package Trash "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Trash;
  package Utilities "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Utilities;
  package Atmosphere "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Atmosphere;
  package Aerodynamics "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Aerodynamics;
  package Transformations "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Transformations;
  package Propulsion "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
  end Propulsion;
  package Examples "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    /* model GeometricTrackingControllerSE3
    extends Modelica.Icons.Example;
    Avionics.Examples.TrajectoriePosition1 trajectorieposition11 annotation(Placement(visible = true, transformation(origin = {-69.3578,-5.87156}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  equation
    connect(expsine1.y,se3dynamics_v11.f_real) annotation(Line(points = {{-61.6624,-8.44037},{-48.4404,-8.44037},{-48.4404,-8.70459},{-48.8807,-8.70459}}));
  end GeometricTrackingControllerSE3;
*/
    model test
      TrajectoriePosition1 Hello1 annotation(Placement(visible = true, transformation(origin = {-52.0262,9.22528}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    end test;
    model SI3MOs
      extends Modelica.Blocks.Icons.Block;
      //  Modelica.Blocks.Interfaces.RealInput u[3] annotation(Placement(visible = true, transformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // annotation(Placement(visible = true, transformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0)), Diagram(), Icon());
      Modelica.Blocks.Interfaces.RealOutput x annotation(Placement(visible = true, transformation(origin = {107.523,64.5872}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {107.523,64.5872}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput y annotation(Placement(visible = true, transformation(origin = {107.156,2.20183}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {107.156,2.20183}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput z annotation(Placement(visible = true, transformation(origin = {108.257,-60.5505}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {108.257,-60.5505}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealVectorInput u[3] annotation(Placement(visible = true, transformation(origin = {-98.7156,6.6055}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-98.7156,6.6055}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // parameter Integer nin = 3 "Number of inputs" annotation(Placement(visible = true, transformation(origin = {0,0.366972}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // Modelica.Blocks.Interfaces.RealInput u[nin] "Connector of Real input signals" annotation(Placement(visible = true, transformation(origin = {-120,0.733945}, extent = {{-20,-20},{20,20}}, rotation = 0)));
    equation
      x = u[1];
      y = u[2];
      z = u[3];
    end SI3MOs;
    model test2
      Examples.TrajectoriePosition1 trajectorieposition11 annotation(Placement(visible = true, transformation(origin = {-74.504,73.4453}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Examples.SI3MOs si3mos1 annotation(Placement(visible = true, transformation(origin = {-34.6715,73.5777}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Dynamics.SE3Dynamics_v1 se3dynamics_v11 annotation(Placement(visible = true, transformation(origin = {-32.2895,46.0522}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    end test2;
    model GeometricTrackingOnSE3
      extends Modelica.Icons.Example;
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-63.8532,38.8991}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.7523,1.10092}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.se3Controller se3controller1 annotation(Placement(visible = true, transformation(origin = {-10.7081,21.0708}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Bodies.SE3Dynamics_v1 se3dynamics_v11 annotation(Placement(visible = true, transformation(origin = {24.4477,20.9407}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      /*  connect(trajectoryposition11.signal,se3controller1.x_d);
  connect(trajectoryheading11.signal,se3controller1.b1_d);
  connect(se3controller1.Laws,se3dynamics_v11.Laws);
  connect(se3dynamics_v11.TV,se3controller1.TV);
*/
      connect(trajectoryposition11.signal,se3controller1.x_d);
      connect(trajectoryheading11.signal,se3controller1.b1_d);
      connect(se3controller1.Laws,se3dynamics_v11.Laws);
      connect(se3dynamics_v11.TV,se3controller1.TV);
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end GeometricTrackingOnSE3;
  end Examples;
  package Bodies "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    record GenericRotational3dDisplacement
      import SI = Modelica.SIunits;
      SI.Angle n[3,1] "Angle";
      SI.AngularVelocity Omega[3,1] "Angular Velocity";
      SI.AngularAcceleration dOmega[3,1];
      //equation
      //algorithm
      // v = der(x);
      // a = der(v);
    end GenericRotational3dDisplacement;
    record GenericLinear3dDisplacement
      import SI = Modelica.SIunits;
      SI.Position x[3,1] "Position";
      SI.Velocity v[3,1] "Velocity";
      SI.Velocity a[3,1] "Accelertion";
    end GenericLinear3dDisplacement;
    class partialModel
      // SE(3)
      // Ronan Cimadure, 2013-03-04
      import SI = Modelica.SIunits;
      //import constants;
      parameter String name = "GTC SE(3) model";
      replaceable parameter SI.Mass m;
      replaceable parameter SI.MomentOfInertia J[3,3] "Inertia matrix with respect of body-fixed frame";
      replaceable parameter SI.Acceleration g = Modelica.Constants.g_n;
      extends GenericLinear3dDisplacement;
      extends GenericRotational3dDisplacement;
      //  final constant SI.Acceleration g_n = 9.80665 "Standard acceleration of gravity on earth";
    end partialModel;
    record GenericDynamics
      // SE(3)
      // Ronan Cimadure, 2013-03-04
      import SI = Modelica.SIunits;
      extends GenericLinear3dDisplacement;
      extends GenericRotational3dDisplacement;
    end GenericDynamics;
    class se3Body
      import SI = Modelica.SIunits;
      replaceable parameter SI.Mass m;
      replaceable parameter SI.MomentOfInertia J[3,3] "Inertia matrix with respect of body-fixed frame";
      replaceable parameter SI.Acceleration g = Modelica.Constants.g_n;
    end se3Body;
    model SE3Dynamics_v1
      import SI = Modelica.SIunits;
      import Constants = Modelica.Constants;
      extends GenericDynamics;
      constant Real e3[3,1] = [0;0;1];
      //
      //
      parameter String name = "GTC SE(3) model";
      parameter SI.Position x_start[3,1] = [0;0;0] "Start Positions";
      parameter SI.Velocity v_start[3,1] = [0;0;0] "Start Velocityies";
      parameter SI.Angle n_start[3,1] = [0;0;0] "Start Angles";
      parameter SI.AngularVelocity Omega_start[3,1] = [0;0;0] "Start Angular Velocity";
      parameter SI.Mass m = 8.34;
      // 4.34
      parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      parameter SI.Acceleration g = Modelica.Constants.g_n;
      //
      input Avionics.Interfaces.se3CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      output Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      //  input Avionics.Interfaces.se3QuadrotorParamsConnector P "quadrotor Parameter" annotation(Placement(visible = true, transformation(origin = {-11.7431,96.1468}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-11.7431,96.1468}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    protected
      SI.Force f;
      SI.MomentOfForce M[3,1];
      SI.Angle R[3,3] "Start Angular Matrix";
    initial equation
      R = Avionics.Functions.Rxyz(n_start[1,1], n_start[2,1], n_start[3,1]);
      Omega = Omega_start;
      x = x_start;
      v = v_start;
    equation
      f = Laws.f;
      vector(M) = Laws.M;
      // LINEAR
      m * a = m * g * e3 - f * R * e3;
      a = der(v);
      v = der(x);
      // ANGULAR
      vector(J * dOmega) + cross(vector(Omega), vector(J * Omega)) = M[:,1];
      dOmega = der(Omega);
      R * skew(Omega[:,1]) = der(R);
      //der(R) = R * skew(Omega); // original
      // Add
      n = Avionics.Functions.vex(R);
      TV.x = x;
      TV.v = v;
      TV.R = R;
      TV.Omega = Omega;
      annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-34.4954,23.8532},{34.8623,-15.0459}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
      annotation(Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-38.5321,61.2844},{33.0275,-33.3945}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
    end SE3Dynamics_v1;
    model dynamics
      import SI = Modelica.SIunits;
      import Constants = Modelica.Constants;
      extends GenericDynamics;
      constant Real e3[3,1] = [0;0;1];
      //
      //
      parameter String name = "GTC SE(3) model";
      parameter SI.Position x_start[3,1] = [0;0;0] "Start Positions";
      parameter SI.Velocity v_start[3,1] = [0;0;0] "Start Velocityies";
      parameter SI.Angle n_start[3,1] = [0;0;0] "Start Angles";
      parameter SI.AngularVelocity Omega_start[3,1] = [0;0;0] "Start Angular Velocity";
      parameter SI.Mass m = 8.34;
      // 4.34
      parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      parameter SI.Acceleration g = Modelica.Constants.g_n;
      //
      /*  input Avionics.Interfaces.se3.commandLaws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
      input Avionics.Interfaces.se3.commandLawsFlow Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      output Avionics.Interfaces.se3.trackVariables TV annotation(Placement(visible = true, transformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    protected
      // SI.Force f;
      SI.MomentOfForce M[3,1];
      SI.Angle R[3,3] "Start Angular Matrix";
    initial equation
      R = Avionics.Functions.Rxyz(n_start[1,1], n_start[2,1], n_start[3,1]);
      Omega = Omega_start;
      x = x_start;
      v = v_start;
    equation
      // flow
      vector(M) = Laws.M;
      //
      zeros(3) = Laws.f + TV.f;
      zeros(3) = Laws.M + TV.M;
      //  TV.f = Laws.f;
      //  TV.M = Laws.M;
      // LINEAR
      m * a = m * g * e3 - Laws.f[3] * R * e3;
      a = der(v);
      v = der(x);
      // ANGULAR
      vector(J * dOmega) + cross(vector(Omega), vector(J * Omega)) = M[:,1];
      dOmega = der(Omega);
      R * skew(Omega[:,1]) = der(R);
      //der(R) = R * skew(Omega); // original
      // Add
      n = Avionics.Functions.vex(R);
      TV.x = x;
      TV.v = v;
      TV.R = R;
      TV.Omega = Omega;
      annotation(Diagram(), Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2}), graphics = {Text(origin = {2.16216,-8.64865}, lineColor = {0,0,255}, extent = {{-38.5321,61.2844},{33.0275,-33.3945}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
    end dynamics;
    model quadcopter
      import SI = Modelica.SIunits;
      import Constants = Modelica.Constants;
      extends GenericDynamics;
      constant Real e3[3,1] = [0;0;1];
      //
      //
      parameter String name = "GTC SE(3) model";
      parameter SI.Position x_start[3,1] = [0;0;0] "Start Positions";
      parameter SI.Velocity v_start[3,1] = [0;0;0] "Start Velocityies";
      parameter SI.Angle n_start[3,1] = [0;0;0] "Start Angles";
      parameter SI.AngularVelocity Omega_start[3,1] = [0;0;0] "Start Angular Velocity";
      parameter SI.Mass m = 8.34;
      // 4.34
      parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      parameter SI.Acceleration g = Modelica.Constants.g_n;
      //
      /*  input Avionics.Interfaces.se3.commandLaws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
      /*  input Avionics.Interfaces.se3.Flange_a_Laws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
      output Avionics.Interfaces.se3.trackVariables TV annotation(Placement(visible = true, transformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Modelica.Mechanics.Translational.Interfaces.Flange_a f_a annotation(Placement(visible = true, transformation(origin = {-95.3405,40.1434}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-95.3405,40.1434}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b M_a annotation(Placement(visible = true, transformation(origin = {-100.358,-45.1613}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-100.358,-45.1613}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    protected
      // SI.Force f;
      SI.MomentOfForce M[3,1];
      SI.Angle R[3,3] "Start Angular Matrix";
    initial equation
      R = Avionics.Functions.Rxyz(M_a.phi[1], M_a.phi[2], M_a.phi[3]);
      Omega = Omega_start;
      x = [f_a.s[1];f_a.s[2];f_a.s[3]];
      v = v_start;
    equation
      // flow
      vector(M) = M_a.tau;
      //
      zeros(3) = f_a.f + TV.f;
      zeros(3) = M_a.tau + TV.tau;
      //  TV.f = Laws.f;
      //  TV.M = Laws.M;
      // LINEAR
      m * a = m * g * e3 - f_a.f[3] * R * e3;
      a = der(v);
      v = der(x);
      // ANGULAR
      vector(J * dOmega) + cross(vector(Omega), vector(J * Omega)) = M[:,1];
      dOmega = der(Omega);
      R * skew(Omega[:,1]) = der(R);
      //der(R) = R * skew(Omega); // original
      // Add
      n = Avionics.Functions.vex(R);
      TV.x = x;
      TV.v = v;
      TV.R = R;
      TV.Omega = Omega;
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
      annotation(Diagram(), Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2}), graphics = {Text(origin = {2.16216,-8.64865}, lineColor = {0,0,255}, extent = {{-38.5321,61.2844},{33.0275,-33.3945}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
    end quadcopter;
  end Bodies;
  package Functions "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    function vex
      /*
%VEX Convert skew-symmetric matrix to vector
%
% V = VEX(S) is the vector (3x1) which has the skew-symmetric matrix S (3x3)
%
%           | 0   -vz  vy|
%           | vz   0  -vx|
%           |-vy   vx  0 |
%
% Notes::
% - This is the inverse of the function SKEW().
% - No checking is done to ensure that the matrix is actually skew-symmetric.
% - The function takes the mean of the two elements that correspond to each unique
%   element of the matrix, ie. vx = 0.5*(S(3,2)-S(2,3))
%
% See also SKEW.

  RONAN 2013-03-11
% Copyright (C) 1993-2011, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for Matlab (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
*/
      input Real S[3,3];
      output Real v[3,1];
    algorithm
      //if isrot(S) || ishomog(S)
      v:=0.5 * [S[3,2] - S[2,3];S[1,3] - S[3,1];S[2,1] - S[1,2]];
      //else
      //    error('argument must be a 3x3 matrix');
      //end
    end vex;
    function Rxyz
      /*
%% Rotation Matrix (3x3)
% X-Y-Z *Euler* angle
%        ----- 
%
%
% A to B R_{XYZ} (lpha,eta,\gamma) =  R_X(lpha) R_Y(eta) R_Z(\gamma)
%
% sources :  p. 375 in Introduction to Robotics 3th Edition, by Jonh J. Greg 
%            http://en.wikipedia.org/wiki/Euler_angles : X1,Y2,Z3
%
% [ r_11, r_12, r_13
%   r_21, r_22, r_23
%   r_31, r_32, r_33 ]
%
% Ronan CIMADURE
% 2012-11-23
%
% -----------------------
% 2012-11-25 : correct Error
% */
      import SI = Modelica.SIunits;
      import SIm = Modelica.Math;
      input SI.Angle a;
      input SI.Angle b;
      input SI.Angle c;
      output SI.Angle R[3,3];
    protected
      SI.Angle ca,cb,cc;
      SI.Angle sa,sb,sc;
      //annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 0.000001));
    algorithm
      ca:=SIm.cos(a);
      sa:=SIm.sin(a);
      cb:=SIm.cos(b);
      sb:=SIm.sin(b);
      cc:=SIm.cos(c);
      sc:=SIm.sin(c);
      R:=[cb * cc,-cb * sc,sb;sa * sb * cc + ca * sc,-sa * sb * sc + ca * cb,-sa * cb;-ca * sb * cc + sa * sc,ca * sb * sc + sa * cc,ca * cb];
    end Rxyz;
    function se3R_d
      import SI = Modelica.SIunits;
      input SI.Angle b1_d[3] "";
      input SI.Angle b3_d[3] "";
      output SI.Angle R[3,3] "Angular Matrix";
    protected
      SI.Angle b31_d[3] "";
      SI.Angle b23_d[3] "";
      SI.Angle b2_d[3] "";
    algorithm
      // Angular
      // R_d : Attitude Tracking
      b31_d:=cross(b3_d, b1_d);
      b2_d:=b31_d / Modelica.Math.Vectors.norm(b31_d);
      b23_d:=cross(b2_d, b3_d);
      //  R := [matrix(b23_d,3,1), matrix(b2_d,3,1), matrix(b3_d,3,1)];
      R:=[b23_d[1],b2_d[1],b3_d[1];b23_d[2],b2_d[2],b3_d[2];b23_d[3],b2_d[3],b3_d[3]];
    end se3R_d;
    function se3Error "e = (TV,W)  [12, e.R==3]"
      import SI = Modelica.SIunits;
      input Avionics.Types.se3TrackingVariable TV;
      input Avionics.Types.se3TrackingVariable W;
      output Avionics.Types.se3Error e "Error Variable";
    algorithm
      e.x:=TV.x - W.x;
      e.v:=TV.v - W.v;
      e.R:=Avionics.Functions.vex(1 / 2 * (transpose(W.R) * TV.R - transpose(TV.R * W.R)));
      e.Omega:=TV.Omega - transpose(TV.R) * W.R * W.Omega;
    end se3Error;
    model R_derivative "Matrix Derivative (for non-linear system)"
      extends Modelica.Icons.Function;
      input Real R[3,3] "Input matrix";
      output Real dotR[3,3];
      //  input Real R[:,:] "Input matrix";
      //  output Real dotR[:,:]
    protected
      Modelica.Blocks.Continuous.Derivative deriv;
      //algorithm
    equation
      for i1 in 1:size(R, 1) loop
      for i2 in 1:size(R, 2) loop
      deriv.u = R[i1,i2];
      dotR[i1,i2] = deriv.y;

      end for;

      end for;
      // deriv.u := R[i1,i2];
      // dotR[i1, i2] := deriv.y;
    end R_derivative;
    block Derivative "Approximated derivative block"
      import Modelica.Blocks.Types.Init;
      parameter Real k(unit = "1") = 1 "Gains";
      parameter Modelica.SIunits.Time T(min = Modelica.Constants.small) = 0.01 "Time constants (T>0 required; T=0 is ideal derivative block)";
      parameter Modelica.Blocks.Types.Init initType = Modelica.Blocks.Types.Init.NoInit "Type of initialization (1: no init, 2: steady state, 3: initial state, 4: initial output)" annotation(Evaluate = true, Dialog(group = "Initialization"));
      parameter Real x_start[3,3] = zeros(3, 3) "Initial or guess value of state" annotation(Dialog(group = "Initialization"));
      parameter Real y_start[3,3] = zeros(3, 3) "Initial value of output (= state)" annotation(Dialog(enable = initType == Init.InitialOutput, group = "Initialization"));
      //  extends Interfaces.SISO;
      //  extends Interfaces.MIMO(nin=9, nout=9);
      input Real u[3,3];
      output Real y[3,3];
      output Real x[3,3](start = x_start) "State of block";
    protected
      parameter Boolean zeroGain = abs(k) < Modelica.Constants.eps;
    initial equation
      if initType == Init.SteadyState then
        der(x) = zeros(3, 3);
      elseif initType == Init.InitialState then
        x = x_start;
      elseif initType == Init.InitialOutput then
        if zeroGain then
          x = u;
        else
          y = y_start;
        end if;
      else
      end if;
      // PROB HERE
    equation
      der(x) = if zeroGain then zeros(3, 3) else (u - x) / T;
      y = if zeroGain then zeros(3, 3) else k / T * (u - x);
      annotation(Documentation(info = "<html>
<p>
This blocks defines the transfer function between the
input u and the output y
(element-wise) as <i>approximated derivative</i>:
</p>
<pre>
             k * s
     y = ------------ * u
            T * s + 1
</pre>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general block <b>TransferFunction</b> instead
and model a derivative block with parameters<br>
b = {k,0}, a = {T, 1}.
</p>

<p>
If k=0, the block reduces to y=0.
</p>
</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {2,2}), graphics = {Line(points = {{-80,78},{-80,-90}}, color = {192,192,192}),Polygon(points = {{-80,90},{-88,68},{-72,68},{-80,90}}, lineColor = {192,192,192}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid),Line(points = {{-90,-80},{82,-80}}, color = {192,192,192}),Polygon(points = {{90,-80},{68,-72},{68,-88},{90,-80}}, lineColor = {192,192,192}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid),Line(points = {{-80,-80},{-80,60},{-70,17.95},{-60,-11.46},{-50,-32.05},{-40,-46.45},{-30,-56.53},{-20,-63.58},{-10,-68.51},{0,-71.96},{10,-74.37},{20,-76.06},{30,-77.25},{40,-78.07},{50,-78.65},{60,-79.06}}, color = {0,0,127}),Text(extent = {{-30,14},{86,60}}, lineColor = {192,192,192}, textString = "DT1"),Text(extent = {{-150,-150},{150,-110}}, lineColor = {0,0,0}, textString = "k=%k")}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {2,2}), graphics = {Text(extent = {{-54,52},{50,10}}, lineColor = {0,0,0}, textString = "k s"),Text(extent = {{-54,-6},{52,-52}}, lineColor = {0,0,0}, textString = "T s + 1"),Line(points = {{-50,0},{50,0}}, color = {0,0,0}),Rectangle(extent = {{-60,60},{60,-60}}, lineColor = {0,0,255}),Line(points = {{-100,0},{-60,0}}, color = {0,0,255}),Line(points = {{60,0},{100,0}}, color = {0,0,255})}));
    end Derivative;
    function dotR
      input Real k;
      input Real T;
      input Real u[3,3];
      output Real y[3,3];
    protected
      Real x[3,3];
      Real dx[3,3];
    algorithm
      y:=k / T * (u - x);
      //y := zeros(3,3);
      dx:=der(x);
      dx:=(u - x) / T;
      //equation
      //  der(x) = (u - x) / T;
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end dotR;
  end Functions;
  package Test "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    function M3MultV3_v1
      Real R[3,3] = [1,2,3;4,5,6;7,8,9];
      Real V[3,1] = [1;3;7];
      Real res[3,1];
    algorithm
      res:=R * V;
    end M3MultV3_v1;
    function M3MultV3_v2
      Real R[3,3] = [1,2,3;4,5,6;7,8,9];
      Real V[3] = [1,3,7];
      Real res[3,1];
    algorithm
      res:=R * V;
    end M3MultV3_v2;
    function M3MultV3_v3
      Real R[3,3] = [1,2,3;4,5,6;7,8,9];
      Real V[3] = [1,3,7];
      Real res[3];
    algorithm
      res:=R * V;
    end M3MultV3_v3;
    model dynamics_v1
      Avionics.Sources.CommandLaws commandlaws1(f = 40.0) annotation(Placement(visible = true, transformation(origin = {-16.1605,52.5697}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Bodies.SE3Dynamics_v1 se3dynamics_v11 annotation(Placement(visible = true, transformation(origin = {12.1097,51.9142}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      connect(commandlaws1.Laws,se3dynamics_v11.Laws);
      annotation(experiment(StartTime = 0.0, StopTime = 50.0, Tolerance = 0.000001), Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-34.4954,23.8532},{34.8623,-15.0459}}, textString = "Dynamics", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
    end dynamics_v1;
    model se3controller2
      extends Modelica.Icons.Example;
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-63.8532,38.8991}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.7523,1.10092}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.se3Controller se3controller1 annotation(Placement(visible = true, transformation(origin = {-10.7081,21.0708}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Types.se3TrackingVariable T annotation(Placement(visible = true, transformation(origin = {-9.84529,-11.5331}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      connect(trajectoryposition11.signal,se3controller1.x_d);
      connect(trajectoryheading11.signal,se3controller1.b1_d);
      connect(se3controller1.Laws,se3dynamics_v11.Laws);
      connect(se3dynamics_v11.TV,se3controller1.TV);
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end se3controller2;
    model se3controller
      extends Modelica.Icons.Example;
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-63.8532,38.8991}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.7523,1.10092}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.se3Controller se3controller1 annotation(Placement(visible = true, transformation(origin = {-10.7081,21.0708}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Types.se3TrackingVariable T annotation(Placement(visible = true, transformation(origin = {-9.84529,-11.5331}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      //  connect(trajectoryposition11.signal,se3controller1.x_d);
      //  connect(trajectoryheading11.signal,se3controller1.b1_d);
      //   connect(T,se3controller1.TV);
      se3controller1.x_d = trajectoryposition11.signal;
      se3controller1.b1_d = trajectoryheading11.signal;
      T.x = [0;0;0];
      T.v = [0;0;0];
      T.R = diagonal({1,1,1});
      T.Omega = [0;0;0];
      se3controller1.TV = T;
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end se3controller;
    model trajectoryTracking
      Avionics.Controller.SE3.se3TrajectoryTracking se3trajectorytracking1 annotation(Placement(visible = true, transformation(origin = {-23.1434,28.3247}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-73.2297,30.3972}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Interfaces.se3TrackConnector T annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    equation
      se3trajectorytracking1.x_d = trajectoryposition11.signal;
      T.x = [0;0;0];
      T.v = [0;0;0];
      T.R = diagonal({1,1,1});
      T.Omega = [0;0;0];
      se3trajectorytracking1.TV = T;
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end trajectoryTracking;
    model transposeMatrix
      Real R[3,3];
      output Real out[3,3];
    equation
      R = [[1,2,3];[4,5,6];[7,8,9]];
      out = transpose(R);
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end transposeMatrix;
    model attitudeTracking
      import SI = Modelica.SIunits;
      SI.Angle b3_d[3];
      Avionics.Controller.SE3.se3AttitudeTracking se3attitudetracking1 annotation(Placement(visible = true, transformation(origin = {-0.0000303352,-26.7229}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-39.0999,0.843849}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Interfaces.se3TrackConnector T annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    equation
      se3attitudetracking1.b1_d = trajectoryheading11.signal;
      b3_d = {0,0,0};
      T.x = [0;0;0];
      T.v = [0;0;0];
      T.R = diagonal({1,1,1});
      T.Omega = [0;0;0];
      se3attitudetracking1.TV = T;
      se3attitudetracking1.b3_d = b3_d;
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end attitudeTracking;
    model derR_d_not_working
      Avionics.Controller.SE3.der_R_d der_r_d1 annotation(Placement(visible = true, transformation(origin = {-45.8509,15.7525}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      der_r_d1.b1_d = {0,0,0};
      der_r_d1.b3_d = {0,0,0};
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end derR_d_not_working;
    model attitudeTracking2
      import SI = Modelica.SIunits;
      import SIm = Modelica.Math;
      SI.Angle b3_d annotation(Placement(visible = true, transformation(origin = {41.3502,3.09423}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-66.1041,75.668}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.SE3.se3TrajectoryTracking se3trajectorytracking1 annotation(Placement(visible = true, transformation(origin = {-28.4107,75.3868}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.SE3.se3AttitudeTracking se3attitudetracking1 annotation(Placement(visible = true, transformation(origin = {1.96903,41.9128}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Interfaces.se3TrackConnector T annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    equation
      /*se3attitudetracking1.b1_d = trajectoryheading11.signal;
  T.x = [0;0;0];
  T.v = [0;0;0];
  T.R = diagonal({1,1,1});
  T.Omega = [0;0;0];
  se3attitudetracking1.TV = T;
  connect(se3trajectorytracking1.b3_d,se3attitudetracking1.b3_d);
*/
    end attitudeTracking2;
    model derivative1_not_working
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.5216,49.0501}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // Modelica.Blocks.Continuous.Derivative derivative1 annotation(Placement(visible = true, transformation(origin = {-17.2712,49.3955}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Functions.R_derivative derivative1;
      Real R[3,3];
      Real dotR[3,3];
    equation
      R = Avionics.Functions.Rxyz(trajectoryheading11.signal[1], trajectoryheading11.signal[2], trajectoryheading11.signal[3]);
      connect(R,derivative1.R);
      connect(derivative1.dotR,dotR);
      //dotR = Avionics.Functions.R_derivative(R);
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end derivative1_not_working;
    model derivation_of_R
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-61.4853,37.3057}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Real R[3,3];
      Real R2[3,3];
      // Real derR;
      //Real derR1;
      //Real a;
      input Real tolerance = 100 * Modelica.Constants.eps "Relative tolerance of solution u";
      Avionics.Functions.Derivative derivative1 annotation(Placement(visible = true, transformation(origin = {52.8497,31.0881}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    initial equation
      R2 = zeros(3, 3);
    equation
      R = Avionics.Functions.Rxyz(trajectoryheading11.signal[1], trajectoryheading11.signal[2], trajectoryheading11.signal[3]);
      //  derR = der(R[2,2]);
      //derR = der(trajectoryheading11.signal[1]);
      // derR1 = der(R[1,3]);
      //R2 = R;
      //a = R2[1,1];
      //derR1 = der(a);
      //R2 = der(R);
      // R2 = Modelica.Math.Nonlinear.solveOneNonlinearEquation(function der(R), -0.5, 10, tolerance);
      connect(R,derivative1.u);
      connect(derivative1.y,R2);
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end derivation_of_R;
    model dotR
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-61.4853,37.3057}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Real R[3,3];
      Real dotR[3,3];
    initial equation
      dotR = zeros(3, 3);
    equation
      R = Avionics.Functions.Rxyz(trajectoryheading11.signal[1], trajectoryheading11.signal[2], trajectoryheading11.signal[3]);
      dotR = Avionics.Functions.dotR(1, 0.01, R);
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
    end dotR;
    model dynamics
      Avionics.Bodies.dynamics dynamics1 annotation(Placement(visible = true, transformation(origin = {-21.6216,50.2703}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Avionics.Sources.CommandLaws commandlaws1(f = {0,0,82.5}) annotation(Placement(visible = true, transformation(origin = {-60,51.8919}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    equation
      connect(commandlaws1.Laws,dynamics1.Laws);
    end dynamics;
    model forces
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-72.5948,41.1079}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Avionics.Sources.laws laws1 annotation(Placement(visible = true, transformation(origin = {19.8251,15.1603}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-73.4694,-20.8322}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    equation
      connect(trajectoryposition11.signal,laws1.M_a) annotation(Line(points = {{-62.4694,-20.8322},{-39.6236,7.7392},{9.56267,11.1662},{9.56267,12.3324}}));
      connect(trajectoryheading11.signal,laws1.f_a) annotation(Line(points = {{-61.5948,41.1079},{-41.1079,41.1079},{-41.1079,16.9096},{-15.7143,18.1924},{9.35857,18.1924}}));
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end forces;
    model force
      Avionics.Sources.chose_f teee1 annotation(Placement(visible = true, transformation(origin = {-17.0686,11.1317}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Blocks.Sources.Sine sine1 annotation(Placement(visible = true, transformation(origin = {-66.0482,10.7607}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass1(m = 1) annotation(Placement(visible = true, transformation(origin = {24.8649,11.7432}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    equation
      connect(teee1.flange_b,mass1.flange_a) annotation(Line(points = {{-7.0394,11.0151},{14.6069,11.0151},{14.8649,10.8108},{14.8649,11.7432}}));
      connect(sine1.y,teee1.u) annotation(Line(points = {{-55.0482,10.7607},{-27.4583,10.3896},{-27.4583,10.986},{-27.1269,10.986}}));
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end force;
    model commandDynamics
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-70.3432,31.8039}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Avionics.Sources.laws laws1 annotation(Placement(visible = true, transformation(origin = {-29.2542,28.5993}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-70.5158,2.61978}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Avionics.Bodies.quadcopter quadcopter1 annotation(Placement(visible = true, transformation(origin = {11.3524,28.5993}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    equation
      connect(trajectoryheading11.signal,laws1.M_a) annotation(Line(points = {{-59.5158,2.61978},{-40.3883,2.61978},{-40.3883,25.7612},{-40.3883,25.7612}}));
      connect(trajectoryposition11.signal,laws1.f_a) annotation(Line(points = {{-59.3432,31.8039},{-39.9517,31.8039},{-39.9517,31.2191},{-39.9517,31.2191}}));
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end commandDynamics;
  end Test;
  package Controller "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    /*  record se3TrackingVariable
    import SI = Modelica.SIunits;
    SI.Position x[3,1] "Position";
    SI.Velocity v[3,1] "Velocity";
    replaceable SI.Angle R[3,3] "Angular Matrix";
    SI.AngularVelocity Omega[3,1] "Angular Velocity";
  end se3TrackingVariable;
  
*/
    class se3parameters
      import SI = Modelica.SIunits;
      //replaceable
      /*	parameter SI.Mass m = 1;
  replaceable Real x = 16 * m;
  replaceable Real v = 5.6 * m;
  replaceable Real R = 8.81;
  replaceable Real Omega = 2.54;
*/
      replaceable parameter SI.Mass m = 1;
      replaceable parameter Real x = 16 * m;
      replaceable parameter Real v = 5.6 * m;
      replaceable parameter Real R = 8.81;
      replaceable parameter Real Omega = 2.54;
      //equation
      //  x = 16 * m;
      //  v = 5.6 * m;
    end se3parameters;
    package SE3
      extends Modelica.Icons.InterfacesPackage;
      //
      //
      model se3Command
        /*
%% SE3THRUST
% 2012 -10-16
% Ronan CIMADURE

%% Tracking Controller Outputs
%the total thrust _f_ and the total moment _M_ can be written as
%U = []; % initilisation
% f  .eq 15

% a_d == W.a , wish acceleration
% e = se3Error

% e_3 : IN LINE [0 0 1] -> R*e_3
function f = se3Thrust(k,e,a_d,R,m,g,e_3)
  Re3 = (R*transpose(e_3)) ;

  coeff = (-k.x .* e.x - k.v .* e.v - m*g*e_3+m*a_d) ;
   
  f = - coeff * Re3;

end
*/
        import SI = Modelica.SIunits;
        se3parameters k;
        /*
  replaceable parameter SI.Mass m = 1;
  parameter Real k_x = 16 * m;
  parameter Real k_v = 5.6 * m;
  parameter Real k_r = 8.81;
  parameter Real k_Omega = 2.54;
*/
        parameter SI.Acceleration g = Modelica.Constants.g_n;
        //
        constant Real e3[3,1] = [0;0;1];
        //
        //  input GenericDynamics;
        input Avionics.SE3Dynamics_v1 B;
        input GenericLinear3dDisplacement e_lin;
        output Avionics.OutCommandLaws laws;
      equation
        //laws.f = (-k.x .* e_lin.x - k.v .* e_lin.v - B.m * B.g * e3 + B.m * B.a) * B.R * e3;
        // SE3Dynamics_v1
      end se3Command;
      model se3TrajectoryTracking
        import SI = Modelica.SIunits;
        se3parameters k "Controller parameters";
        parameter SI.Mass m = 8.34;
        parameter SI.Acceleration g = Modelica.Constants.g_n;
        input Modelica.Blocks.Interfaces.RealInput x_d[3] "Position" annotation(Placement(visible = true, transformation(origin = {-96.9179,14.8435}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-96.9179,14.8435}, extent = {{-12,-12},{12,12}}, rotation = 0)));
        output SI.Angle b3_d[3] "";
        output SI.Force f;
        input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
      protected
        // Vars
        Avionics.Types.se3TrackingVariableLinear e "Error Variable";
        Real track_coeff[3,1];
        SI.Velocity v_d[3,1] "desired Velocity";
        SI.Acceleration a[3,1] "Acceleration";
        constant Real e3[3,1] = [0.0;0.0;1.0];
        Real temp[3,1];
      equation
        vector(v_d) = der(x_d);
        a = der(TV.v);
        vector(e.x) = vector(TV.x) - x_d;
        //e.x = TV.x - transpose(x_d);
        e.v = TV.v - v_d;
        // Tracking Controller
        track_coeff = -k.x * e.x - k.v * e.v - m * g * e3 + m * a;
        b3_d = -vector(track_coeff) / Modelica.Math.Matrices.frobeniusNorm(track_coeff);
        // b3_d = -vector(track_coeff) / Modelica.Math.Vectors.norm(vector(track_coeff));
        // Control output
        // f :
        //  f = scalar(-transpose(track_coeff) * TV.R * e3);
        temp = TV.R * e3;
        f = -sum(transpose(track_coeff) * temp);
        //  b3_d = {1.1,1.2,0.3};
        annotation(Icon(), Diagram());
      end se3TrajectoryTracking;
      model der_R_d
        import SI = Modelica.SIunits;
        Avionics.Controller.SE3.R_d r_d1;
        input SI.Angle b1_d[3];
        input SI.Angle b3_d[3];
        output Real der_R_d[3,3];
      equation
        connect(b1_d,r_d1.b1_d);
        connect(b3_d,r_d1.b3_d);
        der_R_d = der(r_d1.R);
      end der_R_d;
      model R_d
        import SI = Modelica.SIunits;
        input SI.Angle b1_d[3];
        input SI.Angle b3_d[3];
        output Real R[3,3] "Angular Matrix";
        // output Real der_R[3,3];
        // Angular
        // R_d : Attitude Tracking
        //R = Avionics.Functions.se3R_d(b1_d, b3_d);
      protected
        SI.Angle b31_d[3] "";
        SI.Angle b23_d[3] "";
        SI.Angle b2_d[3] "";
      equation
        // Angular
        // R_d : Attitude Tracking
        b31_d = cross(b3_d, b1_d);
        //b2_d = b31_d / Modelica.Math.Matrices.frobeniusNorm(b31_d);
        b2_d = b31_d / Modelica.Math.Vectors.norm(b31_d);
        b23_d = cross(b2_d, b3_d);
        //  R := [matrix(b23_d,3,1), matrix(b2_d,3,1), matrix(b3_d,3,1)];
        R = [b23_d[1],b2_d[1],b3_d[1];b23_d[2],b2_d[2],b3_d[2];b23_d[3],b2_d[3],b3_d[3]];
      end R_d;
      model se3AttitudeTracking
        import SI = Modelica.SIunits;
        input SI.Angle b1_d[3];
        input SI.Angle b3_d[3];
        output SI.MomentOfForce M[3];
        parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
        se3parameters k "Controller parameters";
        input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
      protected
        Avionics.Functions.Derivative dotR;
        // Avionics.Controller.SE3.R_d r_d1 ;
        //
        Real temp_M[3,1];
        Real temp_M_cross[3,1];
        SI.Angle e_R[3,1];
        SI.AngularVelocity e_Omega[3,1] "Angular Velocity";
        Real R_d[3,3] "Angular Matrix";
        Real R_dd[3,3];
        Real der_R_d[3,3];
        Real inv_dotR_d[3,3];
        //   SI.Angle
        SI.AngularVelocity Omega_d[3,1] "Angular Velocity";
        //Real der_Omega_d[3,1];
        Real transp_R_d[3,3];
      equation
        // Angular
        // R_d : Attitude Tracking
        /*  connect(b1_d,r_d1.b1_d);
  connect(b3_d,r_d1.b3_d);
  R_d = r_d1.R;
*/
        R_d = Avionics.Functions.se3R_d(b1_d, b3_d);
        R_dd = R_d;
        dotR.u = R_dd;
        der_R_d = dotR.y;
        //  dotR.u = R_d;
        //  der_R_d = dotR.y;
        //
        //der_R_d = der(R_d);
        //der(R_d) = der_R_d;
        // R_d3 = Modelica.Math.Matrices.inv(R_d2);
        //
        //
        transp_R_d = transpose(R_d);
        inv_dotR_d = Modelica.Math.Matrices.inv(der_R_d);
        // W.Omega
        Omega_d = Avionics.Functions.vex(inv_dotR_d * R_d);
        //Omega_d = Avionics.Functions.vex(Modelica.Math.Matrices.inv(der(R_d)) * R_d);
        e_R = Avionics.Functions.vex(1 / 2 * transp_R_d * TV.R - transpose(TV.R * R_d));
        e_Omega = TV.Omega - transpose(TV.R) * R_d * Omega_d;
        // Control output
        // M :
        vector(temp_M_cross) = cross(vector(TV.Omega), vector(J * TV.Omega));
        temp_M = -k.R * e_R - k.Omega * e_Omega + temp_M_cross - J * (skew(vector(TV.Omega)) * transpose(TV.R) * R_d * Omega_d - transpose(TV.R) * R_d * der(Omega_d));
        //der_Omega_d = der(Omega_d);
        //  temp_M = -k.R * e_R - k.Omega * e_Omega + temp_M_cross;
        //transp_R = transpose(TV.R);
        //transp_R = diagonal({1.1,1.2,0.3});
        //temp_M = -k.R * e_R - k.Omega * e_Omega + temp_M_cross - J * (skew(vector(TV.Omega)) * transp_R * R_d * Omega_d - transp_R * R_d * der_Omega_d);
        M = vector(temp_M);
        // M = {1.1,1.2,0.3};
      end se3AttitudeTracking;
    end SE3;
    model se3QuadrotorParameters
      //  parameter Avionics.Types.se3QuadrotorParameters;
      import SI = Modelica.SIunits;
      parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      parameter SI.Mass m = 8.34 "Mass";
      //
      parameter SI.Distance d = 0.315 "arm lenght";
      //
      parameter SI.Distance c_t_f = 0.0008004 "Coefficient of Thrust Force";
      //;
      //
      parameter SI.Acceleration g = Modelica.Constants.g_n "gravity";
      //  import SI = Modelica.SIunits;
      output Avionics.Interfaces.se3QuadrotorParamsConnector P annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      P.J = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      P.m = 8.34;
      P.d = 0.315;
      P.c_t_f = 0.0008004;
      P.g = Modelica.Constants.g_n;
      annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-59.4495,-49.1743},{62.3853,48.4404}}, textString = "Parameters", textStyle = {TextStyle.Bold})}));
    end se3QuadrotorParameters;
    model ctrlParameters
      import SI = Modelica.SIunits;
      //replaceable
      /*	parameter SI.Mass m = 1;
  replaceable Real x = 16 * m;
  replaceable Real v = 5.6 * m;
  replaceable Real R = 8.81;
  replaceable Real Omega = 2.54;
*/
      replaceable parameter SI.Mass m = 1;
      replaceable Real x;
      // = 16 * m;
      replaceable Real v;
      // = 5.6 * m;
      replaceable Real R;
      replaceable Real Omega;
      parameter Real R_param = 8.81;
      parameter Real Omega_param = 2.54;
    equation
      x = 16 * m;
      v = 5.6 * m;
      R = R_param;
      Omega = Omega_param;
    end ctrlParameters;
    model se3Controller
      parameter String name = "GTC SE(3) controller model";
      import SI = Modelica.SIunits;
      import Interfaces = Avionics.Interfaces;
      //
      // In
      input Modelica.Blocks.Interfaces.RealInput b1_d[3] "Heading" annotation(Placement(visible = true, transformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      input Modelica.Blocks.Interfaces.RealInput x_d[3] "Position" annotation(Placement(visible = true, transformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // Out
      output Avionics.Interfaces.se3CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {100.45,1.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.45,1.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      //
      parameter SI.Mass m = 8.34;
      // 4.34
      parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      parameter SI.Acceleration g = Modelica.Constants.g_n;
      input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    protected
      Avionics.Controller.SE3.se3TrajectoryTracking se3trajectorytracking1 annotation(Placement(visible = true, transformation(origin = {-57.1027,43.3193}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.SE3.se3AttitudeTracking se3attitudetracking1 annotation(Placement(visible = true, transformation(origin = {-12.9395,31.7862}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      //  connect(x_d,se3trajectorytracking1.x_d);
      //  connect(b1_d,se3attitudetracking1.b1_d);
      se3trajectorytracking1.x_d = x_d;
      se3attitudetracking1.b1_d = b1_d;
      connect(se3trajectorytracking1.b3_d,se3attitudetracking1.b3_d);
      connect(TV,se3trajectorytracking1.TV);
      connect(TV,se3attitudetracking1.TV);
      Laws.f = se3trajectorytracking1.f;
      Laws.M = se3attitudetracking1.M;
      //Laws.f = 3.4567;
      //Laws.M = {1.1,1.2,0.3};
      annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-52.844,39.2661},{-82.5688,54.6789}}, textString = "Position", textStyle = {TextStyle.Bold}),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-52.1101,-45.1376},{-78.5321,-38.5321}}, textString = "Heading"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{48.8073,7.33945},{89.5413,-6.6055}}, textString = "Control Laws"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-19.4495,-85.8715},{28.6238,-79.266}}, textString = "Tracking Variables"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-29.7248,83.3028},{8.80734,60.1835}}, textString = "Parameters"),Text(rotation = -180, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{62.3853,2.93578},{-63.1193,-39.633}}, textString = "CONTROLLER", textStyle = {TextStyle.Bold})}));
    end se3Controller;
  end Controller;
  package Sources "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    block TrajectoryHeading1
      /*  // import Interfaces = Modelica.Blocks.Interfaces;
  import SIunits = Modelica.SIunits;
  extends Avionics.Interfaces.MO(nout = 3);
  //
  //  parameter Real amplitude = 1 "Amplitude of sine wave";
  //  parameter SIunits.Frequency freqHz(start = 2) "Frequency of sine wave";
  //  parameter SIunits.Angle phase = 0 "Phase of sine wave";
  //  parameter SIunits.Damping damping(start = 1) "Damping coefficient of sine wave";
  parameter Real offset = 0 "Offset of output signal";
  parameter SIunits.Time startTime = 0 "Output = offset for time < startTime";
protected
  constant Real pi = Modelica.Constants.pi annotation(Placement(visible = true, transformation(origin = {0.366972,0}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  /*  y[1] = offset + (if time < startTime then 0 else Modelica.Math.cos(pi * (time - startTime)));
  y[2] = offset + (if time < startTime then 0 else Modelica.Math.sin(pi * (time - startTime)));
  y[3] = offset + 0;

  signal[1] = offset + (if time < startTime then 0 else Modelica.Math.cos(pi * (time - startTime)));
  signal[2] = offset + (if time < startTime then 0 else Modelica.Math.sin(pi * (time - startTime)));
  signal[3] = offset + 0;

*/
      import SIunits = Modelica.SIunits;
      extends Avionics.Interfaces.MO(nout = 3);
      parameter Real amplitude[3] = {0.4,0.4,0.6} "Amplitude of movement";
      //  parameter SIunits.Frequency freqHz(start = 2) "Frequency of sine wave";
      parameter SIunits.Angle phase = 0 "Phase of movement";
      parameter SIunits.Time startTime = 0 "Output = offset for time < startTime" annotation(Placement(visible = true, transformation(origin = {1.10092,-31.1927}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      parameter Real offset = 0 "Offset of output signal" annotation(Placement(visible = true, transformation(origin = {2.93578,-63.1193}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 0.000001));
    protected
      constant Real pi = Modelica.Constants.pi annotation(Placement(visible = true, transformation(origin = {-0.733945,31.9266}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      /*  Real x;
  Real y;
  Real z;*/
    equation
      signal[1] = offset + (if time < startTime then 0 else amplitude[1] * time);
      signal[2] = offset + (if time < startTime then 0 else amplitude[2] * Modelica.Math.sin(pi * (time - startTime) + phase));
      signal[3] = offset + (if time < startTime then 0 else amplitude[3] * Modelica.Math.cos(pi * (time - startTime) + phase));
      //  signal = {x,y,z};
    end TrajectoryHeading1;
    block TrajectoryPosition1
      import SIunits = Modelica.SIunits;
      extends Avionics.Interfaces.MO(nout = 3);
      parameter Real amplitude[3] = {0.4,0.4,0.6} "Amplitude of movement";
      //  parameter SIunits.Frequency freqHz(start = 2) "Frequency of sine wave";
      parameter SIunits.Angle phase = 0 "Phase of movement";
      parameter SIunits.Time startTime = 0 "Output = offset for time < startTime" annotation(Placement(visible = true, transformation(origin = {1.10092,-31.1927}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      parameter Real offset = 0 "Offset of output signal" annotation(Placement(visible = true, transformation(origin = {2.93578,-63.1193}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 0.000001));
    protected
      constant Real pi = Modelica.Constants.pi annotation(Placement(visible = true, transformation(origin = {-0.733945,31.9266}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      /*  Real x;
  Real y;
  Real z;*/
    equation
      signal[1] = offset + (if time < startTime then 0 else amplitude[1] * time);
      signal[2] = offset + (if time < startTime then 0 else amplitude[2] * Modelica.Math.sin(pi * (time - startTime) + phase));
      signal[3] = offset + (if time < startTime then 0 else amplitude[3] * Modelica.Math.cos(pi * (time - startTime) + phase));
      //  signal = {x,y,z};
    end TrajectoryPosition1;
    block CommandLaws
      import SI = Modelica.SIunits;
      parameter Real f[3] = {0,0,1.0} "thrust";
      parameter Real M[3] = {0,0,0} "Moments";
      /*  output Avionics.Interfaces.se3.commandLaws Laws annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
      output Avionics.Types.se3CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      //zeros(3) = Laws.f + f;
      //zeros(3) = Laws.M + M;
      Laws.f = -f;
      Laws.M = -M;
      // Laws.x = zeros(3);
      // Laws.n = zeros(3);
      annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-59.4495,-49.1743},{62.3853,48.4404}}, textString = "Commmand Laws", textStyle = {TextStyle.Bold})}));
    end CommandLaws;
    model laws
      import SI = Modelica.SIunits;
      Modelica.Mechanics.Rotational.Sources.Torque torque[3] annotation(Placement(visible = true, transformation(origin = {-69.0498,-27.7582}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sources.Force force[3] annotation(Placement(visible = true, transformation(origin = {-67.1947,28.7967}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      parameter SI.Position x[3] = {0,0,0} "Positions";
      parameter SI.Force f[3] = {0,0,0} "Forces";
      parameter SI.Angle n[3] = {0,0,0} "Angles";
      parameter SI.MomentOfForce M[3] = {0,0,0} "Moments";
      Modelica.Blocks.Interfaces.RealInput f_a[3] annotation(Placement(visible = true, transformation(origin = {-104.665,30.3207}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-98.5423,29.7376}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput M_a[3] annotation(Placement(visible = true, transformation(origin = {-102.624,-28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-97.6676,-28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b M_b[3] annotation(Placement(visible = true, transformation(origin = {99.4169,-32.9446}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.4169,-32.9446}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Translational.Interfaces.Flange_b f_b[3] annotation(Placement(visible = true, transformation(origin = {99.4169,28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.4169,28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      // output Avionics.Interfaces.se3.Flange_b_Laws Lawsy;
    equation
      //  connect(fixed1.flange,force1.support) annotation(Line(points = {{12.828,-50.1458},{12.555,-50.1458},{12.555,-39.2898},{12.555,-39.2898}}));
      //connect(fixed_f.flange,force.support);
      connect(f_a,force.f);
      connect(force.flange,f_b);
      //  connect(fixed_f.flange,force.support);
      connect(M_a,torque.tau);
      connect(torque.flange,M_b);
      //
      //	Laws.f = force.f;
      //	Laws.M = torque.tau;
      //	Laws.f = f_b;
      //	Laws.M = M_b;
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end laws;
    model chose_M
      Modelica.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin = {-102.915,4.08163}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-95.9184,2.33236}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.Fixed fixed1 annotation(Placement(visible = true, transformation(origin = {-12.5364,-41.1079}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.Torque torque1 annotation(Placement(visible = true, transformation(origin = {-13.7026,3.207}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b annotation(Placement(visible = true, transformation(origin = {99.7085,1.16618}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.7085,1.16618}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    equation
      connect(torque1.flange,flange_b) annotation(Line(points = {{-3.70262,3.207},{99.4169,3.207},{99.4169,1.45773},{99.4169,1.45773}}));
      connect(fixed1.flange,torque1.support) annotation(Line(points = {{-12.5364,-41.1079},{-14.2857,-41.1079},{-14.2857,-7.58017},{-14.2857,-7.58017}}));
      connect(u,torque1.tau) annotation(Line(points = {{-102.915,4.08163},{-26.2391,4.08163},{-26.2391,3.49854},{-26.2391,3.49854}}));
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end chose_M;
    model chose_f
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(Placement(visible = true, transformation(origin = {100.292,-1.16618}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.4169,-2.91545}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sources.Force force1 annotation(Placement(visible = true, transformation(origin = {-2.01431,0.609594}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin = {-99.0724,1.48423}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-99.0724,1.48423}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    equation
      connect(u,force1.f) annotation(Line(points = {{-99.0724,1.48423},{-14.8423,1.48423},{-14.8423,0.371058},{-14.8423,0.371058}}));
      connect(force1.flange,flange_b) annotation(Line(points = {{7.98569,0.609594},{102.041,0.609594},{102.041,0.371058},{102.041,0.371058}}));
      annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
    end chose_f;
  end Sources;
  package Types "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    connector SO3dim
      Real z annotation(Placement(visible = true, transformation(origin = {4.0367,48.4403}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Real x annotation(Placement(visible = true, transformation(origin = {6.23853,0.733948}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Real y annotation(Placement(visible = true, transformation(origin = {20.1835,-33.0275}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    end SO3dim;
    record Pose
      import SI = Modelica.SIunits;
      SI.Position x[3] "Position";
      SI.Angle n[3] "Angle";
    end Pose;
    record se3Error
      import SI = Modelica.SIunits;
      //equation
      SI.Position x[3,1] "Position";
      SI.Velocity v[3,1] "Velocity";
      SI.Angle R[3,1];
      SI.AngularVelocity Omega[3,1] "Angular Velocity";
    end se3Error;
    record se3QuadrotorParameters "Quadrotor Parameters"
      // IV. Numerical Example
      import SI = Modelica.SIunits;
      SI.MomentOfInertia J[3,3] "Inertia matrix with respect of body-fixed frame";
      SI.Mass m "Mass";
      SI.Distance d "arm lenght";
      SI.Distance c_t_f "Coefficient of Thrust Force";
      SI.Acceleration g;
    end se3QuadrotorParameters;
    record se3CommandLaws "Command Laws"
      import SI = Modelica.SIunits;
      // flow
      SI.Force f[3];
      SI.MomentOfForce M[3];
    end se3CommandLaws;
    record se3TrackingVariableAngular
      import SI = Modelica.SIunits;
      SI.Angle R[3,3] "Angular Matrix";
      SI.AngularVelocity Omega[3,1] "Angular Velocity";
    end se3TrackingVariableAngular;
    record se3TrackingVariable
      extends se3TrackingVariableLinear;
      extends se3TrackingVariableAngular;
    end se3TrackingVariable;
    record se3TrackingVariableLinear
      import SI = Modelica.SIunits;
      SI.Position x[3,1] "Position";
      SI.Velocity v[3,1] "Velocity";
    end se3TrackingVariableLinear;
  end Types;
  package Interfaces "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    import SI = Modelica.SIunits;
    /*  connector XXX
    Real PotentialVariable;
    flow Real FlowVariable;
  end XXX;

  connector Frame2 "Connector for aerodynamic and mechanical purposes"
    import SI = Modelica.SIunits;
    SI.Pressure p "Static pressure";
    SI.Density rho "Atmosphere density";
    SI.Temperature T "Static temperature";
    //Types.TransferMatrix Tmat "Transfer matrix from local to inertial system";
    SI.Position r[3] "Position of frame relative inertial frame";
    SI.AngularVelocity w[3] "Angular velocity";
    SI.Velocity v[3] "Translational velocity in earth frame";
    SI.Acceleration g[3] "Gravitational acceleration";
    flow SI.Force[3] F "Force";
    flow SI.Torque[3] M "Torque";
    flow SI.Mass m "Mass";
    //flow Real[3] md(unit="kg.m") "Product mass*position (m*d) for calc CoG";
    flow SI.MomentOfInertia[3,3] I "Mass moment of inertia";
  end Frame2;
*/
    connector se3TrackConnector
      /* //extends Frame;
  import SI = Modelica.SIunits;
  SI.Position x[3,1] "Position";
  SI.Velocity v[3,1] "Velocity";
  SI.Angle R[3,3] "Angular Matrix";
  SI.AngularVelocity Omega[3,1] "Angular Velocity";
*/
      extends Avionics.Types.se3TrackingVariable;
      //   annotation(Diagram(), Icon());
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
    end se3TrackConnector;
    connector se3PoseConnector
      extends Avionics.Types.Pose;
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
    end se3PoseConnector;
    /*
connector se3HeadingConnector
  extends Avionics.Types.Pose;
  annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
end se3HeadingConnector;
*/
    /*
connector se3PositionConnector
  extends Avionics.Types.Pose;
  annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
end se3PositionConnector;
*/
    //connector bla = input SI.Angle[3]  ;
    //  connector se3AttitudeConnectorOut = input SI.Angle[3];
    //  connector se3PositionConnectorOut = input SI.Position[3];
    connector se3QuadrotorParamsConnector
      extends Avionics.Types.se3QuadrotorParameters;
      //   annotation(Diagram(), Icon());
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
    end se3QuadrotorParamsConnector;
    partial block MO "Multiple Output continuous control block"
      import Interfaces = Modelica.Blocks.Interfaces;
      extends Modelica.Blocks.Icons.Block;
      parameter Integer nout(min = 1) = 1 "Number of outputs";
      Interfaces.RealOutput signal[nout] "Connector of Real output signals" annotation(Placement(transformation(extent = {{100,-10},{120,10}}, rotation = 0)));
      annotation(Documentation(info = "<html>
<p>
Block has one continuous Real output signal vector.
</p>
</html>"));
    end MO;
    connector se3AttitudeConnector
      Modelica.SIunits.Angle n[3];
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
    end se3AttitudeConnector;
    connector se3PositionConnector
      Modelica.SIunits.Position x[3];
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
    end se3PositionConnector;
    partial block MI "Multiple Input continuous control block"
      import Interfaces = Modelica.Blocks.Interfaces;
      extends Modelica.Blocks.Icons.Block;
      parameter Integer nin(min = 1) = 1 "Number of inputs";
      Interfaces.RealInput u[nin] "Connector of Real input signals" annotation(Placement(transformation(extent = {{100,-10},{120,10}}, rotation = 0)));
      annotation(Documentation(info = "<html>
<p>
Block has one continuous Real input signal vector.
</p>
</html>"));
    end MI;
    connector MI3
      import Interfaces = Modelica.Blocks.Interfaces;
      extends Modelica.Blocks.Icons.Block;
      parameter Integer nin(min = 1) = 3 "Number of inputs";
      Interfaces.RealInput u[nin] "Connector of Real input signals" annotation(Placement(transformation(extent = {{100,-10},{120,10}}, rotation = 0)));
    end MI3;
    connector MO3
      import Interfaces = Modelica.Blocks.Interfaces;
      extends Modelica.Blocks.Icons.Block;
      parameter Integer nout(min = 1) = 3 "Number of inputs";
      Interfaces.RealOutput s[nout] "Connector of Real input signals" annotation(Placement(visible = true, transformation(origin = {-90,1.46789}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-90,1.46789}, extent = {{-10,-10},{10,10}}, rotation = 0)));
      annotation(Diagram(), Icon());
    end MO3;
    connector se3CommandLaws
      extends Avionics.Types.se3CommandLaws;
      //   annotation(Diagram(), Icon());
      annotation(defaultComponentName = "Command Laws (f,M)", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
    end se3CommandLaws;
    package se3
      connector trackVariables
        extends Avionics.Types.se3TrackingVariable;
        /*  SI.Position s[3] "Absolute position of flange";
  flow SI.Force f[3] "Cut force directed into flange";
  SI.Angle phi[3] "Absolute rotation angle of flange";
  flow SI.Torque tau[3] "Cut torque in the flange";
  */
        // flow SI.Force f[3];
        // flow SI.MomentOfForce M[3];
      end trackVariables;
      connector commandLaws
        import SI = Modelica.SIunits;
        SI.Force f[3];
        SI.MomentOfForce M[3];
      end commandLaws;
      connector commandLawsFlow
        import SI = Modelica.SIunits;
        /*SI.Position xx[3] "Position";
  SI.Angle nn[3] "Angle";
  flow SI.Force f[3];
  flow SI.MomentOfForce M[3];
*/
        extends Flange_b_Laws;
      end commandLawsFlow;
      model Flange_b_Laws
        // s,f ; phi,tau
        //extends Modelica.Mechanics.Translational.Interfaces.Flange_b[3];
        //f[3];
        //extends Modelica.Mechanics.Rotational.Interfaces.Flange_b;
        //M[3];
        SI.Position s[3] "Absolute position of flange";
        flow SI.Force f[3] "Cut force directed into flange";
        SI.Angle phi[3] "Absolute rotation angle of flange";
        flow SI.Torque tau[3] "Cut torque in the flange";
        annotation(defaultComponentName = "Command Laws (f,M)", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
      end Flange_b_Laws;
      model Flange_a_Laws
        // s,f ; phi,tau
        // extends Modelica.Mechanics.Translational.Interfaces.Flange_a;
        //f[3];
        // extends Modelica.Mechanics.Rotational.Interfaces.Flange_a;
        //M[3];
        SI.Position s[3] "Absolute position of flange";
        flow SI.Force f[3] "Cut force directed into flange";
        SI.Angle phi[3] "Absolute rotation angle of flange";
        flow SI.Torque tau[3] "Cut torque in the flange";
        annotation(defaultComponentName = "Command Laws (f,M)", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
      end Flange_a_Laws;
    end se3;
  end Interfaces;
end Avionics;


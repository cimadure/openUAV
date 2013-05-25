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
      annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
      Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-63.8532,38.8991}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.7523,1.10092}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.se3Controller se3controller1 annotation(Placement(visible = true, transformation(origin = {-11.0092,20.3694}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Bodies.SE3Dynamics_v1 se3dynamics_v11 annotation(Placement(visible = true, transformation(origin = {33.3945,21.4703}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      Avionics.Controller.se3QuadrotorParameters se3quadrotorparameters1 annotation(Placement(visible = true, transformation(origin = {-12.3898,47.7571}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      connect(trajectoryposition11.signal,se3controller1.x_d);
      connect(trajectoryheading11.signal,se3controller1.b1_d);
      connect(se3controller1.Laws,se3dynamics_v11.Laws);
      connect(se3dynamics_v11.TV,se3controller1.TV);
      connect(se3quadrotorparameters1.P,se3controller1.P);
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
      //  extends partialModel(m = m_param, g = g_param, J = J_param, x(start = x_start), v(start = v_start, fixed = true), Omega(start = Omega_start));
      //extends partialModel;
      extends GenericDynamics;
      //extends Avionics.Bodies.se3Body;
      //(x(start = x_start), v(start = v_start, fixed = true), Omega(start = Omega_start));
      constant Real e3[3,1] = [0;0;1];
      SI.Force f;
      //(start=f_start);
      SI.MomentOfForce M[3,1];
      //(start=M_start);
      SI.Angle R[3,3](start = Avionics.Functions.Rxyz(n_start[1,1], n_start[2,1], n_start[3,1]));
      // "Start Angular Matrix"
      annotation(experiment(StartTime = 0.0, StopTime = 50.0, Tolerance = 0.000001), Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-34.4954,23.8532},{34.8623,-15.0459}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
      annotation(Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-38.5321,61.2844},{33.0275,-33.3945}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
      //
      parameter String name = "GTC SE(3) model";
      parameter SI.Position x_start[3,1] = [0;0;0] "Start Positions";
      parameter SI.Velocity v_start[3,1] = [0;0;0] "Start Velocityies";
      //  SI.Angle R_start[3,3] =
      parameter SI.Angle n_start[3,1] = [0;0;0] "Start Angles";
      parameter SI.AngularVelocity Omega_start[3,1] = [0;0;0] "Start Angular Velocity";
      parameter SI.Mass m_param = 8.34;
      // 4.34
      parameter SI.MomentOfInertia J_param[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      parameter SI.Acceleration g_param = Modelica.Constants.g_n;
      //  parameter SI.Force f_start = 0 "initial force";
      //  parameter SI.MomentOfForce M_start[3,1]  =[0;0;0] "initial moments";
      input Avionics.Interfaces.CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      output Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      input Avionics.Interfaces.se3QuadrotorParamsConnector P "quadrotor Parameter" annotation(Placement(visible = true, transformation(origin = {-11.7431,96.1468}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-11.7431,96.1468}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    equation
      //  m = P.m;
      //  g = P.g;
      //  J = P.J;
      //
      f = Laws.f;
      vector(M) = Laws.M;
      //
      // LINEAR
      P.m * a = P.m * P.g * e3 - f * R * e3;
      a = der(v);
      v = der(x);
      // ANGULAR
      vector(P.J * dOmega) + cross(vector(Omega), vector(P.J * Omega)) = M[:,1];
      dOmega = der(Omega);
      R * skew(Omega[:,1]) = der(R);
      //der(R) = R * skew(Omega); // original
      // Ajout
      n = Avionics.Functions.vex(R);
      TV.x = x;
      TV.v = v;
      TV.R = R;
      TV.Omega = Omega;
    end SE3Dynamics_v1;
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
      model se3AttitudeTracking
        import SI = Modelica.SIunits;
        SI.Angle b1_d[3,1];
        SI.Angle b2_d[3,1];
        input SI.Angle b3_d[3,1];
        SI.Angle b23_d[3,1];
        SI.Angle R_d[3,3] "Angular Matrix";
        //  SI.Position e_x[3,1] "Position";
        // SI.Velocity e_v[3,1] "Velocity";
        se3Error e;
      equation
        vector(b23_d) = cross(vector(b2_d), vector(b3_d));
        vector(b2_d) = cross(vector(b3_d), vector(b1_d)) / Modelica.Math.Vectors.norm(cross(vector(b3_d), vector(b1_d)), 2);
        R_d = [b23_d,b2_d,b3_d];
      end se3AttitudeTracking;
      block TrackingController
        TrackingError e "Error Variable";
        se3parameters k "Controller parameters";
      end TrackingController;
      //
      model se3TrajectoryTracking
        import SI = Modelica.SIunits;
        output SI.Angle b3_d[3,1];
        SI.Angle temp[3,1];
        TrackingError e "Error Variable";
        se3parameters k "Controller parameters";
      equation
        temp = -k.x * e.x - k.v * e.v - m * g * e3 + m * a;
        b3_d = temp / Modelica.Math.Vectors.norm(temp);
      end se3TrajectoryTracking;
      //
      block TrackingError
        //  loadFile(Avionics);
        // se3TrackingVariable e;
        // faire un replaceable ou un redefine pour e.R er pouvoir mettre directement ca
        extends se3Error annotation(Placement(visible = true, transformation(origin = {54.7009,2.13675}, extent = {{-12,-12},{12,12}}, rotation = 0)));
        replaceable input se3TrackingVariable B annotation(Placement(visible = true, transformation(origin = {-36.3248,34.6154}, extent = {{-12,-12},{12,12}}, rotation = 0)));
        replaceable input se3TrackingVariable W annotation(Placement(visible = true, transformation(origin = {-38.4615,-40.5983}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      equation
        x = B.x - W.x;
        v = B.v - W.v;
        R = Avionics.vex(1 / 2 * (transpose(W.R) * B.R - transpose(B.R * W.R)));
        Omega = B.Omega - transpose(B.R) * W.R * W.Omega;
      end TrackingError;
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
    end SE3;
    model se3Controller
      parameter String name = "GTC SE(3) model";
      import SI = Modelica.SIunits;
      import Interfaces = Avionics.Interfaces;
      // In
      //  input Avionics.Interfaces.se3AttitudeConnector b1_d "Heading" annotation(Placement(visible = true, transformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      //  input Avionics.Interfaces.se3PositionConnector x_d "Position" annotation(Placement(visible = true, transformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      input Modelica.Blocks.Interfaces.RealInput b1_d[3] "Heading" annotation(Placement(visible = true, transformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      input Modelica.Blocks.Interfaces.RealInput x_d[3] "Position" annotation(Placement(visible = true, transformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // input Modelica.Blocks.Interfaces.RealVectorInput x_d[3] "Position" annotation(Placement(visible = true, transformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // Out
      output Avionics.Interfaces.CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {100.45,1.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.45,1.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      // params
      parameter se3parameters k;
      // Constantes
      Avionics.Types.se3Error e "Error Variable";
      // 12
      Avionics.Types.se3TrackingVariable W;
      // 18
      annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-52.844,39.2661},{-82.5688,54.6789}}, textString = "Position", textStyle = {TextStyle.Bold}),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-52.1101,-45.1376},{-78.5321,-38.5321}}, textString = "Heading"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{48.8073,7.33945},{89.5413,-6.6055}}, textString = "Control Laws"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-19.4495,-85.8715},{28.6238,-79.266}}, textString = "Tracking Variables"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-29.7248,83.3028},{8.80734,60.1835}}, textString = "Parameters"),Text(rotation = -180, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{62.3853,2.93578},{-63.1193,-39.633}}, textString = "CONTROLLER", textStyle = {TextStyle.Bold})}));
      //  input Avionics.Interfaces.se3PoseConnector Pose_d "Pose" annotation(Placement(visible = true, transformation(origin = {-99.0826,8.44037}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.0826,8.44037}, extent = {{-12,-12},{12,12}}, rotation = 0)));
      constant Real e3[3,1] = [0.0;0.0;1.0];
      input Avionics.Interfaces.se3QuadrotorParamsConnector P "quadrotor Parameter" annotation(Placement(visible = true, transformation(origin = {-11.7431,96.1468}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-11.7431,96.1468}, extent = {{12,-12},{-12,12}}, rotation = 90)));
      input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
    protected
      // Vars
      Real temp_M[3,1];
      Real temp_M_cross[3,1];
      SI.Acceleration a[3,1] "Acceleration";
      SI.Angle b3_d[3] "";
      //SI.Angle b1_d[3] "";
      Real track_coeff[3,1];
    equation
      // k.m = 0;
      //  k.m = P.m;
      // Error
      e = Avionics.Functions.se3Error(TV, W);
      a = der(TV.v);
      //
      // Tracking Controller
      track_coeff = -k.x * e.x - k.v * e.v - P.m * P.g * e3 + P.m * a;
      b3_d = -vector(track_coeff) / Modelica.Math.Vectors.norm(vector(track_coeff));
      // Linear
      //  vector(W.x) = Pose_d.x;
      //vector(W.x) = x_d.x;
      vector(W.x) = x_d;
      W.v = der(W.x);
      //
      // Angular
      //		b1_d =
      // b1_d = Pose_d.n;
      // R_d : Attitude Tracking
      //W.R = Avionics.Functions.se3R_d(b1_d.n, b3_d);
      W.R = Avionics.Functions.se3R_d(b1_d, b3_d);
      // W.Omega
      W.Omega = Avionics.Functions.vex(Modelica.Math.Matrices.inv(der(W.R)) * W.R);
      // Control output
      // f :
      Laws.f = scalar(-transpose(track_coeff) * TV.R * e3);
      // M :
      vector(temp_M_cross) = cross(vector(TV.Omega), vector(P.J * TV.Omega));
      //Laws.M
      temp_M = -k.R * e.R - k.Omega * e.Omega + temp_M_cross - P.J * (skew(vector(TV.Omega)) * transpose(TV.R) * W.R * W.Omega - transpose(TV.R) * W.R * der(W.Omega));
      Laws.M = vector(temp_M);
    end se3Controller;
    model se3QuadrotorParameters
      //  extends Avionics.Types.se3QuadrotorParameters;
      //  import SI = Modelica.SIunits;
      annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-59.4495,-49.1743},{62.3853,48.4404}}, textString = "Parameters", textStyle = {TextStyle.Bold})}));
      output Avionics.Interfaces.se3QuadrotorParamsConnector P annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
    equation
      P.J = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      P.m = 8.34;
      P.d = 0.315;
      P.c_t_f = 0.0008004;
      P.g = Modelica.Constants.g_n;
    end se3QuadrotorParameters;
  end Controller;
  package Sources "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    model TrajectoryHeading1
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
    model TrajectoryPosition1
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
    record se3TrackingVariable
      import SI = Modelica.SIunits;
      SI.Position x[3,1] "Position";
      SI.Velocity v[3,1] "Velocity";
      // replaceable
      SI.Angle R[3,3] "Angular Matrix";
      SI.AngularVelocity Omega[3,1] "Angular Velocity";
    end se3TrackingVariable;
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
      SI.MomentOfInertia J[3,3];
      //= diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
      SI.Mass m "Mass";
      //= 8.34;
      SI.Distance d "arm lenght";
      //= 0.315;
      SI.Distance c_t_f "Coefficient of Thrust Force";
      // = 0.0008004;
      //
      SI.Acceleration g;
      // = Modelica.Constants.g_n;
    end se3QuadrotorParameters;
  end Types;
  package Interfaces "Interface definitions for the Hydraulics library"
    extends Modelica.Icons.InterfacesPackage;
    import SI = Modelica.SIunits;
    connector XXX
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
    connector se3TrackConnector
      /* //extends Frame;
  import SI = Modelica.SIunits;
  SI.Position x[3,1] "Position";
  SI.Velocity v[3,1] "Velocity";
  SI.Angle R[3,3] "Angular Matrix";
  SI.AngularVelocity Omega[3,1] "Angular Velocity";
*/
      extends Avionics.Types.se3TrackingVariable;
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
      //   annotation(Diagram(), Icon());
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
    connector se3AttitudeConnectorOut = input SI.Angle[3];
    connector se3PositionConnectorOut = input SI.Position[3];
    connector se3QuadrotorParamsConnector
      extends Avionics.Types.se3QuadrotorParameters;
      annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
      //   annotation(Diagram(), Icon());
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
    connector CommandLaws
      import SI = Modelica.SIunits;
      // flow
      SI.Force f;
      SI.MomentOfForce M[3];
      annotation(defaultComponentName = "Command Laws (f,M)", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
      //   annotation(Diagram(), Icon());
    end CommandLaws;
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
      annotation(Diagram(), Icon());
      Interfaces.RealOutput s[nout] "Connector of Real input signals" annotation(Placement(visible = true, transformation(origin = {-90,1.46789}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-90,1.46789}, extent = {{-10,-10},{10,10}}, rotation = 0)));
    end MO3;
  end Interfaces;
end Avionics;


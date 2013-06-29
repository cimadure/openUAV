within Avionics.Controller.SE3;
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


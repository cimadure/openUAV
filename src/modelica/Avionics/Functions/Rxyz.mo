within Avionics.Functions;
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


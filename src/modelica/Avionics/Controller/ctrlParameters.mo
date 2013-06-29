within Avionics.Controller;
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


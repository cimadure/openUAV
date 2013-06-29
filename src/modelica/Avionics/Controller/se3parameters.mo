within Avionics.Controller;
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


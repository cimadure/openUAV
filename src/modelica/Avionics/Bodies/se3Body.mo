within Avionics.Bodies;
class se3Body
  import SI = Modelica.SIunits;
  replaceable parameter SI.Mass m;
  replaceable parameter SI.MomentOfInertia J[3,3] "Inertia matrix with respect of body-fixed frame";
  replaceable parameter SI.Acceleration g = Modelica.Constants.g_n;
end se3Body;


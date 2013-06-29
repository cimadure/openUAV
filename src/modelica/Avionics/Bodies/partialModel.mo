within Avionics.Bodies;
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


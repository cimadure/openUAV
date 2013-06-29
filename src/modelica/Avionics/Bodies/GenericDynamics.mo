within Avionics.Bodies;
record GenericDynamics
  // SE(3)
  // Ronan Cimadure, 2013-03-04
  import SI = Modelica.SIunits;
  extends GenericLinear3dDisplacement;
  extends GenericRotational3dDisplacement;
end GenericDynamics;


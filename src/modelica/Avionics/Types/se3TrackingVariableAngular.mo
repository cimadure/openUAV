within Avionics.Types;
record se3TrackingVariableAngular
  import SI = Modelica.SIunits;
  SI.Angle R[3,3] "Angular Matrix";
  SI.AngularVelocity Omega[3,1] "Angular Velocity";
end se3TrackingVariableAngular;


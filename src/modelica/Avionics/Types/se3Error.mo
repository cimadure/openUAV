within Avionics.Types;
record se3Error
  import SI = Modelica.SIunits;
  //equation
  SI.Position x[3,1] "Position";
  SI.Velocity v[3,1] "Velocity";
  SI.Angle R[3,1];
  SI.AngularVelocity Omega[3,1] "Angular Velocity";
end se3Error;


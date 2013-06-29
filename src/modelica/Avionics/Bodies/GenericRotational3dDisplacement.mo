within Avionics.Bodies;
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


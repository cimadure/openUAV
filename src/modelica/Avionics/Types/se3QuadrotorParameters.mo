within Avionics.Types;
record se3QuadrotorParameters "Quadrotor Parameters"
  // IV. Numerical Example
  import SI = Modelica.SIunits;
  SI.MomentOfInertia J[3,3] "Inertia matrix with respect of body-fixed frame";
  SI.Mass m "Mass";
  SI.Distance d "arm lenght";
  SI.Distance c_t_f "Coefficient of Thrust Force";
  SI.Acceleration g;
end se3QuadrotorParameters;


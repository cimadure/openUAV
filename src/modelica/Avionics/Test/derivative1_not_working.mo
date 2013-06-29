within Avionics.Test;
model derivative1_not_working
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.5216,49.0501}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  // Modelica.Blocks.Continuous.Derivative derivative1 annotation(Placement(visible = true, transformation(origin = {-17.2712,49.3955}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Functions.R_derivative derivative1;
  Real R[3,3];
  Real dotR[3,3];
equation
  R = Avionics.Functions.Rxyz(trajectoryheading11.signal[1], trajectoryheading11.signal[2], trajectoryheading11.signal[3]);
  connect(R,derivative1.R);
  connect(derivative1.dotR,dotR);
  //dotR = Avionics.Functions.R_derivative(R);
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end derivative1_not_working;


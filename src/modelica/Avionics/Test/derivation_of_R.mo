within Avionics.Test;
model derivation_of_R
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-61.4853,37.3057}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Real R[3,3];
  Real R2[3,3];
  // Real derR;
  //Real derR1;
  //Real a;
  input Real tolerance = 100 * Modelica.Constants.eps "Relative tolerance of solution u";
  Avionics.Functions.Derivative derivative1 annotation(Placement(visible = true, transformation(origin = {52.8497,31.0881}, extent = {{-12,-12},{12,12}}, rotation = 0)));
initial equation
  R2 = zeros(3, 3);
equation
  R = Avionics.Functions.Rxyz(trajectoryheading11.signal[1], trajectoryheading11.signal[2], trajectoryheading11.signal[3]);
  //  derR = der(R[2,2]);
  //derR = der(trajectoryheading11.signal[1]);
  // derR1 = der(R[1,3]);
  //R2 = R;
  //a = R2[1,1];
  //derR1 = der(a);
  //R2 = der(R);
  // R2 = Modelica.Math.Nonlinear.solveOneNonlinearEquation(function der(R), -0.5, 10, tolerance);
  connect(R,derivative1.u);
  connect(derivative1.y,R2);
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end derivation_of_R;


within Avionics.Test;
model dotR
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-61.4853,37.3057}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Real R[3,3];
  Real dotR[3,3];
initial equation
  dotR = zeros(3, 3);
equation
  R = Avionics.Functions.Rxyz(trajectoryheading11.signal[1], trajectoryheading11.signal[2], trajectoryheading11.signal[3]);
  dotR = Avionics.Functions.dotR(1, 0.01, R);
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end dotR;


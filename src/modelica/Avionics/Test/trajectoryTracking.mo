within Avionics.Test;
model trajectoryTracking
  Avionics.Controller.SE3.se3TrajectoryTracking se3trajectorytracking1 annotation(Placement(visible = true, transformation(origin = {-23.1434,28.3247}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-73.2297,30.3972}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Interfaces.se3TrackConnector T annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
equation
  se3trajectorytracking1.x_d = trajectoryposition11.signal;
  T.x = [0;0;0];
  T.v = [0;0;0];
  T.R = diagonal({1,1,1});
  T.Omega = [0;0;0];
  se3trajectorytracking1.TV = T;
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end trajectoryTracking;


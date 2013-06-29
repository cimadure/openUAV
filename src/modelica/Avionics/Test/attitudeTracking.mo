within Avionics.Test;
model attitudeTracking
  import SI = Modelica.SIunits;
  SI.Angle b3_d[3];
  Avionics.Controller.SE3.se3AttitudeTracking se3attitudetracking1 annotation(Placement(visible = true, transformation(origin = {-0.0000303352,-26.7229}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-39.0999,0.843849}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Interfaces.se3TrackConnector T annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
equation
  se3attitudetracking1.b1_d = trajectoryheading11.signal;
  b3_d = {0,0,0};
  T.x = [0;0;0];
  T.v = [0;0;0];
  T.R = diagonal({1,1,1});
  T.Omega = [0;0;0];
  se3attitudetracking1.TV = T;
  se3attitudetracking1.b3_d = b3_d;
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end attitudeTracking;


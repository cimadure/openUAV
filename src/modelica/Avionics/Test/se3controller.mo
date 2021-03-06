within Avionics.Test;
model se3controller
  extends Modelica.Icons.Example;
  Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-63.8532,38.8991}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.7523,1.10092}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Controller.se3Controller se3controller1 annotation(Placement(visible = true, transformation(origin = {-10.7081,21.0708}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Types.se3TrackingVariable T annotation(Placement(visible = true, transformation(origin = {-9.84529,-11.5331}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  //  connect(trajectoryposition11.signal,se3controller1.x_d);
  //  connect(trajectoryheading11.signal,se3controller1.b1_d);
  //   connect(T,se3controller1.TV);
  se3controller1.x_d = trajectoryposition11.signal;
  se3controller1.b1_d = trajectoryheading11.signal;
  T.x = [0;0;0];
  T.v = [0;0;0];
  T.R = diagonal({1,1,1});
  T.Omega = [0;0;0];
  se3controller1.TV = T;
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end se3controller;


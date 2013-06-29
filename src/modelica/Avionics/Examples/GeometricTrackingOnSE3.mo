within Avionics.Examples;
model GeometricTrackingOnSE3
  extends Modelica.Icons.Example;
  Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-63.8532,38.8991}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-62.7523,1.10092}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Controller.se3Controller se3controller1 annotation(Placement(visible = true, transformation(origin = {-10.7081,21.0708}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Bodies.SE3Dynamics_v1 se3dynamics_v11 annotation(Placement(visible = true, transformation(origin = {24.4477,20.9407}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  /*  connect(trajectoryposition11.signal,se3controller1.x_d);
  connect(trajectoryheading11.signal,se3controller1.b1_d);
  connect(se3controller1.Laws,se3dynamics_v11.Laws);
  connect(se3dynamics_v11.TV,se3controller1.TV);
*/
  connect(trajectoryposition11.signal,se3controller1.x_d);
  connect(trajectoryheading11.signal,se3controller1.b1_d);
  connect(se3controller1.Laws,se3dynamics_v11.Laws);
  connect(se3dynamics_v11.TV,se3controller1.TV);
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end GeometricTrackingOnSE3;


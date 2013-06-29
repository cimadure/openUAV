within Avionics.Test;
model attitudeTracking2
  import SI = Modelica.SIunits;
  import SIm = Modelica.Math;
  SI.Angle b3_d annotation(Placement(visible = true, transformation(origin = {41.3502,3.09423}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-66.1041,75.668}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Controller.SE3.se3TrajectoryTracking se3trajectorytracking1 annotation(Placement(visible = true, transformation(origin = {-28.4107,75.3868}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Controller.SE3.se3AttitudeTracking se3attitudetracking1 annotation(Placement(visible = true, transformation(origin = {1.96903,41.9128}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Interfaces.se3TrackConnector T annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
equation
  /*se3attitudetracking1.b1_d = trajectoryheading11.signal;
  T.x = [0;0;0];
  T.v = [0;0;0];
  T.R = diagonal({1,1,1});
  T.Omega = [0;0;0];
  se3attitudetracking1.TV = T;
  connect(se3trajectorytracking1.b3_d,se3attitudetracking1.b3_d);
*/
end attitudeTracking2;


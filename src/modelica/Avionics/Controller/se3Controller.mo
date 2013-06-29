within Avionics.Controller;
model se3Controller
  parameter String name = "GTC SE(3) controller model";
  import SI = Modelica.SIunits;
  import Interfaces = Avionics.Interfaces;
  //
  // In
  input Modelica.Blocks.Interfaces.RealInput b1_d[3] "Heading" annotation(Placement(visible = true, transformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.8165,-42.2018}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  input Modelica.Blocks.Interfaces.RealInput x_d[3] "Position" annotation(Placement(visible = true, transformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.4495,45.5045}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  // Out
  output Avionics.Interfaces.se3CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {100.45,1.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.45,1.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  //
  parameter SI.Mass m = 8.34;
  // 4.34
  parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
  parameter SI.Acceleration g = Modelica.Constants.g_n;
  input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
protected
  Avionics.Controller.SE3.se3TrajectoryTracking se3trajectorytracking1 annotation(Placement(visible = true, transformation(origin = {-57.1027,43.3193}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Controller.SE3.se3AttitudeTracking se3attitudetracking1 annotation(Placement(visible = true, transformation(origin = {-12.9395,31.7862}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  //  connect(x_d,se3trajectorytracking1.x_d);
  //  connect(b1_d,se3attitudetracking1.b1_d);
  se3trajectorytracking1.x_d = x_d;
  se3attitudetracking1.b1_d = b1_d;
  connect(se3trajectorytracking1.b3_d,se3attitudetracking1.b3_d);
  connect(TV,se3trajectorytracking1.TV);
  connect(TV,se3attitudetracking1.TV);
  Laws.f = se3trajectorytracking1.f;
  Laws.M = se3attitudetracking1.M;
  //Laws.f = 3.4567;
  //Laws.M = {1.1,1.2,0.3};
  annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-52.844,39.2661},{-82.5688,54.6789}}, textString = "Position", textStyle = {TextStyle.Bold}),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-52.1101,-45.1376},{-78.5321,-38.5321}}, textString = "Heading"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{48.8073,7.33945},{89.5413,-6.6055}}, textString = "Control Laws"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-19.4495,-85.8715},{28.6238,-79.266}}, textString = "Tracking Variables"),Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-29.7248,83.3028},{8.80734,60.1835}}, textString = "Parameters"),Text(rotation = -180, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{62.3853,2.93578},{-63.1193,-39.633}}, textString = "CONTROLLER", textStyle = {TextStyle.Bold})}));
end se3Controller;


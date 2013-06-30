within Avionics.Test;
model forces2
  Avionics.Sources.TrajectoryPosition1 trajectoryposition11(amplitude = {0.0,0.0,0}) annotation(Placement(visible = true, transformation(origin = {-73.4694,-20.8322}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.laws laws1 annotation(Placement(visible = true, transformation(origin = {20.3656,15.1603}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Bodies.dynamics dynamics1(m = 8) annotation(Placement(visible = true, transformation(origin = {53.3229,15.3829}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.sinForce sinforce1(offset = -dynamics1.m * Modelica.Constants.g_n, phase = Modelica.Constants.pi, amplitude = {0,0,100}) annotation(Placement(visible = true, transformation(origin = {-72.4324,42.7027}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(sinforce1.signal,laws1.f_a) annotation(Line(points = {{-61.4324,42.7027},{10.2703,42.7027},{10.2703,18.9189},{10.2703,18.9189}}));
  connect(laws1.M_b,dynamics1.M_a) annotation(Line(points = {{30.3073,11.8658},{40.5512,11.8658},{63.2646,12.7542},{63.2646,12.0884}}));
  connect(laws1.f_b,dynamics1.f_a) annotation(Line(points = {{30.3073,17.9883},{41.5955,17.9883},{63.2646,18.0659},{63.2646,18.2109}}));
  connect(trajectoryposition11.signal,laws1.M_a) annotation(Line(points = {{-62.4694,-20.8322},{-39.6236,7.7392},{10.1032,11.1662},{10.1032,12.3323}}));
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end forces2;


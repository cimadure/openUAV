within Avionics.Test;
model forces
  Avionics.Sources.TrajectoryPosition1 trajectoryposition11(amplitude = {0.0,0.0,0}) annotation(Placement(visible = true, transformation(origin = {-73.4694,-20.8322}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.laws laws1 annotation(Placement(visible = true, transformation(origin = {20.3656,15.1603}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Bodies.dynamics dynamics1 annotation(Placement(visible = true, transformation(origin = {52.2418,14.3018}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.stepForce stepforce1(amplitude = {0.4,100,-1}, offset = -80) annotation(Placement(visible = true, transformation(origin = {-74.0541,24.8649}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(stepforce1.signal,laws1.f_a) annotation(Line(points = {{-63.0541,24.8649},{9.18919,24.8649},{9.18919,18.1341},{10.5114,18.1341}}));
  connect(laws1.f_b,dynamics1.f_a) annotation(Line(points = {{30.3073,17.9883},{41.5955,17.9883},{62.1835,18.0659},{62.1835,17.1298}}));
  connect(laws1.M_b,dynamics1.M_a) annotation(Line(points = {{30.3073,11.8658},{40.5512,11.8658},{62.1835,12.7542},{62.1835,11.0073}}));
  connect(trajectoryposition11.signal,laws1.M_a) annotation(Line(points = {{-62.4694,-20.8322},{-39.6236,7.7392},{10.1032,11.1662},{10.1032,12.3323}}));
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end forces;


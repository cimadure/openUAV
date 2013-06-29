within Avionics.Test;
model commandDynamics
  Avionics.Sources.TrajectoryPosition1 trajectoryposition11 annotation(Placement(visible = true, transformation(origin = {-70.3432,31.8039}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.TrajectoryHeading1 trajectoryheading11 annotation(Placement(visible = true, transformation(origin = {-70.5158,2.61978}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Bodies.quadcopter quadcopter1 annotation(Placement(visible = true, transformation(origin = {11.3524,28.5993}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.laws laws1 annotation(Placement(visible = true, transformation(origin = {-21.8663,32.2933}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(trajectoryheading11.signal,laws1.M_a) annotation(Line(points = {{-59.5158,2.61978},{-40.3883,2.61978},{-32.1287,25.7612},{-32.1287,29.4653}}));
  connect(trajectoryposition11.signal,laws1.f_a) annotation(Line(points = {{-59.3432,31.8039},{-39.9517,31.8039},{-32.3328,31.2191},{-32.3328,35.3254}}));
  connect(laws1.f_b,quadcopter1.f_a);
  connect(laws1.M_b,quadcopter1.M_a);
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end commandDynamics;


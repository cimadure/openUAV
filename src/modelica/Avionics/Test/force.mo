within Avionics.Test;
model force
  Avionics.Sources.chose_f teee1 annotation(Placement(visible = true, transformation(origin = {-17.0686,11.1317}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Blocks.Sources.Sine sine1 annotation(Placement(visible = true, transformation(origin = {-66.0482,10.7607}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass mass1(m = 1) annotation(Placement(visible = true, transformation(origin = {24.8649,11.7432}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(teee1.flange_b,mass1.flange_a) annotation(Line(points = {{-7.0394,11.0151},{14.6069,11.0151},{14.8649,10.8108},{14.8649,11.7432}}));
  connect(sine1.y,teee1.u) annotation(Line(points = {{-55.0482,10.7607},{-27.4583,10.3896},{-27.4583,10.986},{-27.1269,10.986}}));
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end force;


within Avionics.Sources;
model chose_f
  Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(Placement(visible = true, transformation(origin = {100.292,-1.16618}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.4169,-2.91545}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sources.Force force1 annotation(Placement(visible = true, transformation(origin = {-2.01431,0.609594}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin = {-99.0724,1.48423}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-99.0724,1.48423}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(u,force1.f) annotation(Line(points = {{-99.0724,1.48423},{-14.8423,1.48423},{-14.8423,0.371058},{-14.8423,0.371058}}));
  connect(force1.flange,flange_b) annotation(Line(points = {{7.98569,0.609594},{102.041,0.609594},{102.041,0.371058},{102.041,0.371058}}));
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end chose_f;


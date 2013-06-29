within Avionics.Sources;
model chose_M
  Modelica.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin = {-102.915,4.08163}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-95.9184,2.33236}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Components.Fixed fixed1 annotation(Placement(visible = true, transformation(origin = {-12.5364,-41.1079}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.Torque torque1 annotation(Placement(visible = true, transformation(origin = {-13.7026,3.207}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b annotation(Placement(visible = true, transformation(origin = {99.7085,1.16618}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.7085,1.16618}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(torque1.flange,flange_b) annotation(Line(points = {{-3.70262,3.207},{99.4169,3.207},{99.4169,1.45773},{99.4169,1.45773}}));
  connect(fixed1.flange,torque1.support) annotation(Line(points = {{-12.5364,-41.1079},{-14.2857,-41.1079},{-14.2857,-7.58017},{-14.2857,-7.58017}}));
  connect(u,torque1.tau) annotation(Line(points = {{-102.915,4.08163},{-26.2391,4.08163},{-26.2391,3.49854},{-26.2391,3.49854}}));
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end chose_M;


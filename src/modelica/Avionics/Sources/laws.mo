within Avionics.Sources;
model laws
  import SI = Modelica.SIunits;
  Modelica.Mechanics.Rotational.Sources.Torque torque[3] annotation(Placement(visible = true, transformation(origin = {-69.0498,-27.7582}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sources.Force force[3] annotation(Placement(visible = true, transformation(origin = {-67.1947,28.7967}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  parameter SI.Position x[3] = {0,0,0} "Positions";
  parameter SI.Force f[3] = {0,0,0} "Forces";
  parameter SI.Angle n[3] = {0,0,0} "Angles";
  parameter SI.MomentOfForce M[3] = {0,0,0} "Moments";
  Modelica.Blocks.Interfaces.RealInput f_a[3] annotation(Placement(visible = true, transformation(origin = {-104.665,30.3207}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-98.5423,29.7376}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput M_a[3] annotation(Placement(visible = true, transformation(origin = {-102.624,-28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-97.6676,-28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b M_b[3] annotation(Placement(visible = true, transformation(origin = {99.4169,-32.9446}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.4169,-32.9446}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Interfaces.Flange_b f_b[3] annotation(Placement(visible = true, transformation(origin = {99.4169,28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {99.4169,28.2799}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  // output Avionics.Interfaces.se3.Flange_b_Laws Lawsy;
equation
  //  connect(fixed1.flange,force1.support) annotation(Line(points = {{12.828,-50.1458},{12.555,-50.1458},{12.555,-39.2898},{12.555,-39.2898}}));
  //connect(fixed_f.flange,force.support);
  connect(f_a,force.f);
  connect(force.flange,f_b);
  //  connect(fixed_f.flange,force.support);
  connect(M_a,torque.tau);
  connect(torque.flange,M_b);
  //
  //	Laws.f = force.f;
  //	Laws.M = torque.tau;
  //	Laws.f = f_b;
  //	Laws.M = M_b;
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end laws;


within Avionics.Bodies;
model quadcopter
  import SI = Modelica.SIunits;
  import Constants = Modelica.Constants;
  extends GenericDynamics;
  constant Real e3[3,1] = [0;0;1];
  //
  //
  parameter String name = "GTC SE(3) model";
  parameter SI.Position x_start[3,1] = [0;0;0] "Start Positions";
  parameter SI.Velocity v_start[3,1] = [0;0;0] "Start Velocityies";
  parameter SI.Angle n_start[3,1] = [0;0;0] "Start Angles";
  parameter SI.AngularVelocity Omega_start[3,1] = [0;0;0] "Start Angular Velocity";
  parameter SI.Mass m = 8.34;
  // 4.34
  parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
  parameter SI.Acceleration g = Modelica.Constants.g_n;
  //
  /*  input Avionics.Interfaces.se3.commandLaws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
  /*  input Avionics.Interfaces.se3.Flange_a_Laws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
  output Avionics.Interfaces.se3.trackVariables TV annotation(Placement(visible = true, transformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Modelica.Mechanics.Translational.Interfaces.Flange_a f_a[3] annotation(Placement(visible = true, transformation(origin = {-95.3405,40.1434}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-95.3405,40.1434}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b M_a[3] annotation(Placement(visible = true, transformation(origin = {-100.358,-45.1613}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-100.358,-45.1613}, extent = {{-10,-10},{10,10}}, rotation = 0)));
protected
  // SI.Force f;
  //SI.MomentOfForce M[3,1];
  SI.Angle R[3,3] "Start Angular Matrix";
initial equation
  R = Avionics.Functions.Rxyz(M_a[1].phi, M_a[2].phi, M_a[3].phi);
  Omega = Omega_start;
  x = [f_a[1].s;f_a[2].s;f_a[3].s];
  v = v_start;
equation
  // flow
  //vector(M) = M_a.tau;
  //M = [M_a[1].tau;M_a[2].tau;M_a[3].tau];
  //
  //  zeros(3) = f_a.f + TV.f;
  //  zeros(3) = M_a.tau + vector(M);
  // f_a[1] = 0; f
  //  TV.f = -f_a.f;
  //  TV.M = -M_a.tau;
  // LINEAR
  m * a = m * g * e3 - f_a[3].f * R * e3;
  a = der(v);
  v = der(x);
  // ANGULAR
  vector(J * dOmega) + cross(vector(Omega), vector(J * Omega)) = M_a.tau;
  dOmega = der(Omega);
  R * skew(Omega[:,1]) = der(R);
  //der(R) = R * skew(Omega); // original
  // Add
  n = Avionics.Functions.vex(R);
  TV.x = x;
  TV.v = v;
  TV.R = R;
  TV.Omega = Omega;
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
  annotation(Diagram(), Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2}), graphics = {Text(origin = {2.16216,-8.64865}, lineColor = {0,0,255}, extent = {{-38.5321,61.2844},{33.0275,-33.3945}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
end quadcopter;


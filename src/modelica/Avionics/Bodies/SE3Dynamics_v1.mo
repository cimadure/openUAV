within Avionics.Bodies;
model SE3Dynamics_v1
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
  input Avionics.Interfaces.se3CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-99.55,3.77064}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  output Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {100.55,1.83487}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  //  input Avionics.Interfaces.se3QuadrotorParamsConnector P "quadrotor Parameter" annotation(Placement(visible = true, transformation(origin = {-11.7431,96.1468}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-11.7431,96.1468}, extent = {{12,-12},{-12,12}}, rotation = 90)));
protected
  SI.Force f;
  SI.MomentOfForce M[3,1];
  SI.Angle R[3,3] "Start Angular Matrix";
initial equation
  R = Avionics.Functions.Rxyz(n_start[1,1], n_start[2,1], n_start[3,1]);
  Omega = Omega_start;
  x = x_start;
  v = v_start;
equation
  f = Laws.f;
  vector(M) = Laws.M;
  // LINEAR
  m * a = m * g * e3 - f * R * e3;
  a = der(v);
  v = der(x);
  // ANGULAR
  vector(J * dOmega) + cross(vector(Omega), vector(J * Omega)) = M[:,1];
  dOmega = der(Omega);
  R * skew(Omega[:,1]) = der(R);
  //der(R) = R * skew(Omega); // original
  // Add
  n = Avionics.Functions.vex(R);
  TV.x = x;
  TV.v = v;
  TV.R = R;
  TV.Omega = Omega;
  annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-34.4954,23.8532},{34.8623,-15.0459}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
  annotation(Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-38.5321,61.2844},{33.0275,-33.3945}}, textString = "Dynamics T. Lee", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
end SE3Dynamics_v1;


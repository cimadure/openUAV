within Avionics.Controller.SE3;
model se3TrajectoryTracking
  import SI = Modelica.SIunits;
  se3parameters k "Controller parameters";
  parameter SI.Mass m = 8.34;
  parameter SI.Acceleration g = Modelica.Constants.g_n;
  input Modelica.Blocks.Interfaces.RealInput x_d[3] "Position" annotation(Placement(visible = true, transformation(origin = {-96.9179,14.8435}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-96.9179,14.8435}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  output SI.Angle b3_d[3] "";
  output SI.Force f;
  input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
protected
  // Vars
  Avionics.Types.se3TrackingVariableLinear e "Error Variable";
  Real track_coeff[3,1];
  SI.Velocity v_d[3,1] "desired Velocity";
  SI.Acceleration a[3,1] "Acceleration";
  constant Real e3[3,1] = [0.0;0.0;1.0];
  Real temp[3,1];
equation
  vector(v_d) = der(x_d);
  a = der(TV.v);
  vector(e.x) = vector(TV.x) - x_d;
  //e.x = TV.x - transpose(x_d);
  e.v = TV.v - v_d;
  // Tracking Controller
  track_coeff = -k.x * e.x - k.v * e.v - m * g * e3 + m * a;
  b3_d = -vector(track_coeff) / Modelica.Math.Matrices.frobeniusNorm(track_coeff);
  // b3_d = -vector(track_coeff) / Modelica.Math.Vectors.norm(vector(track_coeff));
  // Control output
  // f :
  //  f = scalar(-transpose(track_coeff) * TV.R * e3);
  temp = TV.R * e3;
  f = -sum(transpose(track_coeff) * temp);
  //  b3_d = {1.1,1.2,0.3};
  annotation(Icon(), Diagram());
end se3TrajectoryTracking;


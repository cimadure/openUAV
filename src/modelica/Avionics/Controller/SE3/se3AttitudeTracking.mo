within Avionics.Controller.SE3;
model se3AttitudeTracking
  import SI = Modelica.SIunits;
  input SI.Angle b1_d[3];
  input SI.Angle b3_d[3];
  output SI.MomentOfForce M[3];
  parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
  se3parameters k "Controller parameters";
  input Avionics.Interfaces.se3TrackConnector TV annotation(Placement(visible = true, transformation(origin = {-0.834862,-97.8803}, extent = {{-12,12},{12,-12}}, rotation = -90), iconTransformation(origin = {-0.834862,-97.8803}, extent = {{12,-12},{-12,12}}, rotation = 90)));
protected
  Avionics.Functions.Derivative dotR;
  // Avionics.Controller.SE3.R_d r_d1 ;
  //
  Real temp_M[3,1];
  Real temp_M_cross[3,1];
  SI.Angle e_R[3,1];
  SI.AngularVelocity e_Omega[3,1] "Angular Velocity";
  Real R_d[3,3] "Angular Matrix";
  Real R_dd[3,3];
  Real der_R_d[3,3];
  Real inv_dotR_d[3,3];
  //   SI.Angle
  SI.AngularVelocity Omega_d[3,1] "Angular Velocity";
  //Real der_Omega_d[3,1];
  Real transp_R_d[3,3];
equation
  // Angular
  // R_d : Attitude Tracking
  /*  connect(b1_d,r_d1.b1_d);
  connect(b3_d,r_d1.b3_d);
  R_d = r_d1.R;
*/
  R_d = Avionics.Functions.se3R_d(b1_d, b3_d);
  R_dd = R_d;
  dotR.u = R_dd;
  der_R_d = dotR.y;
  //  dotR.u = R_d;
  //  der_R_d = dotR.y;
  //
  //der_R_d = der(R_d);
  //der(R_d) = der_R_d;
  // R_d3 = Modelica.Math.Matrices.inv(R_d2);
  //
  //
  transp_R_d = transpose(R_d);
  inv_dotR_d = Modelica.Math.Matrices.inv(der_R_d);
  // W.Omega
  Omega_d = Avionics.Functions.vex(inv_dotR_d * R_d);
  //Omega_d = Avionics.Functions.vex(Modelica.Math.Matrices.inv(der(R_d)) * R_d);
  e_R = Avionics.Functions.vex(1 / 2 * transp_R_d * TV.R - transpose(TV.R * R_d));
  e_Omega = TV.Omega - transpose(TV.R) * R_d * Omega_d;
  // Control output
  // M :
  vector(temp_M_cross) = cross(vector(TV.Omega), vector(J * TV.Omega));
  temp_M = -k.R * e_R - k.Omega * e_Omega + temp_M_cross - J * (skew(vector(TV.Omega)) * transpose(TV.R) * R_d * Omega_d - transpose(TV.R) * R_d * der(Omega_d));
  //der_Omega_d = der(Omega_d);
  //  temp_M = -k.R * e_R - k.Omega * e_Omega + temp_M_cross;
  //transp_R = transpose(TV.R);
  //transp_R = diagonal({1.1,1.2,0.3});
  //temp_M = -k.R * e_R - k.Omega * e_Omega + temp_M_cross - J * (skew(vector(TV.Omega)) * transp_R * R_d * Omega_d - transp_R * R_d * der_Omega_d);
  M = vector(temp_M);
  // M = {1.1,1.2,0.3};
end se3AttitudeTracking;


within Avionics.Controller.SE3;
model R_d
  import SI = Modelica.SIunits;
  input SI.Angle b1_d[3];
  input SI.Angle b3_d[3];
  output Real R[3,3] "Angular Matrix";
  // output Real der_R[3,3];
  // Angular
  // R_d : Attitude Tracking
  //R = Avionics.Functions.se3R_d(b1_d, b3_d);
protected
  SI.Angle b31_d[3] "";
  SI.Angle b23_d[3] "";
  SI.Angle b2_d[3] "";
equation
  // Angular
  // R_d : Attitude Tracking
  b31_d = cross(b3_d, b1_d);
  //b2_d = b31_d / Modelica.Math.Matrices.frobeniusNorm(b31_d);
  b2_d = b31_d / Modelica.Math.Vectors.norm(b31_d);
  b23_d = cross(b2_d, b3_d);
  //  R := [matrix(b23_d,3,1), matrix(b2_d,3,1), matrix(b3_d,3,1)];
  R = [b23_d[1],b2_d[1],b3_d[1];b23_d[2],b2_d[2],b3_d[2];b23_d[3],b2_d[3],b3_d[3]];
end R_d;


within Avionics.Functions;
function se3R_d
  import SI = Modelica.SIunits;
  input SI.Angle b1_d[3] "";
  input SI.Angle b3_d[3] "";
  output SI.Angle R[3,3] "Angular Matrix";
protected
  SI.Angle b31_d[3] "";
  SI.Angle b23_d[3] "";
  SI.Angle b2_d[3] "";
algorithm
  // Angular
  // R_d : Attitude Tracking
  b31_d:=cross(b3_d, b1_d);
  b2_d:=b31_d / Modelica.Math.Vectors.norm(b31_d);
  b23_d:=cross(b2_d, b3_d);
  //  R := [matrix(b23_d,3,1), matrix(b2_d,3,1), matrix(b3_d,3,1)];
  R:=[b23_d[1],b2_d[1],b3_d[1];b23_d[2],b2_d[2],b3_d[2];b23_d[3],b2_d[3],b3_d[3]];
end se3R_d;


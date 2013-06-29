within Avionics.Controller.SE3;
model der_R_d
  import SI = Modelica.SIunits;
  Avionics.Controller.SE3.R_d r_d1;
  input SI.Angle b1_d[3];
  input SI.Angle b3_d[3];
  output Real der_R_d[3,3];
equation
  connect(b1_d,r_d1.b1_d);
  connect(b3_d,r_d1.b3_d);
  der_R_d = der(r_d1.R);
end der_R_d;


within Avionics.Test;
model derR_d_not_working
  Avionics.Controller.SE3.der_R_d der_r_d1 annotation(Placement(visible = true, transformation(origin = {-45.8509,15.7525}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  der_r_d1.b1_d = {0,0,0};
  der_r_d1.b3_d = {0,0,0};
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end derR_d_not_working;


within Avionics.Controller;
model se3QuadrotorParameters
  //  parameter Avionics.Types.se3QuadrotorParameters;
  import SI = Modelica.SIunits;
  parameter SI.MomentOfInertia J[3,3] = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
  parameter SI.Mass m = 8.34 "Mass";
  //
  parameter SI.Distance d = 0.315 "arm lenght";
  //
  parameter SI.Distance c_t_f = 0.0008004 "Coefficient of Thrust Force";
  //;
  //
  parameter SI.Acceleration g = Modelica.Constants.g_n "gravity";
  //  import SI = Modelica.SIunits;
  output Avionics.Interfaces.se3QuadrotorParamsConnector P annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  P.J = diagonal({0.082,0.0845,0.1377}) "Inertia matrix with respect of body-fixed frame";
  P.m = 8.34;
  P.d = 0.315;
  P.c_t_f = 0.0008004;
  P.g = Modelica.Constants.g_n;
  annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-59.4495,-49.1743},{62.3853,48.4404}}, textString = "Parameters", textStyle = {TextStyle.Bold})}));
end se3QuadrotorParameters;


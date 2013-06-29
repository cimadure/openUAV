within Avionics.Test;
model dynamics_v1
  Avionics.Sources.CommandLaws commandlaws1(f = 40.0) annotation(Placement(visible = true, transformation(origin = {-16.1605,52.5697}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Avionics.Bodies.SE3Dynamics_v1 se3dynamics_v11 annotation(Placement(visible = true, transformation(origin = {12.1097,51.9142}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  connect(commandlaws1.Laws,se3dynamics_v11.Laws);
  annotation(experiment(StartTime = 0.0, StopTime = 50.0, Tolerance = 0.000001), Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-34.4954,23.8532},{34.8623,-15.0459}}, textString = "Dynamics", fontSize = 16, fontName = "Times New Roman", textStyle = {TextStyle.Bold})}));
end dynamics_v1;


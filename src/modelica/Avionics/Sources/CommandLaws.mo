within Avionics.Sources;
block CommandLaws
  import SI = Modelica.SIunits;
  parameter Real f[3] = {0,0,1.0} "thrust";
  parameter Real M[3] = {0,0,0} "Moments";
  /*  output Avionics.Interfaces.se3.commandLaws Laws annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
*/
  output Avionics.Types.se3CommandLaws Laws annotation(Placement(visible = true, transformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {99.4495,-2.93578}, extent = {{-12,-12},{12,12}}, rotation = 0)));
equation
  //zeros(3) = Laws.f + f;
  //zeros(3) = Laws.M + M;
  Laws.f = -f;
  Laws.M = -M;
  // Laws.x = zeros(3);
  // Laws.n = zeros(3);
  annotation(Diagram(), Icon(graphics = {Text(rotation = 0, lineColor = {0,0,255}, fillColor = {0,0,0}, pattern = LinePattern.Solid, fillPattern = FillPattern.None, lineThickness = 0.25, extent = {{-59.4495,-49.1743},{62.3853,48.4404}}, textString = "Commmand Laws", textStyle = {TextStyle.Bold})}));
end CommandLaws;


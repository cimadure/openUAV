within Avionics.Examples;
model SI3MOs
  extends Modelica.Blocks.Icons.Block;
  //  Modelica.Blocks.Interfaces.RealInput u[3] annotation(Placement(visible = true, transformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  // annotation(Placement(visible = true, transformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-92.844,7.33945}, extent = {{-12,-12},{12,12}}, rotation = 0)), Diagram(), Icon());
  Modelica.Blocks.Interfaces.RealOutput x annotation(Placement(visible = true, transformation(origin = {107.523,64.5872}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {107.523,64.5872}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput y annotation(Placement(visible = true, transformation(origin = {107.156,2.20183}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {107.156,2.20183}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput z annotation(Placement(visible = true, transformation(origin = {108.257,-60.5505}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {108.257,-60.5505}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealVectorInput u[3] annotation(Placement(visible = true, transformation(origin = {-98.7156,6.6055}, extent = {{-12,-12},{12,12}}, rotation = 0), iconTransformation(origin = {-98.7156,6.6055}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  // parameter Integer nin = 3 "Number of inputs" annotation(Placement(visible = true, transformation(origin = {0,0.366972}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  // Modelica.Blocks.Interfaces.RealInput u[nin] "Connector of Real input signals" annotation(Placement(visible = true, transformation(origin = {-120,0.733945}, extent = {{-20,-20},{20,20}}, rotation = 0)));
equation
  x = u[1];
  y = u[2];
  z = u[3];
end SI3MOs;


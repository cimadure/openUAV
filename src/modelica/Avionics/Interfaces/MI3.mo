within Avionics.Interfaces;
connector MI3
  import Interfaces = Modelica.Blocks.Interfaces;
  extends Modelica.Blocks.Icons.Block;
  parameter Integer nin(min = 1) = 3 "Number of inputs";
  Interfaces.RealInput u[nin] "Connector of Real input signals" annotation(Placement(transformation(extent = {{100,-10},{120,10}}, rotation = 0)));
end MI3;


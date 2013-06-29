within Avionics.Interfaces;
connector MO3
  import Interfaces = Modelica.Blocks.Interfaces;
  extends Modelica.Blocks.Icons.Block;
  parameter Integer nout(min = 1) = 3 "Number of inputs";
  Interfaces.RealOutput s[nout] "Connector of Real input signals" annotation(Placement(visible = true, transformation(origin = {-90,1.46789}, extent = {{-10,-10},{10,10}}, rotation = 0), iconTransformation(origin = {-90,1.46789}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  annotation(Diagram(), Icon());
end MO3;


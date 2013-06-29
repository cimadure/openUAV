within Avionics.Interfaces;
partial block MI "Multiple Input continuous control block"
  import Interfaces = Modelica.Blocks.Interfaces;
  extends Modelica.Blocks.Icons.Block;
  parameter Integer nin(min = 1) = 1 "Number of inputs";
  Interfaces.RealInput u[nin] "Connector of Real input signals" annotation(Placement(transformation(extent = {{100,-10},{120,10}}, rotation = 0)));
  annotation(Documentation(info = "<html>
<p>
Block has one continuous Real input signal vector.
</p>
</html>"));
end MI;


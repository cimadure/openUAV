within Avionics.Interfaces;
partial block MO "Multiple Output continuous control block"
  import Interfaces = Modelica.Blocks.Interfaces;
  extends Modelica.Blocks.Icons.Block;
  parameter Integer nout(min = 1) = 1 "Number of outputs";
  Interfaces.RealOutput signal[nout] "Connector of Real output signals" annotation(Placement(transformation(extent = {{100,-10},{120,10}}, rotation = 0)));
  annotation(Documentation(info = "<html>
<p>
Block has one continuous Real output signal vector.
</p>
</html>"));
end MO;


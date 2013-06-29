within Avionics.Types;
record Pose
  import SI = Modelica.SIunits;
  SI.Position x[3] "Position";
  SI.Angle n[3] "Angle";
end Pose;


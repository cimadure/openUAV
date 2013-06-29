within Avionics.Bodies;
record GenericLinear3dDisplacement
  import SI = Modelica.SIunits;
  SI.Position x[3,1] "Position";
  SI.Velocity v[3,1] "Velocity";
  SI.Velocity a[3,1] "Accelertion";
end GenericLinear3dDisplacement;


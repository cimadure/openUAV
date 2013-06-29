within Avionics.Functions;
function se3Error "e = (TV,W)  [12, e.R==3]"
  import SI = Modelica.SIunits;
  input Avionics.Types.se3TrackingVariable TV;
  input Avionics.Types.se3TrackingVariable W;
  output Avionics.Types.se3Error e "Error Variable";
algorithm
  e.x:=TV.x - W.x;
  e.v:=TV.v - W.v;
  e.R:=Avionics.Functions.vex(1 / 2 * (transpose(W.R) * TV.R - transpose(TV.R * W.R)));
  e.Omega:=TV.Omega - transpose(TV.R) * W.R * W.Omega;
end se3Error;


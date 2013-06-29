within Avionics.Functions;
function dotR
  input Real k;
  input Real T;
  input Real u[3,3];
  output Real y[3,3];
protected
  Real x[3,3];
  Real dx[3,3];
algorithm
  y:=k / T * (u - x);
  //y := zeros(3,3);
  dx:=der(x);
  dx:=(u - x) / T;
  //equation
  //  der(x) = (u - x) / T;
  annotation(Icon(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})), Diagram(coordinateSystem(extent = {{-100,-100},{100,100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2,2})));
end dotR;


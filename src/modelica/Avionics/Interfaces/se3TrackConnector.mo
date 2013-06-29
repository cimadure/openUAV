within Avionics.Interfaces;
connector se3TrackConnector
  /* //extends Frame;
  import SI = Modelica.SIunits;
  SI.Position x[3,1] "Position";
  SI.Velocity v[3,1] "Velocity";
  SI.Angle R[3,3] "Angular Matrix";
  SI.AngularVelocity Omega[3,1] "Angular Velocity";
*/
  extends Avionics.Types.se3TrackingVariable;
  //   annotation(Diagram(), Icon());
  annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
end se3TrackConnector;


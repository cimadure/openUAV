within Avionics.Interfaces.se3;
model Flange_b_Laws
  // s,f ; phi,tau
  //extends Modelica.Mechanics.Translational.Interfaces.Flange_b[3];
  //f[3];
  //extends Modelica.Mechanics.Rotational.Interfaces.Flange_b;
  //M[3];
  SI.Position s[3] "Absolute position of flange";
  flow SI.Force f[3] "Cut force directed into flange";
  SI.Angle phi[3] "Absolute rotation angle of flange";
  flow SI.Torque tau[3] "Cut torque in the flange";
  annotation(defaultComponentName = "Command Laws (f,M)", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
end Flange_b_Laws;


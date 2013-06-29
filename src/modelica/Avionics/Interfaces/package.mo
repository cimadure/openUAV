within Avionics;
package Interfaces "Interface definitions for the Hydraulics library"
  extends Modelica.Icons.InterfacesPackage;
  import SI = Modelica.SIunits;
  /*  connector XXX
    Real PotentialVariable;
    flow Real FlowVariable;
  end XXX;

  connector Frame2 "Connector for aerodynamic and mechanical purposes"
    import SI = Modelica.SIunits;
    SI.Pressure p "Static pressure";
    SI.Density rho "Atmosphere density";
    SI.Temperature T "Static temperature";
    //Types.TransferMatrix Tmat "Transfer matrix from local to inertial system";
    SI.Position r[3] "Position of frame relative inertial frame";
    SI.AngularVelocity w[3] "Angular velocity";
    SI.Velocity v[3] "Translational velocity in earth frame";
    SI.Acceleration g[3] "Gravitational acceleration";
    flow SI.Force[3] F "Force";
    flow SI.Torque[3] M "Torque";
    flow SI.Mass m "Mass";
    //flow Real[3] md(unit="kg.m") "Product mass*position (m*d) for calc CoG";
    flow SI.MomentOfInertia[3,3] I "Mass moment of inertia";
  end Frame2;
*/
  /*
connector se3HeadingConnector
  extends Avionics.Types.Pose;
  annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
end se3HeadingConnector;
*/
  /*
connector se3PositionConnector
  extends Avionics.Types.Pose;
  annotation(defaultComponentName = "Tracking Variable", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Rectangle(extent = {{-10,10},{10,-10}}, lineColor = {95,95,95}, lineThickness = 0.5),Rectangle(extent = {{-30,100},{30,-100}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}, grid = {1,1}, initialScale = 0.16), graphics = {Text(extent = {{-140,-50},{140,-88}}, lineColor = {0,0,0}, textString = "%name"),Rectangle(extent = {{-12,40},{12,-40}}, lineColor = {0,0,0}, fillColor = {192,192,192}, fillPattern = FillPattern.Solid)}));
end se3PositionConnector;
*/
  //connector bla = input SI.Angle[3]  ;
  //  connector se3AttitudeConnectorOut = input SI.Angle[3];
  //  connector se3PositionConnectorOut = input SI.Position[3];
end Interfaces;


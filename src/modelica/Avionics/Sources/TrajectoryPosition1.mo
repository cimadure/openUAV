within Avionics.Sources.stepForce;
block TrajectoryPosition1
  import SIunits = Modelica.SIunits;
  extends Avionics.Interfaces.MO(nout = 3);
  parameter Real amplitude[3] = {0.4,0.4,0.6} "Amplitude of movement";
  //  parameter SIunits.Frequency freqHz(start = 2) "Frequency of sine wave";
  parameter SIunits.Angle phase = 0 "Phase of movement";
  parameter SIunits.Time startTime = 0 "Output = offset for time < startTime" annotation(Placement(visible = true, transformation(origin = {1.10092,-31.1927}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  parameter Real offset = 0 "Offset of output signal" annotation(Placement(visible = true, transformation(origin = {2.93578,-63.1193}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  // annotation(experiment(StartTime = 0.0, StopTime = 1.0, Tolerance = 0.000001));
protected
  constant Real pi = Modelica.Constants.pi annotation(Placement(visible = true, transformation(origin = {-0.733945,31.9266}, extent = {{-12,-12},{12,12}}, rotation = 0)));
  /*  Real x;
  Real y;
  Real z;*/
equation
  signal[1] = offset + (if time < startTime then 0 else amplitude[1] * time);
  signal[2] = offset + (if time < startTime then 0 else amplitude[2] * Modelica.Math.sin(pi * (time - startTime) + phase));
  signal[3] = offset + (if time < startTime then 0 else amplitude[3] * Modelica.Math.cos(pi * (time - startTime) + phase));
  //  signal = {x,y,z};
end TrajectoryPosition1;


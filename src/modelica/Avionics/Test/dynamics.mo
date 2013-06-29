within Avionics.Test;
model dynamics
  Avionics.Bodies.dynamics dynamics1 annotation(Placement(visible = true, transformation(origin = {-21.6216,50.2703}, extent = {{-10,-10},{10,10}}, rotation = 0)));
  Avionics.Sources.CommandLaws commandlaws1(f = {0,0,82.5}) annotation(Placement(visible = true, transformation(origin = {-60,52.0868}, extent = {{-10,-10},{10,10}}, rotation = 0)));
equation
  connect(commandlaws1.Laws,dynamics1.Laws);
end dynamics;


within Avionics.Functions;
model R_derivative "Matrix Derivative (for non-linear system)"
  extends Modelica.Icons.Function;
  input Real R[3,3] "Input matrix";
  output Real dotR[3,3];
  //  input Real R[:,:] "Input matrix";
  //  output Real dotR[:,:]
protected
  Modelica.Blocks.Continuous.Derivative deriv;
  //algorithm
equation
  for i1 in 1:size(R, 1) loop
  for i2 in 1:size(R, 2) loop
  deriv.u = R[i1,i2];
  dotR[i1,i2] = deriv.y;

  end for;

  end for;
  // deriv.u := R[i1,i2];
  // dotR[i1, i2] := deriv.y;
end R_derivative;


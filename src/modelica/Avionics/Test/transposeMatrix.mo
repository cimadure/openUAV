within Avionics.Test;
model transposeMatrix
  Real R[3,3];
  output Real out[3,3];
equation
  R = [[1,2,3];[4,5,6];[7,8,9]];
  out = transpose(R);
  annotation(experiment(StartTime = 0.0, StopTime = 10.0, Tolerance = 0.000001));
end transposeMatrix;


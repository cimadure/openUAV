within Avionics.Test;
function M3MultV3_v1
  Real R[3,3] = [1,2,3;4,5,6;7,8,9];
  Real V[3,1] = [1;3;7];
  Real res[3,1];
algorithm
  res:=R * V;
end M3MultV3_v1;


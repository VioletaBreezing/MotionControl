function S = skew(phi)
S = [   0   -phi(3)  phi(2);
      phi(3)    0   -phi(1);
     -phi(2)  phi(1)    0 ];
end
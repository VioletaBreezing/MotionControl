clear;


RH = norm([32.01e-3, 55.44e-3]);
RV = norm([125.57e-3, 72.5e-3]);

GBh_inv = [1, -1/2, -1/2;
           0, -sqrt(3)/2, sqrt(3)/2;
           -RH, -RH, -RH];

GBv_inv= [1, 1, 1;
          RV, -1/2*RV, -1/2*RV;
          0, -sqrt(3)/2*RV, sqrt(3)/2*RV];

GB_h = inv(GBh_inv)
GB_v = inv(GBv_inv)
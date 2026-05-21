% DOI: 10.23919/ACC.2004.1384042
% Title: Trajectory Planning and Feedforward Design for High Performance
% Motion Systems
% Authors: Paul Lambrechts, Matthijs Boerlage, Maarten Steinbuch

configure.d = 1500;
configure.j = 80;
configure.a = 8;
configure.v = 2;
configure.p = 1;

configure.td = calculate_td(configure);
configure.tj = calculate_tj(configure);

function td = calculate_td(configure)
    d = configure.d;
    j = configure.j;
    a = configure.a;
    v = configure.v;
    x = configure.p;

    td = sqrt(sqrt(x/(8*x)));
    td = min(td, power(v/(2*d), 1/3));
    td = min(td, sqrt(a/d));
    td = min(td, j/d); 
end

function tj = calculate_tj(configure)
    d = configure.d;
    j = configure.j;
    a = configure.a;
    v = configure.v;
    x = configure.p;
    td = configure.td;

    pol3 = [1, 0, 0, 0];
    pol3(2) = 5*td;
    pol3(3) = 8*td^2;
    pol3(4) = 4*td^3 - x/(2*d*td);

    r = roots(pol3);
    tj = r(imag(r) == 0);

    pol2 = [1, 0, 0];
    pol2(2) = 3*td;
    pol2(3) = 2*td^2 - v/(d*td);
    r = roots(pol2);
    tj_v = r((imag(r) == 0 & r > 0));
    
end

function ta = calculate_ta(configure)
    d = configure.d;
    j = configure.j;
    a = configure.a;
    v = configure.v;
    x = configure.p;
end

function tv = calculate_tv(configure)
    d = configure.d;
    j = configure.j;
    a = configure.a;
    v = configure.v;
    x = configure.p;
end

function find_real_root()
end


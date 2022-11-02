function yy = quadbasis(w)
% 3 node quadratic basis set
%
%
%
x = linspace(w(1),w(end),5);
y = [1 .5 0  0 0;
     0 .5 1 .5 0;
     0  0 0 .5 1];

p1 = fit(x.',y(1,:).','poly2');
p2 = fit(x.',y(2,:).','poly2');
p3 = fit(x.',y(3,:).','poly2');

c1 = smooth(p1(x));
c2 = smooth(p2(x));
c3 = smooth(p3(x));

yy(1,:) = spline(x,c1,w);
yy(2,:) = spline(x,c2,w);
yy(3,:) = spline(x,c3,w);

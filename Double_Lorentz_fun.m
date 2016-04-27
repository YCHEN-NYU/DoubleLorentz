% single peak Lorentz profile

function y = Double_Lorentz_fun(p,x)

y = zeros(length(x),2); % allocate yout
% p(1): a11, p(2): w11, p(3): theta11, p(4): x10, p(5): a12, p(6): w2, p(7):
% theta2, p(8): x20, p(9)
y(:,1) = p(1)*(p(2)*cos(p(3))+(x-p(4))*sin(p(3)))./(p(2)^2+(x-p(4)).^2)+...
    + p(5)*(p(6)*cos(p(7))+(x-p(8))*sin(p(7)))./(p(6)^2+(x-p(8)).^2) + ...
    p(9)+p(10)*x+p(11)*x.^2;

y(:,2) = p(15)*(p(2)*cos(p(3)+pi/2)+(x-p(4))*sin(p(3)+pi/2))./(p(2)^2+(x-p(4)).^2)+...
    + p(16)*(p(6)*cos(p(7)+pi/2)+(x-p(8))*sin(p(7)+pi/2))./(p(6)^2+(x-p(8)).^2) + ...
p(12)+p(13)*x+p(14)*x.^2;

end
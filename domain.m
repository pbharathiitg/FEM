function [ke,fe] = domain(xi, eta, coord, k, Q)

% This function calculate element Stiffness Matrix

% INPUT:
% ======
% k = Thermal_conductivity
% Q = Heat generation per unit volume
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);

alpha = 1 - xi - eta ;
ddmat = [4*xi - 1, 0, 1 - 4*alpha, 4*eta, -4*eta, 4-8*xi-4*eta;
         0,   4*eta - 1,  1-4*alpha, 4*xi, 4-8*eta-4*xi, -4*xi];

J = ddmat*coord;

N = [xi*(2*xi -1), eta*(2*eta - 1), alpha*(2*alpha - 1), 4*xi*eta, 4*alpha*eta, 4*alpha*xi];

B = inv(J)*ddmat;
ke = 0.5*k*det(J)*(B'*B);

fe = N'*0.5*Q*det(J);

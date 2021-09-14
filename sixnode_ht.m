clear all

k = 30;  h = 60; Tinf = 25; Qele12 = 0; Qele34 = 1e6;qn = 2e5;


% for triangular area integration
% -------------------------------------------------------
xi1 = 1/3;  eta1 = 1/3;  w1 = -27/48;
xi2 = 0.6;  eta2 = 0.2;  w2 = 25/48;
xi3 = 0.2;  eta3 = 0.6;  w3 = 25/48;
xi4 = 0.2;  eta4 = 0.2;  w4 = 25/48;

% Gausspoints for line integration
% --------------------------------
beta1 = -0.774597;   ew1 = 5/9;
beta2 = 0.774597;    ew2 = 5/9;
beta3 = 0;           ew3 = 8/9;

% Element 1
coord = [0.0, 0.0;
         0.5, 0.0;
         0, 0.6;
         0.25, 0.0;
         0.25, 0.3;
         0, 0.3];
[ke1dg1,fe1dg1] = domain(xi1, eta1, coord, k, Qele12);
[ke1dg2,fe1dg2] = domain(xi2, eta2, coord, k, Qele12);
[ke1dg3,fe1dg3] = domain(xi3, eta3, coord, k, Qele12);
[ke1dg4,fe1dg4] = domain(xi4, eta4, coord, k, Qele12);
ke1d = ke1dg1*w1 + ke1dg2*w2 + ke1dg3*w3 + ke1dg4*w4;
fe1d = fe1dg1*w1 + fe1dg2*w2 + fe1dg3*w3 + fe1dg4*w4;

[ke1hg1,fe1hg1] = gamah(beta1, coord, h, Tinf, 1);
[ke1hg2,fe1hg2] = gamah(beta2, coord, h, Tinf, 1);
[ke1hg3,fe1hg3] = gamah(beta3, coord, h, Tinf, 1);
ke1h = ke1hg1*ew1 + ke1hg2*ew2 + ke1hg3*ew3;
fe1h = fe1hg1*ew1 + fe1hg2*ew2 + fe1hg3*ew3;

ke1 = ke1d + ke1h;

fe1 = fe1d + fe1h;


% Element 2
coord = [0.5, 0.0;
         0.5, 0.6;
         0.0, 0.6;
         0.5, 0.3;
         0.25, 0.6;
         0.25, 0.3];
     
[ke2dg1,fe2dg1] = domain(xi1, eta1, coord, k, Qele12);
[ke2dg2,fe2dg2] = domain(xi2, eta2, coord, k, Qele12);
[ke2dg3,fe2dg3] = domain(xi3, eta3, coord, k, Qele12);
[ke2dg4,fe2dg4] = domain(xi4, eta4, coord, k, Qele12);
ke2d = ke2dg1*w1 + ke2dg2*w2 + ke2dg3*w3 + ke2dg4*w4;
fe2d = fe2dg1*w1 + fe2dg2*w2 + fe2dg3*w3 + fe2dg4*w4;

[ke2hg1,fe2hg1] = gamah(beta1, coord, h, Tinf, 2);
[ke2hg2,fe2hg2] = gamah(beta2, coord, h, Tinf, 2);
[ke2hg3,fe2hg3] = gamah(beta3, coord, h, Tinf, 2);
ke2h = ke2hg1*ew1 + ke2hg2*ew2 + ke2hg3*ew3;
fe2h = fe2hg1*ew1 + fe2hg2*ew2 + fe2hg3*ew3;

ke2 = ke2d + ke2h ;
fe2 = fe2d + fe2h ;

% Element 3
coord = [0.5, 0.0;
         0.8, 0.3;
         0.5, 0.3;
         0.65, 0.15;
         0.65, 0.3;
         0.5, 0.15];
[ke3dg1,fe3dg1] = domain(xi1, eta1, coord, k, Qele34);
[ke3dg2,fe3dg2] = domain(xi2, eta2, coord, k, Qele34);
[ke3dg3,fe3dg3] = domain(xi3, eta3, coord, k, Qele34);
[ke3dg4,fe3dg4] = domain(xi4, eta4, coord, k, Qele34);
ke3d = ke3dg1*w1 + ke3dg2*w2 + ke3dg3*w3 + ke3dg4*w4;
fe3d = fe3dg1*w1 + fe3dg2*w2 + fe3dg3*w3 + fe3dg4*w4;

fe3qg1 = gamaq(beta1, coord, qn, 2);
fe3qg2 = gamaq(beta2, coord, qn, 2);
fe3qg3 = gamaq(beta3, coord, qn, 2);

fe3q = fe3qg1*ew1 + fe3qg2*ew2 + fe3qg3*ew3;

ke3 = ke3d ;

fe3 = fe3d + fe3q;

% Element 4
coord = [0.5, 0.0;
         0.8, 0.0;
         0.8, 0.3;
         0.65, 0.0;
         0.8, 0.15;
         0.65, 0.15];
[ke4dg1,fe4dg1] = domain(xi1, eta1, coord, k, Qele34);
[ke4dg2,fe4dg2] = domain(xi2, eta2, coord, k, Qele34);
[ke4dg3,fe4dg3] = domain(xi3, eta3, coord, k, Qele34);
[ke4dg4,fe4dg4] = domain(xi4, eta4, coord, k, Qele34);
ke4d = ke4dg1*w1 + ke4dg2*w2 + ke4dg3*w3 + ke4dg4*w4;
fe4d = fe4dg1*w1 + fe4dg2*w2 + fe4dg3*w3 + fe4dg4*w4;

fe4qg1 = gamaq(beta1, coord, qn, 2);
fe4qg2 = gamaq(beta2, coord, qn, 2);
fe4qg3 = gamaq(beta3, coord, qn, 2);

fe4q = fe4qg1*ew1 + fe4qg2*ew2 + fe4qg3*ew3;

[ke4hg1,fe4hg1] = gamah(beta1, coord, h, Tinf, 1);
[ke4hg2,fe4hg2] = gamah(beta2, coord, h, Tinf, 1);
[ke4hg3,fe4hg3] = gamah(beta3, coord, h, Tinf, 1);
ke4h = ke4hg1*ew1 + ke4hg2*ew2 + ke4hg3*ew3;
fe4h = fe4hg1*ew1 + fe4hg2*ew2 + fe4hg3*ew3;

ke4 = ke4d + ke4h ;
fe4 = fe4d + fe4h + fe4q;

% Assembly
K = zeros(16,16);
F = zeros(16,1);

K([1:2,7,8,14,13],[1:2,7,8,14,13]) = ke1(1:6,1:6);
K([2,6,7,5,12,14],[2,6,7,5,12,14]) = K([2,6,7,5,12,14],[2,6,7,5,12,14]) + ke2(1:6,1:6);
K([2,4,5,16,11,15],[2,4,5,16,11,15]) = K([2,4,5,16,11,15],[2,4,5,16,11,15]) + ke3(1:6,1:6);
K([2,3,4,9,10,16],[2,3,4,9,10,16]) = K([2,3,4,9,10,16],[2,3,4,9,10,16]) + ke4(1:6,1:6);

F([1:2,7,8,14,13]) = fe1(1:6);
F([2,6,7,5,12,14]) = F([2,6,7,5,12,14]) + fe2(1:6);
F([2,4,5,16,11,15]) = F([2,4,5,16,11,15]) + fe3(1:6);
F([2,3,4,9,10,16]) = F([2,3,4,9,10,16]) + fe4(1:6);



% Imposition of B.C.
Kreduce = K([2:6,8:12,14:16],[2:6,8:12,14:16]);
Freduce = F([2:6,8:12,14:16]) - (K(1,[2:6,8:12,14:16])*200 + K(7,[2:6,8:12,14:16])*200 + K(13,[2:6,8:12,14:16])*200)';

% Finding Solution
ureduce = inv(Kreduce)*Freduce;
un = [200;ureduce(1:5);200;ureduce(6:10);200;ureduce(11:13)] ;




% % ==============================================
% %% Printing Intermediate Result to The Output File
% % ------------------------------------------------
fid=fopen('Steps_sixnode_traingular_element','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'k = %12.4e, h = %12.4e, Tinf = %12.4e, Qele12 = %12.4e, Qele34 = %12.4e, qn = %12.4e\n\n',k,h,Tinf,Qele12,Qele34,qn);
fprintf(fid,'Ke1 \n');
fprintf(fid,'----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',ke1(i,1:6));
end
fprintf(fid,'Ke2 \n');
fprintf(fid,'----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',ke2(i,1:6));
end
fprintf(fid,'Ke3 \n');
fprintf(fid,'----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',ke3(i,1:6));
end
fprintf(fid,'Ke4 \n');
fprintf(fid,'----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',ke4(i,1:6));
end



fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'fe1 \n');
fprintf(fid,'-----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\n',fe1(i));
end
fprintf(fid,'\nfe2 \n');
fprintf(fid,'-----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\n',fe2(i));
end
fprintf(fid,'\nfe3 \n');
fprintf(fid,'-----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\n',fe3(i));
end
fprintf(fid,'\nfe4 \n');
fprintf(fid,'-----\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\n',fe4(i));
end

fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:16
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\n\n',K(i,1:16));
end
fprintf(fid,'\n\n\n');



fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:16
   fprintf(fid,'%12.4e\n',F(i));
end

fprintf(fid,'\n\nImposition of Boundary Condition\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K u = F\n');
fprintf(fid,'--------\n');
for i = 1:16
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:16),F(i));
end

fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:13
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',Kreduce(i,1:13),Freduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'Tn\n');
fprintf(fid,'--\n');
for i = 1:16
   fprintf(fid,'T%d = %12.4e\n\n',i,un(i));
end


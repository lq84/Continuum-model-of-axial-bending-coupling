function [B, J_det] = BMatrixQ(xi,eta,coord)
%This sets up the matrix that helps transform the displacement into
%the strain and then the stress.
%The inputs are: 
%xi and eta, scalar values.
%coord a 8x2 matrix containing the (x,y) locations of the element's nodes.
%The outputs are the B matrix, a 6x24 matrix, and J_det the determinant of
%the Jacobian.
N = [1/4*(1-xi)*(1-eta)*(-1-xi-eta),...
     1/4*(1+xi)*(1-eta)*(-1+xi-eta),...
     1/4*(1+xi)*(1+eta)*(-1+xi+eta),...
     1/4*(1-xi)*(1+eta)*(-1-xi+eta),...
     1/2*(1-xi^2)*(1-eta),...
     1/2*(1-eta^2)*(1+xi),...
     1/2*(1-xi^2)*(1+eta),...
     1/2*(1-eta^2)*(1-xi)];
B_hat = [-1/4*(1-eta)*(-xi-eta-1)-1/4*(1-xi)*(1-eta),...
         -1/4*(1-xi )*(-xi-eta-1)-1/4*(1-xi)*(1-eta);
          1/4*(1-eta)*(xi-eta-1 )+1/4*(1+xi)*(1-eta),...
         -1/4*(1+xi)*(xi-eta-1)-1/4*(1+xi)*(1-eta);
          1/4*(1+eta)*(xi+eta-1)+1/4*(1+xi)*(1+eta),...
          1/4*(1+xi)*(xi+eta-1)+1/4*(1+xi)*(1+eta);
         -1/4*(1+eta)*(-xi+eta-1)-1/4*(1-xi)*(1+eta),...
          1/4*(1-xi )*(-xi+eta-1)+1/4*(1-xi)*(1+eta);
          -xi*(1-eta),-1/2+1/2*xi^2;
          1/2-1/2*eta^2,-eta*(1+xi);
          -xi*(1+eta),1/2-1/2*xi^2;
          -1/2+1/2*eta^2,-eta*(1-xi)]';

% Calculate the Jacobian
J=B_hat*coord;

% Calculate the determinate

J_det = det(J);
% J_inverse = inv(J);

% Form the B matrix
% B_2by4 = J_inverse*B_hat;
%dbstop if warning
B_2by8 = J\B_hat;

% Form B, which is an 6 by 24
% The form of this matrix was derived from
%Quadrilateral isoparametric finite elements for plane elastic
%Cosserat bodies
%My notes on this are written up in a Maple document of the same
%name.

%% Micropolar theory
B = sparse(6,8*3);
%corresponds to epsilon_xx = ux,x
B(1,1:3:24) = B_2by8(1,:);

%corresponds to epsilon_yy = uy,y
B(2,2:3:24) = B_2by8(2,:);

%Corresponds to epsilon_xy = uy,x - phi
B(3,2:3:24) = B_2by8(1,:);
B(3,3:3:24) = -N;

%Corresponds to epsilon_yx = ux,y + phi
B(4,1:3:24) = B_2by8(2,:);
B(4,3:3:24)      = N;
B(5:6,3:3:24)    = B_2by8; %Corresponds to k13, k23.



end
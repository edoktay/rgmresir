% GMRES with double/quadruple precision for use with GCRODR and GMRESIR
%
% Solves A z = r for z, then returns x + z
% Generates Arnoldi relation A V(:,1:m) = V(:,1:m+1) H
%
% INPUT:  A      N-by-N matrix
%         X      current solution vector
%         R      N-by-1 preconditioned residual vector
%         M      number of GMRES iterations to perform
%         M1     left preconditioner for A
%         M2     right preconditioner for A
%         tol    specifies the tolerance of the method
% OUTPUT: x      updated solution vector
%         R      preconditioned residual vector
%         V      N-by-M+1 matrix containing orthogonal basis for Krylov subspace 
%         H      M+1-by-M upper Hessenburg reduction of matrix operator
%         K      number of GMRES iterations actually performed
%         RESVEC vector containing norm of residual at each iteration of GMRES
%
% Adapted from the codes developed in Parks, Michael L., et al. "Recycling
% Krylov subspaces for sequences of linear systems.", SIAM Journal on
% Scientific Computing 28.5 (2006): 1651-1674.
%
% The library and associated functions of GCRO-DR are available at 
% https://www.sandia.gov/-mlparks/software/.
function [x,r,V,H,k,resvec] = gmres1_dq(A,x,r,m,M1,M2,tol)

%Ensure double working precision
A = double(A);
x = double(x);
M1 = double(M1);
M2 = double(M2);
r = double(r);

% Initialize V

size(r);
V(:,1) = r / norm(r);
for k = 1:m
   % Find w using preconditioning if available.
   w = double(mp(M2,34)\(mp(M1,34)\(mp(A,34)*mp(V(:,k),34))));

   % Create next column of V and H
   for j = 1:k
      H(j,k) = V(:,j)' * w;
      w = w - H(j,k) * V(:,j);
   end

   H(k+1,k) = norm(w);
   V(:,k+1) = w / H(k+1,k);

   % Initialize right hand side of least-squares system
   rhs = zeros(k+1,1);
   rhs(1) = norm(r);

   % Solve least squares system; Calculate residual norm
   y = H \ rhs;
   res = rhs - H * y;
   resvec(k) = norm(res);
   if resvec(k) <= tol
      % Calculate solution and residual and return
      x = x + V(:,1:k) * y;
      r = V * res;
      return
   end
end

% Calculate solution and residual.
x = x + V(:,1:m) * y;
r = V * res;
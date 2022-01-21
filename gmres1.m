% GMRES for use with GCRODR and GMRESIR
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
% The library and associated functions of GCRO-DR are available at 
% https://www.sandia.gov/-mlparks/software/.
%
function [x,r,V,H,k,resvec] = gmres1(A,x,r,m,M1,M2,tol)

% Initialize V

size(r);
V(:,1) = r / norm(r);

for k = 1:m
   % Find w using preconditioning if available.
   w = single(double(M2)\(double(M1)\(double(A)*double(V(:,k)))));

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
   if resvec(k) < tol
      % Calculate solution and residual and return
      x = x + V(:,1:k) * y;
      r = V * res;
      return
   end
end

% Calculate solution and residual.
x = x + V(:,1:m) * y;
r = V * res;
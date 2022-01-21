% GCRODR   GCRO with Deflated Restarting with single/double precision, 
%     the maximum subspace dimension m and the number k of approximate 
%     eigenvectors kept from one cycle to the next.
%
% Adapted from the codes developed in Parks, Michael L., et al. "Recycling
% Krylov subspaces for sequences of linear systems.", SIAM Journal on
% Scientific Computing 28.5 (2006): 1651-1674.
%
% The library and associated functions of GCRO-DR are available at 
% https://www.sandia.gov/-mlparks/software/.
%
% Note: Run 'clear all' between calls to GCRO-DR to ensure no subspace is recycled

function [x,resvec,r,nmv,relres] = gcrodr(A,b,m,k,x0,tol,M1,M2,reuse_name)

% Initialize optional variables.
if(nargin < 5 | isempty(x0))
   x0 = zeros(size(b));
end
if(nargin < 6 | isempty(tol))
   tol = 1e-6;
end
if(nargin < 7 | isempty(M1))
   existM1 = 0;
   M1 = [];
else
   existM1 = 1;
end
if(nargin < 8 | isempty(M2))
   existM2 = 0;
   M2 = [];
else
   existM2 = 1;
end
if(nargin < 9 | isempty(reuse_name))
   reuse_name = 'default';
end

% initialize solution vector
x = zeros(size(x0));

% initialize matvec count
nmv = 1;

% Calculate initial preconditioned residual.
r = b - A*x0;
r = single(double(M2)\(double(M1)\double(r)));

% Calculate initial preconditioned residual norm.
resvec = zeros(2,1);
resvec(1) = norm(r);
disp(sprintf('||r|| = %e\t\tnmv = %d',resvec(nmv),nmv-1));

% Precondition rhs if available.
if(existM1)
   bnorm = norm(single(double(M2)\(double(M1)\double(b))));
else
   bnorm = norm(b);
end

persistent U_persist;
%%%%%%%%%%%%%%%%%%%  Initialize U (Recycled subspace) %%%%%%%%%%%%%%%%%%%
if(isfield(U_persist,reuse_name))
   % Initialize U with information recycled from previous call to solver.
   eval(sprintf('U = U_persist.%s;',reuse_name));

   % C = A * U (with preconditioning)
   % We can frequently represent A_new = A_old + deltaA. Note that
   % A_new * U = A_old*U + delta_A*U
   %           = C_old   + delta_A*U
   % where we already have C_old. Computing deltaA*U is generally much less 
   % expensive than computing A_new*U, so we do not record these (k) matvecs
   C = single(double(M2)\(double(M1)\(double(A)*double(U))));
 
   % Orthonormalize C and adjust U accordingly so that C = A*U
   [C,R] = qr(C,0);
   U = single(double(U) / double(R));

   % Set residual norm to be the initial residual norm for this case.
   x = x + U*C'*r;
   r = r - C*C'*r;
   resvec(1) = norm(r);
else
   % I have no subspace to recycle from a previous call to the solver
   % So, do one cycle of GMRES to "prime the pump"

   % Perform m GMRES iterations to produce Krylov subspace of dimension m
   [x,r,V,H,p,resvec_inner] = gmres1(A,x,r,m,M1,M2,tol*bnorm);

   % Record residual norms and increment matvec count
   resvec(2:p+1) = resvec_inner;
   nmv = nmv + p;
   disp(sprintf('||r|| = %e\t\tnmv = %d',resvec(nmv),nmv-1));

   % Find the k smallest harmonic Ritz vectors.
   % Check to be sure GMRES took more than k iterations. Else, I can't compute
   % k harmonic Ritz vectors
   if k < p
      P = getHarmVecs1(p,k,H);

     % Form U (the subspace to recycle)
      U = V(:,1:p) * P;
  
      % Form orthonormalized C and adjust U accordingly so that C = A*U
      [C,R] = qr(H*P,0);
      C = V * C;
      U = single(double(U) / double(R));
   end

   % If p < m, early convergence of GMRES
   if p < m
      % Assign all (unassigned) outgoing values and return
      x = x0 + x;  
      nmv = nmv - 1;
      relres = resvec(p+1) / bnorm;
      % Save information to be carried to next call to solver.
      if k < p
         eval(sprintf('U_persist.%s = U;',reuse_name));
      end
      return
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%  Main Body of Solver  %%%%%%%%%%%%%%%%%%%%%%%%%

while(resvec(nmv) / bnorm > tol) 
   % Do m-k steps of Arnoldi
   [V,H,B,p,resvec_inner] = gmres2(A,r,m-k,M1,M2,C,tol*bnorm);
   resvec(nmv+1:nmv+p) = resvec_inner;
   nmv = nmv + p;
   disp(sprintf('||r|| = %e\t\tnmv = %d',resvec(nmv),nmv-1));

   % Rescale U; Store inverses of the norms of columns of U in diagonal matrix D
   for i = 1:k
      d(i) = norm(U(:,i));
      U(:,i) = U(:,i) / d(i);
   end
   D = diag(1 ./ d);

   % Form large H
   H2 = zeros(p+k+1,p+k);
   H2(1:k,1:k) = D;
   H2(1:k,k+1:p+k) = B;
   H2(k+1:p+k+1,k+1:p+k) = H;

   % Calculate solution update
   rhs = [C V]' * r;
   y = H2 \ rhs;
   x = x + [U V(:,1:p)] * y;

   % Calculate new residual
   r = r - [C V] * (H2 * y);

   % If p < m-k, early convergence of GMRES
   if p < m-k
      x = x0 + x;
      relres = resvec(nmv) / bnorm;
      nmv = nmv - 1;
      % Save information to be carried to next call to solver.
      eval(sprintf('U_persist.%s = U;',reuse_name));
      return
   end

   % Calculate Harmonic Ritz vectors.
   P = getHarmVecs2(p+k,k,H2,V,U,C);

   % Form new U and C.
   U = [U V(:,1:p)] * P;

   % Form orthonormalized C and adjust U accordingly so that C = A*U
   [Q,R] = qr(H2*P,0);
   C = [C V] * Q;
   U = U / R; 

end

% Save information to be carried to next call to solver.
eval(sprintf('U_persist.%s = U;',reuse_name));

% Calculate final solution and residual.
 x = x0 + x;

% Calculate relative residual.
relres = resvec(nmv) / bnorm;

% Correct matvec count
nmv = nmv - 1;
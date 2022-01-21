function [x, error, its, flag] = gmres_dq( A, x, b, L, U, restrt, tol)
%GMRES_DQ   Left-preconditioned GMRES in double/quadruple precision
%   Solves Ax=b by solving the preconditioned linear system (LU)^{-1}Ax=(LU)^{-1}b
%   using the Generalized Minimal residual ( GMRES ) method.
%   Currently uses (preconditioned) relative residual norm to check for convergence 
%   (same as Matlab GMRES)
%   Single precision used throughout, except in applying (U\L\A) to a vector 
%   which is done in double precision
%
%   input   A        REAL nonsymmetric positive definite matrix
%           x        REAL initial guess vector
%           b        REAL right hand side vector
%           L        REAL L factor of lu(A)
%           U        REAL U factor of lu(A)
%           restrt   INTEGER number of iterations between restarts
%           tol      REAL error tolerance
%
%   output  x        REAL solution vector
%           error    REAL error norm
%           iter     INTEGER number of (inner) iterations performed
%           flag     INTEGER: 0 = solution found to tolerance
%                             1 = no convergence given max_it

flag = 0;
its = 0;

%Ensure single working precision
A = double(A);
b = double(b);
x = double(x);

%Cast half precision L and U factors as doubles
L = double(L);
U = double(U);
n = size(A,1);

rtmp = b-A*x;
r = mp(double(L),34)\mp(double(rtmp),34);
r = mp(double(U),34)\r;
r = double(r);
resvec(1) = norm(r);
nmv = 1;
disp(sprintf('||r|| = %e\t\tnmv = %d',resvec(nmv),nmv-1));
M1 = L; M2 = U;
init = 1;

% Precondition rhs
bnorm = norm(double(mp(M2,34)\(mp(M1,34)\mp(b,34))));
[x,r,V,H,p,resvec1] = gmres1_dq(A,x,r,restrt,M1,M2,tol*bnorm);
resvec(2:p+1) = resvec1;
nmv = nmv+p;
disp(sprintf('||r|| = %e\t\tnmv = %d',resvec(nmv),nmv-1));

while(resvec(nmv)/bnorm > tol)
    init = init+1;
    [x,r,V,H,p,resvecp] = gmres1_dq(A,x,r,restrt,M1,M2,tol*bnorm);
    resvec(nmv+1:nmv+p) = resvecp;
    nmv = nmv + p;
    disp(sprintf('||r|| = %e\t\tnmv = %d',resvec(nmv),nmv-1));
    if p<restrt
        error = resvec;
        nmv = nmv - 1;
        its = nmv;
        return
    end
end

error = resvec;
nmv = nmv - 1;
its = nmv;

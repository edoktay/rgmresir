% Run tests for matrices from mode 3 randsvd matrices
% Note : run examples one by one, clear all the workspace before running
% rgmresir again.

close all
clear all
clc
warning off

n = 100;
maxit = 1000;
condnums = [1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12];
uf = 0; u =2; ur = 4;
mm =40; kk =20;
restrt = mm;

i = 6;
fprintf('\nRunning HDQ test for mode 3 random matrix with condition number 1e%s\n',num2str(log10(condnums(i))));

rng(1);
A = gallery('randsvd',n,condnums(i),3);
b = randn(n,1);

snbase = strcat('figs/mode3_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(restrt),num2str(uf),num2str(u),num2str(ur)),restrt);
rgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_GCRODR_',num2str(mm),'_',num2str(kk),'_',num2str(uf),num2str(u),num2str(ur)),mm,kk);

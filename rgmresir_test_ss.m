% Run tests for matrices from SuiteSparse
% Note: requires ssget MATLAB interface
% Note 2: run examples one by one, clear all the workspace before running
% rgmresir again.

clc
close all
clear all
warning off

matids = [907, 1641, 293, 906, 462, 464, 1199, 253, 2338,168,172,232,262,293,318,463,1326,1487,1943];
maxit = 1000;
mm = 40; kk = 20;
restrt = mm;

i = 1;
Problem = ssget(matids(i));
A = full(Problem.A);
b = ones(size(A,2),1);

uf = 0; u = 2; ur = 4;

namearr = split(Problem.name,'/');
name = namearr{2};
fprintf('Running HDQ tests for matrix %s\n', name);

snbase = strcat('figs/',name,'_');
gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(restrt),num2str(uf),num2str(u),num2str(ur)),restrt);
rgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'RGMRESIR_',num2str(mm),'_',num2str(kk),'_',num2str(uf),num2str(u),num2str(ur)),mm,kk);



       


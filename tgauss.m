function [A,B,tc]=tgauss(n)
%eval(EXPRESSION) evaluates the MATLAB code in EXPRESSION. Specify
%    EXPRESSION as a character vector or string scalar.
eval(strcat('tgauss',num2str(n)));

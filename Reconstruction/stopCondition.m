function tolTyp=stopCondition(tolType,solv)

%STOPCONDITION Decides on the stop condition of CGSOLVER and ADMMSOLVER
%   TOLTYP=STOPCONDITION(TOLTYPE,SOLV)
%   * TOLTYPE is a string denoting the stopping condition
%   * SOLV is a string denoting the solver
%   * TOLTYP is an integer describing the stopping condition
%

tolTyp=-1;
if strcmp(tolType,'NormwiseBackward2Error');tolTyp=0;%|r_k|/(|A||x_k|+|b|) from P Tichy, "On error estimation in the conjugate gradient method: normwise backward error," Proc of ALGORITMY, 323-332, 2016.
elseif strcmp(tolType,'RelativeResidual2Error');tolTyp=1;%|r_k|/|b|
elseif strcmp(tolType,'DataUpdate2Error');tolTyp=2;%|x_k-x_{k-1}|/|b|
elseif strcmp(tolType,'RelativeResidualIError');tolTyp=3;
elseif strcmp(tolType,'DataUpdateIError');tolTyp=4;
elseif strcmp(tolType,'Convergent2Variables');tolTyp=5;
elseif strcmp(tolType,'Energy');tolTyp=6;
elseif strcmp(tolType,'None');tolTyp=7;
end

assert(~((strcmp(solv,'CG') && ismember(tolTyp,[-1 5])) || (strcmp(solv,'ADMM') && ismember(tolTyp,[-1 0 1 3])) || (strcmp(solv,'IRWLS') && ismember(tolTyp,[-1 0 1 2 3 4 5]))),'The %s stopping criterion is not implemented for %s solver',tolType,solv);

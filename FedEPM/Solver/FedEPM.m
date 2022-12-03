function out = FedEPM(di,n,A,b,k0,pars)
% This solver solves logistic regression problem in the following form:
%
%         min_{x_i,x\in\R^n}  sum_{i=1}^m f_i(x_i;(A_i,b_i))  
%            s.t.             x_i=x, i=1,2,...,m
%
% where 
%      f_i(x;(A_i,b_i)) = (mu/2)*||x||^2 
%                       + sum_{j=1}^{d_i} (log(1+exp(<x,(a^i)_j>))-(b^i)_j*<x,(a^i)_j> )
%      (A_i,b_i) is the data for node/client i
%      A_i = [(a^i)_1, (a^i)_2, ..., (a^i)_{d_i}]^T \in\R^{d_i-by-n} 
%      b_i = [(b^i)_1, (b^i)_2, ..., (b^i)_{d_i}]^T \in\R^{d_i-by-1} 
% =========================================================================
% Inputs:
%   di      : 1-by-m row vector, di = (d_1, d_2, ..., d_m)        (REQUIRED)
%             d_i is the number of rows of A_i
%             Let d = d_1+d_2+...+d_m
%   n       : Dimension of solution x                             (REQUIRED)
%   A       : A=[A_1; A_2; ...; A_m]\in\R^{d-by-n}                (REQUIRED)
%   b       : b=[b_1; b_2; ...; b_m]\in\R^{d-by-1}                (REQUIRED)
%   k0      : A positive integer controlling communication rounds (REQUIRED)
%             The larger k0 is the fewer communication rounds are 
%   pars  :   All parameters are OPTIONAL                                  
%             pars.rho   -- Rate of devices selected to train, (default,0.5)
%                        -- The larger the more devices selected
%             pars.noise -- Degree of privacy preserving, (default,0.1)   
%                           The smaller the stronger privacy maintained
%             pars.tol   -- Tolerance of the halting condition (default,1e-9*n)
%             pars.maxit -- Maximum number of iterations (default,500*k0) 
% =========================================================================
% Outputs:
%     out.sol:      The solution x
%     out.obj:      Objective function value at out.sol
%     out.time:     CPU time
%     out.iter:     Number of iterations 
%     out.cr:       Number of communication rounds
%     out.snr:      Signal-to-noise ratio of out.sol
% =========================================================================
% Written by Shenglong Zhou on 23/08/2022 based on the algorithm proposed in
%     Shenglong Zhou,  Geoffrey Ye Li,
%     Exact Penalty Method for Federated Learning,
%     arXiv:2110.15318, 2022    	
% Send your comments and suggestions to <<< slzhou2021@163.com >>>                                  
% WARNING: Accuracy may not be guaranteed!!!!!  
% =========================================================================
warning off; rng('shuffle'); 


if  nargin < 5
    disp(' No enough inputs. No problems will be solverd!'); return;
elseif nargin<6 
    pars   = [];  
end

[m,tol,alpha,rho,maxit,noise] = set_parameters(di,n,k0,pars); 
[fun,lam,eta,r0,fAvg] = def_funs(m,n,A,b,di,pars);

Fnorm      = @(x)norm(x,'fro')^2;
snr        = zeros(m,1);
muj0       = r0.*ones(m,1);
objx       = zeros(1,maxit);
obj        = zeros(1,ceil(maxit/k0));
err        = zeros(1,ceil(maxit/k0));
X          = zeros(n,m);   
Z          = X;
for j      = 1:m
    Z(:,j) = X(:,j) + randl(n,1)/muj0(j); 
end             
m0         = ceil(rho*m);
M          = 1:m;
x          = zeros(n,1);
 
ct         = 0;
fprintf(' Start to run the solver -- FedEPM \n');
fprintf(' -----------------------------------------------------------\n');
fprintf('                          Iter        f(y)          Time  \n');  
fprintf(' -----------------------------------------------------------\n');
t0         = tic; 

% main body ------------------------------------------------ 
for iter = 0 : maxit
       
    ct0 = tic;
    if  mod(iter, k0)==0
        x       = fAvg(Z);  
        [fx,gY] = fun(x);  
        M       = randperm(m);
        M       = sort(M(1:m0)); 
        s       = iter/k0+1; 
        obj(s)  = fx;
        err(s)  = Fnorm(sum(gY,2));
    end
     
    objx(iter+1) = fx;   
    if mod(iter, k0)==0    
        fprintf(' Communication at iter = %4d      %9.4f      %6.3fsec\n',iter, fx, toc(t0));  
        stop1 = s > 3 && var(obj(s-3:s))<tol/(1+fx);
        stop2 = s > 3 && fx > max(obj(s-3:s-1)) && err(s) > max(err(s-3:s-1));
        stop3 = err(s)<1e-6;
        if stop1 || stop2  || stop3
          % fprintf(' Communication at iter = %4d      %9.4f      %6.3fsec\n',iter, fx, toc(t0));  
           break;  
        end    
    end 
    
    ct = ct + toc(ct0);
     
    for j = 1:m
        if  ismember(j,M)     
            gap        = X(:,j)-x;
            muj        = muj0(j) *(1+1e-8*Fnorm(gap))*alpha^(iter+1);  
            mud        = muj /di(j);
            tmp        = mud*gap-gY(:,j);
            X(:,j)     = x + sign(tmp).*max(abs(tmp)-lam,0)/(mud + eta);   
            if  mod(iter+1, k0)==0   
                v      = 2*norm(gY(:,j),1)/muj/noise; 
                Z(:,j) = X(:,j) + randl(n,1) *v;   
                snr(j) = log10(norm(X(:,j))/norm( Z(:,j) - X(:,j)));
            end  
        end
    end         
end

out.sol    = x;
out.objx   = obj(1:s);
out.obj    = fx;
out.snr    = min(snr);
out.acc    = 1-nnz(b-max(0,sign(A*x)))/length(b); 
out.iter   = iter+1;
out.time   = toc(t0);
out.cr     = ceil(iter/k0);
%out.ct     = (out.time-ct)/out.cr;
fprintf(' -----------------------------------------------------------\n');

end

%--------------------------------------------------------------------------
function [m,tol,alpha,rho,maxit,noise] = set_parameters(di,n,k0,pars) 
    m       = length(di);
    maxit   = 500*k0;
    tol     = 1e-9*n;  
    rho     = 0.5;
    alpha   = 1.001;  
    noise   = 0.2;
    if isfield(pars,'noise');   noise   = pars.noise;   end
    if isfield(pars,'tol');     tol     = pars.tol;     end
    if isfield(pars,'alpha');   alpha   = pars.alpha;   end
    if isfield(pars,'maxit');   maxit   = pars.maxit;   end
    if isfield(pars,'rho');     rho     = pars.rho;     end
end

%--------------------------------------------------------------------------
function [fun,lam,eta,r0,fAvg] = def_funs(m,n,A,b,di,pars)
    I       = zeros(m+1,1);
    I(1)    = 0;
    for j   = 1 : m, I(j+1) = I(j)+di(j); end
    Ai      = cell(1,m);
    bi      = cell(1,m);
    for j   = 1 : m 
        indj   = I(j)+1:I(j+1);
        Ai{j}  = A(indj,:);  
        bi{j}  = b(indj);
    end 

    eta  = 2e-5*0.6*(m/100+0.5);
    lam  = eta /2;
    r0   = 0.05;
    fun  = @(x)funcLogist(x,Ai,bi,m,n,1./di);
    if isfield(pars,'r0');  r0  = pars.r0;  end
    if isfield(pars,'lam'); lam = pars.lam; end
    if isfield(pars,'eta'); eta = pars.eta; end
    
    favg     = @(var)avg(var,lam,eta,m);  
    fAvg     = @(var)Avg(var,n,favg);

end

%--------------------------------------------------------------------------
function  [objX,gradX]  = funcLogist(x,Ai,bi,m,n,w) 
       
    objX   = 0; 
    gradX  = zeros(n,m);
    for i  = 1:m
        Ax   = Ai{i}*x;  
        eAx  = 1 + exp(Ax);
        objX = objX + w(i)* ( sum( log(eAx)-bi{i}.*Ax ) + 5e-4*norm(x,'fro')^2 ); 
        if nargout   == 2 
           gradX(:,i) =  w(i)*( ((1-bi{i}-1./eAx)'*Ai{i})'+1e-3*x);
        end
    end
end
% 
% %--------------------------------------------------------------------------
% function  [objX,gradX]  = funcLinear(x,Ai,bi,m,n,w) 
%      
%     objX     = 0; 
%     gradX    = zeros(n,m);
%     for j    = 1:m  
%         tmp  = Ai{j}*x-bi{j};
%         objX = objX  + norm( tmp )^2*w(j); 
%         if nargout   == 2
%            gradX(:,j) =  w(j)*(tmp'* Ai{j} )';
%         end
%     end
%     objX = objX/2;
% end

%--------------------------------------------------------------------------
function x = randl(m, n)
% function: x  = randl(m, n)
% m - number of matrix rows
% n - number of matrix columns
% x - matrix with Laplacian distributed random numbers 
%     with mean mu = 0 and std sigma = 1 (columnwise)
% generation of two i.i.d. sequences
u1 = rand(m, n);
u2 = rand(m, n);
% generation of a matrix with Laplacian
% distributed random numbers (columwise)
x = log(u1./u2);
x = bsxfun(@minus, x, mean(x));
x = bsxfun(@rdivide, x, std(x));
end

%--------------------------------------------------------------------------
function  mxt   = avg(x,lam,eta,m)  
          w     = sort(x,'descend');
          tmp   = lam/eta;
          mx    = mean(x)-tmp;
          for t = 1: m 
              mxt     = mx + 2*tmp*t/m;
              if t < m && w(t) >  mxt &&  mxt>w(t+1) 
                  break; 
              end
          end
          if t == m
             fun     = @(v)sum(lam*abs(v-w)+(eta/2)*(v-w).^2);
             [~,ind] = min(arrayfun(fun, w));
             mxt     = w(ind(1));
          end
end

%--------------------------------------------------------------------------
function  avx   = Avg(X,n,avg)  
          avx   = zeros(n,1);
          for i = 1:n
              avx(i) = avg(X(i,:));  
          end
 
end

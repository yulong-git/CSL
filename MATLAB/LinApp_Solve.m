function [PP, QQ, UU, RR, SS, VVV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WWW,TT,NN,Z0,Sylv)
% Written by Kerk Phillips, November 2013
% Version 1.1 updated April 1, 2014
% This code takes Uhlig's original code and puts it in the form of a
% function.  This version outputs the policy function coefficients: PP,
% QQ and UU for X, and RR, SS and VV for Y.
%
% This function take the follwoing as inputs:
%  The matrices of deriviatives: FF - TT.
%  The autoregression coefficient matrix NN from the law of motion for Z.
%  Z0 is the Z-point about which the linearization is taken.  For
%   linearizing about the steady state this is Zbar and normally Zbar = 0.
%   QQ if true.
%  Sylv is an indicator variable telling the program to use the built-in
%   function sylvester() to solve for QQ and SS, if possible.  Default is
%   to use Sylv=1.
% Source: R. Evans and K. Phillips (2014) "Linearization about the Current
% State: A Computational Method for Approximating Nonlinear Policy 
% Functions during Simulation," mimeo, Brigham Young University Department
% of Economics.
%
% Copyright: K. Phillips.  Feel free to copy, modify and use at your own 
% risk.  However, you are not allowed to sell this software or otherwise 
% impinge on its free distribution.

% Use Sylvester equation solution for QQ and SS by default
if (~exist('Sylv', 'var'))
    Sylv = 1;
end

% Find the values of nx, ny and nz from the input matrices
[nx,~] = size(FF);
[ny,~] = size(AA);
[nz,~] = size(NN);

message = '                                                                       ';
warnings = [];
DISPLAY_IMMEDIATELY = 0;
TOL = .000001; % Roots smaller than TOL are regarded as zero.
               % Complex numbers with distance less than TOL 
               % are regarded as equal.
if exist('MANUAL_ROOTS')~=1,
   MANUAL_ROOTS = 0; % = 1, if you want to choose your own
                               % roots, otherwise set = 0. 
                               % See SOLVE.M for instructions.
end;
if exist('IGNORE_VV_SING')~=1,
   IGNORE_VV_SING = 1; % =1: Ignores, if VV is singular.  
                       % Sometimes useful for sunspots.  
                       % Cross your fingers...  
end;
DISPLAY_ROOTS = 0;  % Set = 1, if you want to see the roots.

% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% SOLVE.M solves for the decision rules in a linear system,
% which is assumed to be of the form
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) +...
%   LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,
% where it is assumed that x(t) is the endogenous state vector,
% y(t) the other endogenous variables and z(t) the exogenous state
% vector.  It is assumed that the row dimension of AA is at least as large
% as the dimensionality of the endogenous state vector x(t).  
% The program solves for the equilibrium law of motion
% x(t) = PP x(t-1) + QQ z(t)
% y(t) = RR x(t-1) + SS z(t).
% To use this program, define the matrices AA, BB, .., NN.
% SOLVE.M then calculates PP, QQ, RR and SS.  It also calculates
% WW with the property [x(t)',y(t)',z(t)']=WW [x(t)',z(t)'].
% 
% A few additional variables are used
% overwriting variables with the same names
% that might have been used before.  They are:
% sumimag, sumabs, message, warnings, CC_plus, CC_0, Psi_mat, Gamma_mat,
% Theta_mat, Xi_mat, Delta_mat, Xi_eigvec, Xi_eigval, Xi_sortabs, 
% Xi_sortindex, Xi_sortvec, Xi_sortval, Xi_select, drop_index, Omega_mat,
% Lambda_mat, PP_imag, VV, LLNN_plus_MM, QQSS_vec
% 
% Source: H. Uhlig (1995) "A Toolkit for Solving Nonlinear Dynamic
% Stochastic Models Easily," Discussion Paper, Institute for
% Empirical Macroeconomis, Federal Reserve Bank of Minneapolis #101 or
% Tilburg University, CentER DP 9597.
%
% This update includes the suggestion by Andrew Atkeson to use generalized
% eigenvalues to perform the computations.  IN PARTICULAR, THE DEFINITION
% OF PSI_MAT, GAMMA_MAT and THETA_MAT have been changed!
%
% You can also select roots manually.  For the manual selection procedure,
% see the instructions upon inspecting this file (filename: solve.m)
%
% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  MANUAL SELECTION OF ROOTS PROCEDURE.  INSTRUCTIONS:     %
      % For manual selection, set MANUAL_ROOTS = 1 and %
      % define Xi_manual somewhere earlier in your calculations. %
      % Xi_manual should be a vector of length m_states with     %
      % distinct integer entries between 1 and 2*m_states. The   %
      % program then uses the roots Xi_sortval(Xi_manual) and    %
      % the corresponding eigenvectors Xi_sortvec(Xi_manual).    %
      % Thus, to choose the desired roots, run the program once  %
      % with automatic root selection, take a look at Xi_sortval,%
      % and Xi_sortvec and write down the indices of the desired %  
      % roots.                                                   %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
k_exog = min(size(NN));
sumimag = sum(sum(abs(imag(AA))))+sum(sum(abs(imag(BB))))...
    +sum(sum(abs(imag(CC))));
sumabs  = sum(sum(abs(AA)))   +sum(sum(abs(BB)))   +sum(sum(abs(CC)));
if sumimag / sumabs > .000001,
  message = ...
['SOLVE.M: I hate to point this out to you, but some of your matrices    '  
'         contain complex numbers, which does not make much sense. You  '
'         should check your steady state parameters and calculations.   '
'         I will proceed anyhow, but you will probably get nonsense.    '];
  if DISPLAY_IMMEDIATELY, disp(message); end;
  warnings = [warnings;message];
end;
if rank(CC)<n_endog,
  message = ...
'SOLVE.M: Sorry!  Rank(CC) needs to be at least n! Cannot solve for PP. ';
  if DISPLAY_IMMEDIATELY, disp(message); end;
  warnings = [warnings;message];
else
  CC_plus = pinv(CC);
  CC_0 = (null(CC'))';
  Psi_mat   = [ zeros(l_equ-n_endog,m_states)
              FF - JJ*CC_plus*AA           ];
  Gamma_mat = [ CC_0 * AA
                JJ*CC_plus*BB - GG + KK*CC_plus*AA ];
  Theta_mat = [ CC_0 * BB
                KK*CC_plus*BB - HH                 ];
  Xi_mat    = [ Gamma_mat,     Theta_mat
                eye(m_states), zeros(m_states) ];
  Delta_mat = [ Psi_mat,       zeros(m_states)
                zeros(m_states), eye(m_states) ];
  [Xi_eigvec,Xi_eigval] = eig(Xi_mat,Delta_mat);
  if rank(Xi_eigvec)<m_states,
     message = ...
 'SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP.         ';
     if DISPLAY_IMMEDIATELY, disp(message); end;
     warnings = [warnings;message];
  else
    [Xi_sortabs,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
    Xi_sortvec = Xi_eigvec(1:2*m_states,Xi_sortindex);
    Xi_sortval = diag(Xi_eigval(Xi_sortindex,Xi_sortindex));
    Xi_select = 1 : m_states;
    if imag(Xi_sortval(m_states))~=0,
      if (abs( Xi_sortval(m_states) - conj(Xi_sortval(m_states+1)) )...
              < TOL),
      % NOTE: THIS LAST LINE MIGHT CREATE PROBLEMS, IF THIS EIGENVALUE 
      % OCCURS MORE THAN ONCE!!
      % IF YOU HAVE THAT PROBLEM, PLEASE TRY MANUAL ROOT SELECTION.  
        drop_index = 1;
        while (abs(imag(Xi_sortval(drop_index)))>TOL) && ...
                (drop_index < m_states),
          drop_index = drop_index + 1;
        end;
        if drop_index >= m_states,
          message = ...
['SOLVE.M: You are in trouble. You have complex eigenvalues, and I cannot'
'   find a real eigenvalue to drop to only have conjugate-complex pairs.'
'   Put differently: your PP matrix will contain complex numbers. Sorry!'];
          if DISPLAY_IMMEDIATELY, disp(message); end;
          warnings = [warnings;message];
          if m_states == 1,
            message = ...
['   TRY INCREASING THE DIMENSION OF YOUR STATE SPACE BY ONE!            '
'   WATCH SUNSPOTS!                                                     '];                     
            if DISPLAY_IMMEDIATELY, disp(message); end;
            warnings = [warnings;message];
          end;
        else
          message = ...
['SOLVE.M: I will drop the lowest real eigenvalue to get real PP.        '
'         I hope that is ok. You may have sunspots.                     ']; 
          if DISPLAY_IMMEDIATELY, disp(message); end;
          warnings = [warnings;message];
          Xi_select = [ 1: (drop_index-1), (drop_index+1):(m_states+1)];
        end;
      end;
    end;
  if MANUAL_ROOTS,
      message = ...
['SOLVE.M: You have chosen to select roots manually.  I am crossing my   '
'         fingers that you are doing it correctly.  In particular,      '
'         you should have defined Xi_manual.  Type help solve           '
'         and inspect SOLVE.M to get further information on how to do it'];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
      if exist('Xi_manual'),
         Xi_select = Xi_manual;
      else
         message = ...
['SOLVE.M: You have not defined Xi_manual.  Either define it or turn off '
'         the manual roots selection procedure with                     '
'         MANUAL_ROOTS = 0                                              '
'         Right now, I better let your calculations crash - sorry!      '
'         If you get results, they are based on previous calculations.  '];
         disp(message);
         warnings = [warnings;message];
      end;
  else
      if max(Xi_select) < 2*m_states,
        if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
          message = ...
['SOLVE.M: You may be in trouble. There are stable roots NOT used for PP.'
'         I have used the smallest roots: I hope that is ok.            '  
'         If not, try manually selecting your favourite roots.          '
'         For manual root selection, take a look at the file solve.m    '
'         Watch out for sunspot solutions.                              '];
          if DISPLAY_IMMEDIATELY, disp(message); end;
          warnings = [warnings;message];
        end; 
      end;
    if max(abs(Xi_sortval(Xi_select)))  > 1 + TOL,
      message = ...
['SOLVE.M: You may be in trouble.  There are unstable roots used for PP. '
'         Keep your fingers crossed or change your model.               '];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
    end;
    if abs( max(abs(Xi_sortval(Xi_select))) - 1  ) < TOL,
      message = ...
['SOLVE.M: Your matrix PP contains a unit root. You probably do not have '
'         a unique steady state, do you?  Should not be a problem, but  '
'         you do not have convergence back to steady state after a shock'
'         and you should better not trust long simulations.             '];
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
    end;

    Lambda_mat = diag(Xi_sortval(Xi_select));
    Omega_mat  = [Xi_sortvec((m_states+1):(2*m_states),Xi_select)];
    if rank(Omega_mat)<m_states,
      message = ...
'SOLVE.M: Sorry! Omega is not invertible. Cannot solve for PP.          ';
      if DISPLAY_IMMEDIATELY, disp(message); end;
      warnings = [warnings;message];
    else
      PP = Omega_mat*Lambda_mat/Omega_mat;
      PP_imag = imag(PP);
      PP = real(PP);
      if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
        message = ...
['SOLVE.M: PP is complex.  I proceed with the real part only.            '  
'         Hope that is ok, but you are probably really in trouble!!     '
'         You should better check everything carefully and be           '
'         distrustful of all results which follow now.                  '];
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
      end;
      RR = - CC_plus*(AA*PP+BB);
      if ny>0
          PM = (FF-JJ*CC\AA);
          if rank(PM)<Nx+ny
              Sylv = 0;
          end
      else
          if rank(FF) < nx
              Sylv = 0;
          end
      end
      if Sylv
          if ny>0
            Anew = PM \ (FF*PP+GG+JJ*RR-KK*CC\AA);
            Bnew = NN;
            Cnew = PM \ (JJ*CC\DD*N+KK*CC\DD-LL*NN-MM);
            QQ = sylvester(Anew,Bnew,Cnew);
            SS = -CC \(AA*QQ+DD);
          else
            Anew = FF \ (FF*PP+GG);
            Bnew = NN;
            Cnew = FF \ (-LL*NN-MM);
            QQ = sylvester(Anew,Bnew,Cnew);
            SS = [];         
          end
      else
          % EVERYTHING FROM HERE DOWN IS SOLVING FOR QQ & SS
          VV = [ kron(eye(k_exog),AA),   kron(eye(k_exog),CC)
                 kron(NN',FF)+kron(eye(k_exog),(FF*PP+JJ*RR+GG)), ...
                                         kron(NN',JJ)+kron(eye(k_exog),KK) ];
          if ( (rank(VV) < k_exog*(m_states+n_endog)) && ...
                                           (~IGNORE_VV_SING) ),
            message = ...
    ['SOLVE.M: Sorry! V is not invertible.  Cannot solve for QQ and SS. You  '
    '         can try setting IGNORE_VV_SING = 1 and wish for the best...   '];
            if DISPLAY_IMMEDIATELY, disp(message); end;
            warnings = [warnings;message];
          else
            if ( (rank(VV) < k_exog*(m_states+n_endog)) ),
              message = ...
    ['SOLVE.M: Warning! V is not invertible.  However, you have set          '
    '         IGNORE_VV_SING = 1, and thus, since you have told me to       '
    '         ignore this, I will proceed.  Keep your fingers crossed...    '];
              if DISPLAY_IMMEDIATELY, disp(message); end;
              warnings = [warnings;message];
            end;
            LLNN_plus_MM = LL*NN + MM;
            QQSS_vec = - VV \ [ DD(:)
                                LLNN_plus_MM(:) ];
            if max(abs(QQSS_vec)) == Inf,
               message = ...
    ['SOLVE.M: You probably are in trouble!  QQ or SS contain undefined      '
    '         entries! Most likely, the matrix VV is not invertible.        '];
               if DISPLAY_IMMEDIATELY, disp(message); end;
               warnings = [warnings;message];
            end;           
            QQ = reshape(QQSS_vec(1:m_states*k_exog),m_states,k_exog);
            SS = reshape(QQSS_vec((m_states*k_exog+1):...
                ((m_states+n_endog)*k_exog)),n_endog,k_exog);
            WW = [ eye(m_states)         , zeros(m_states,k_exog)
                   RR*pinv(PP)           , (SS-RR*pinv(PP)*QQ) 
                   zeros(k_exog,m_states), eye(k_exog)            ];
          end;       
      end;
    end;
  end;
end;
  
% find constant terms
%  construct aggregated matrices
if ny>0
    UU = -(FF*PP+GG+JJ*RR+FF-(JJ+KK)*(CC\AA)) \ ...
         (TT+(FF*QQ+JJ*SS+LL)*(NN*Z0-Z0)-(JJ+KK)*(CC\WWW));
    VVV = - CC \ (WWW+AA*UU);
else
    UU = -(FF*PP+FF+GG) \ (TT+(FF*QQ+LL)*(NN*Z0-Z0));
    VVV = [];
end

end;










































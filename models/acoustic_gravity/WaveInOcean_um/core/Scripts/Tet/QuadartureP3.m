%% Enriched P3
syms alpha;
 
%alpha = (3 - sqrt(3*(sqrt(sym(2))-1)))/6

%beta = (2+sqrt(sym(2)))/6; %Geevers

beta = 1/3 + 2*sym(sqrt(7))/21; %Joly - Cohen - Roberts - Tordjman

disp('value of beta')
eval(beta)

w_beta = 2 / (90* beta*(1-beta)^2);

R = beta*(1-beta)*(beta-alpha)*(beta-1+alpha) + (1-beta)*(1+beta)*((1-beta)/2-alpha)*((1-beta)/2-1+alpha)/2;
R = simplify(w_beta*expand(R));
L = simplify(expand((alpha - alpha^2)/3 + (alpha^2 - alpha -1)/6 + 1/5 - 1/15));

P = simplify(expand(L-R));

%Alpha solution to L-R = 0 
root_alpha = solve(P,alpha);

alpha_true = simplify(expand(min(root_alpha)));

disp('value of alpha')
eval(alpha_true)

%On calcul w_alpha avec la fonction x (1-x) =  (lambda_1)(1-lambda_1)
w_alpha = simplify(expand((1/6 - w_beta*(beta*(1-beta) + (1-beta)*(1+beta)/2))/(4*alpha_true*(1-alpha_true))));

w_s = simplify(expand(1 - 6*w_alpha-3*w_beta))/3;

disp('value of ws, walpha, wbeta')
eval(w_s)
eval(w_alpha)
eval(w_beta)


%  beta = sym(0.5853);
%  alpha_true = sym(0.2935);
%  w_s = sym(0.0148);
%  w_alpha = sym(0.0488);
%  w_beta = sym(0.2208);

 
disp('integral check')

%Verification sur (lambda_1)(1-lambda_1)lambda_2 (should be 1/20)
approximate_integral_1 = ( w_beta*(beta*(1-beta)*(1/2-beta/2) ...
                              + beta*(1/2-beta/2)*(1/2+beta/2) ...
                              + (1/2-beta/2)*(1/2-beta/2)*(1/2+beta/2)) ... 
                        + w_alpha*(alpha_true*(1-alpha_true)^2 ...
                                  + alpha_true^2*(1-alpha_true)));

 
%Verification sur lambda_1  (should be 1/3)
approximate_integral_2 = ( w_beta*(beta ...
                              + (1/2-beta/2) ...
                              + (1/2-beta/2)) ... 
                        + w_alpha*(2*alpha_true ...
                                  + 2*(1-alpha_true)) ...
                        + w_s);

%Verification sur lambda_1 lambda_2 (should be 1/12)
approximate_integral_3 = ( w_beta*(2*beta*(1/2-beta/2) ...
                              + (1/2-beta/2)^2) ... 
                        + w_alpha*(2*alpha_true*(1-alpha_true)));


%Verification sur (lambda_1)^2 (should be 1/6)
approximate_integral_4 = ( w_beta*(beta^2 ...
                              + (1/2-beta/2)^2 ...
                              + (1/2-beta/2)^2) ... 
                        + w_alpha*(2*alpha_true^2+2*(1-alpha_true)^2) ...
                        + w_s);

%Verification sur lambda_1^2 * lambda_2 * lambda_3  (should be 1/180)
approximate_integral_5=  w_beta*(beta^2*(1/2-beta/2)^2 ...
                              + 2*beta*(1/2-beta/2)^2*(1/2-beta/2));


%Verification sur (lambda_1)^3 (should be 1/10)
approximate_integral_6 = ( w_beta*(beta^3 ...
                              + (1/2-beta/2)^3 ...
                              + (1/2-beta/2)^3) ... 
                        + w_alpha*(2*alpha_true^3+2*(1-alpha_true)^3) ...
                        + w_s);

%Verification sur (lambda_1)^4 (should be 1/15)
approximate_integral_7 = ( w_beta*(beta^4 ...
                              + (1/2-beta/2)^4 ...
                              + (1/2-beta/2)^4) ... 
                        + w_alpha*(2*alpha_true^4+2*(1-alpha_true)^4) ...
                        + w_s);


%Verification sur lambda_1^3 * lambda_2 * lambda_3  (should be 2*3!/7!)
approximate_integral_8=  w_beta*(beta^3*(1/2-beta/2)*(1/2-beta/2)   ...
                              + 2*(1/2-beta/2)^3*beta*(1/2-beta/2));


%Verification sur (lambda_1)^5 (should be 2/42)
approximate_integral_10 = ( w_beta*(beta^5 ...
                              + (1/2-beta/2)^5 ...
                              + (1/2-beta/2)^5) ... 
                        + w_alpha*(2*alpha_true^5+2*(1-alpha_true)^5) ...
                        + w_s);


eval(approximate_integral_1 - 1/20)
eval(approximate_integral_2 - 1/3)
eval(approximate_integral_3 - 1/12)
eval(approximate_integral_4 - 1/6)
eval(approximate_integral_5 - 1/180)
eval(approximate_integral_6 - 1/10)
eval(approximate_integral_7 - 1/15)
eval(approximate_integral_8 - 2*factorial(3)/factorial(7))
eval(approximate_integral_10 - 1/21)

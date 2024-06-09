/*
 * Example 1 from Tuhin G M Al Mamun
 *
 * Tuhin G M Al Mamun
 * Ph.D. student at Hannam University
 * Specializing in macroeconomics with a focus on monetary and fiscal policy dynamics, as well as environmental economics.
 * Proficient in Python, Julia, MATLAB-Dynare, VSCode, EViews, and Stata.
 *
 * Research explores the interactions between monetary and fiscal policies and their impacts on economic stability and growth.
 * Passionate about understanding the relationship between economic policies and environmental sustainability.
 * Employing advanced statistical methods and computational models to provide valuable insights for policymakers.
 */

var Y, C, K, IK, IH, H, HE, LE, E, W, R, A, B,L;
varexo e, v;

parameters alpha, beta, deltak, deltah, deltahe, gamma, theta, rhoA, rhoB, eta, omega;

alpha = 0.35;
beta = 0.97;
deltak = 0.06;
deltah = 0.01;
deltahe = 0.01; // depreciation rate for educated human capital
gamma = 0.40;
theta = 0.80;
rhoA = 0.95;
rhoB = 0.95;
eta = 0.5;
omega = 0.5; // weight parameter for educated human capital

model;
C = (gamma/(1-gamma))*(1-L-E)*H*W;
1 = beta*((C/C(+1))*(R(+1)+(1-deltak)));
Y = A*(K(-1)^alpha)*((L*H + LE*HE)^(1-alpha)); // Updated production function
K = (Y-C)+(1-deltak)*K(-1);
IK = Y-C;
H = IH+(1-deltah)*H(-1);
IH = B*(E)^theta;
HE = IH+(1-deltahe)*HE(-1); // Equation for educated human capital
LE = omega * L; // Assuming a fraction of labor is educated
(1-gamma)/((1-L-E)*theta*B*(E)^(theta-1)) = beta*((gamma*W(+1)*L(+1))/(C(+1)*((1-gamma)*(1-deltah)/(1-L(+1)-E(+1)*theta*B*(E(+1))^(theta-1)))));
W = (1-alpha)*A*(K(-1)^alpha)*((L*H + LE*HE)^(-alpha)); // Updated wage equation
R = alpha*A*(K(-1)^(alpha-1))*((L*H + LE*HE)^(1-alpha)); // Updated rental rate of capital
log(A) = rhoA*log(A(-1)) + e;
log(B) = rhoB*log(B(-1)) + v;

end;

initval;
Y = 1;
C = 0.8;
L = 0.3;
K = 3.5;
IK = 0.2;
E = 0.15;
IK = 0.15^0.8;
H = IK/deltah;
HE = IK/deltahe; // Initialize educated human capital
LE = omega * L; // Initialize educated labor
W = (1-alpha)*Y/L;
R = alpha*Y/K;
A = 1;
B = 1;
e = 0;
v = 0;
end;

steady;

// Blanchard-Kahn Conditions
check;

shocks;
var e; stderr 0.009;
var v; stderr 0.009;
end;

//Observable Variables
varobs C ;

calib_smoother(diffuse_filter,datafile=vars, filtered_vars, filter_step_ahead = [3:4])Y, K, IK, IH, H, L, E, W, R, A, B,L;
// Estimation command
estimation(datafile=vars, filtered_vars, mh_replic= 2000, forecast=5, mh_conf_sig = 0.95, mh_nblocks = 2, mh_jscale = 3.5);

// Stochastic simulation
stoch_simul;
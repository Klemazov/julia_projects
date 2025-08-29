using DifferentialEquations
using .kinetic_base
using Plots

#example solving diff eqa

function F1d(u, p, t)
    Eₐ, A, β, = p;
    K(A, Eₐ,β*t)*F2(u)     
end

function F1B1(du, u, p, t)
    Eₐ, A, β, = p;
    du[1] = K(A, Eₐ, β*t)*F1(u[1]); 
    du[2] = K(A, Eₐ, β*t)*B1(du[1]);    
end

tspan = (0,10);
u0 = [1e-2, 1e-2];
p = [1.0, 1, 1];
prob = ODEProblem(F1B1, u0, tspan, p);
sol = solve(prob);
sol
plot(sol)
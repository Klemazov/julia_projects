using DifferentialEquations
using .kinetic_base


#example solving diff eqa

function F1d(u, p, t)
    Eₐ, A, β, = p;
    K(A, Eₐ,β*t)*F1(u)     
end

tspan = (0,1);
u0 = 1e-2;
p = [1.0, 100.0, 20.];
prob = ODEProblem(F1d, u0,tspan, p);
sol = solve(prob);

module kinetic_base
export R, K, F1

const R = 8.314;

# Constant of reaction
K(A::Number,Eₐ::Number,T::Number) = A*exp(-Eₐ/(R*T));



# for solving via diff equations
F1(α) = 1- α; #first order

F2(α) = (1- α)^2; #second order

Fn(α,n) = (1-α)^n; #n-th order


B1(α) = (1-α)*α; # Prout-Tompkins equation

Bna(α,n,m) = (1-α)^n*α^m; # the extended Prout-Tompkins equation

C1(α, Kcat) = F1(α)*(1+Kcat*α) #Kama-Sourour auto-catalysis reaction n = 1, m = 1

Cn(α,n, Kcat) = Fn(α,n)*(1+Kcat*α) #Kama-Sourour auto-catalysis reaction m = 1

Cnm(α,n,m, Kcat) = Fn(α,n)*(1+Kcat*α^m) #Kama-Sourour auto-catalysis reaction




end
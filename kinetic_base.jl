const R = 8.314;

# Constant of reaction
K(A,Eₐ,T) = A*exp(-Eₐ/(R*T));


# functions of kinetic models


# SIMPLE ORDER EQUATIONS

"first order equation F1(α) = 1 - α"
function  F1(α::T) where {T} 
@. (one(eltype(α)) - α);
end

"second order equation F2(α) = (1 - α)^2"
function  F2(α::T) where {T}  
@. (one(eltype(α)) - α)^2;
end

"n-th order equation Fn(α,n) = (1 - α)^n"
function  Fn(α::T, n) where {T} 
@. (one(eltype(α)) - α)^n;
end



# PROUT-TOMPKINS EQUATIONS

"Prout-Tompkins equation B1(α) = (1-α)⋅α"
function  B1(α::T) where {T} 
@. (one(eltype(α)) - α)*α;
end


"Extended Prout-Tompkins equation Bnm(α,n,m) = (1-α)ⁿ ⋅ αᵐ"
function  B1(α::T,n,m) where {T} 
@. (one(eltype(α)) - α)^n*α^m;
end


#AUTOCATALYS EQUATIONS

"Kama-Sourour auto-catalysis reaction n = 1, m = 1  C1(α) = F1(α)⋅(1+Kcat⋅α)"
function C1(α::T, Kcat) where T
    @. F1(α)*(one(eltype(α)+Kcat*α))
end

"Kama-Sourour auto-catalysis reaction m = 1 Cn(α) = Fn(α)⋅(1+Kcat⋅α)"
function Cn(α::T, n, Kcat) where T
    @. Fn(α,n)*(one(eltype(α)+Kcat*α))
end

"Kama-Sourour auto-catalysis reaction Cnm(α) = Fn(α)⋅(1+Kcat⋅α^m)"
function Cnm(α::T, n, m, Kcat) where T
    @. Fn(α,n)*(one(eltype(α)+Kcat*α^m))
end


#LINEAR REGRESSION

## some test data\
function cumtrapz(X::T, Y::T) where {T <: AbstractVector}
  # Check matching vector length
  @assert length(X) == length(Y)
  # Initialize Output
  out = similar(X)
  out[1] = 0
  # Iterate over arrays
  for i in 2:length(X)
    out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
  end
  # Return output
  out
end




using Plots;
using CSV, DataFrames, StringEncodings
df = CSV.read("test_data\\km5.txt", DataFrame; normalizenames=true, header = 27);
df = rename(df, [:temperature, :time, :mass]);

time = df.time
temperature = df.temperature
mass = df.mass

function fraction(x)
    return @. (abs(x-x[1]))/x[1]
end

alpha = fraction(mass);

dadt = abs.(diff(alpha)); #TODO rework abs mock

push!(dadt, dadt[end]);

falpha = B1(alpha);
# utils functions for linear regression for 1 step kinetic
# ln(dadt/f(a))  = lnA-Ea/RT -> Y = a+bX where Y = ln(dadt/f(a)); a = ln(A); b = Ea; X = -1/RT
# try to find a and b


function Y(dadt, falpha)
    return @. log(dadt/falpha)
end

function minusreverseRT(T)
    return @. (-one(eltype(T))/(R*T))
end


Yalpha = Y(dadt, falpha)
X = minusreverseRT(temperature);

## USING GENERAL LINEAR MODEL PACKAGE


using GLM
data = DataFrame(X = X, Y = Yalpha);
ols = lm(@formula(Y ~ X), data)

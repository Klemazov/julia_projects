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
using Plots;
using CSV, DataFrames, StringEncodings
df = CSV.read("test_data\\PLA5MIN.txt", DataFrame; normalizenames=true, header = 48, delim = "\t");
df = rename(df, [:time, :temperature, :heatflow]);

time = df.time[6400:11250]
temperature = df.temperature[6400:11250]
heatflow = df.heatflow[6400:11250] #kind of dadt




# utils functions for linear regression for 1 step kinetic
# ln(dadt/f(a))  = lnA-Ea/RT -> Y = a+bX where Y = ln(dadt/f(a)); a = ln(A); b = Ea; X = -1/RT
# try to find a and b


function Y(dadt, falpha)
    return @. log(dadt/falpha)
end

function minusreverseRT(T)
    return @. (-one(eltype(T))/(R*T))
end

## USING GENERAL LINEAR MODEL PACKAGE


using GLM


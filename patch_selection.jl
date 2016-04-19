
# patches are characterized the by following four parameters
type Patch
    β::Float64 # MC's probability of predation
    λ::Float64 # MC's probability of finding food
    α::Int64 # MC's cost to visit
    Y::Int64 # MC's value of food
end

# an Ecosystem is essentially an array of patches
typealias Ecosystem Array{Patch,1};

# the patches from page 54
myplace = Ecosystem([
    Patch(0., 0., 1, 0),
    Patch(0.004, 0.4, 1, 3), 
    Patch(0.02, 0.6, 1, 5)
    ]);

x_c = 3::Int64; # critical energy threshold
C = 10::Int64; # max capacity
T = 20::Int64; # time horizon

function xprime(x::Int64, α::Int64, y::Int64, x_critical::Int64, C::Int64)
    return clamp(x - α + y, x_critical, C)
end;

function xprimeprime(x::Int64, α::Int64, x_critical::Int64, C::Int64)
    return clamp(x - α, x_critical, C)
end;


function findOptimalSelections(
    aplace::Ecosystem,
    C::Int64,
    T::Int64,
    x_c::Int64)
    Fitness = zeros(Float64, C+1, T)
    Selection = zeros(Int64, C+1, T)
    for t = T:-1:1
        for x = 0:C
            if t == T
                Fitness[x+1, t] = x > x_c ? 1.0 : 0.0
            elseif x <= x_c 
                Fitness[x+1, t] = 0.0
            else 
                patchoptions = [(1-P.β)*
                    (P.λ*Fitness[xprime(x, P.α, P.Y, x_c, C)+1,t+1] +
                    (1.0-P.λ)*Fitness[xprimeprime(x, P.α, x_c, C)+1,t+1])
                    for P in aplace]
                (Fitness[x+1, t], Selection[x+1, t]) = findmax(patchoptions)
            end
        end
    end
    return Selection, Fitness
end;

(Selection, Fitness) = findOptimalSelections(myplace, C, T, x_c);

## compare to Mangel and Clark page 55
Fitness[5:11,19:-1:1]



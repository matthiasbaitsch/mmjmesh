module MMJBase

# Exports
## dimensions.jl
export D0, D1, D2, D3
## functions.jl
export atol
## macros.jl
export @abstractmethod
export @notimplemented
export @unreachable

# Parts
include("functions.jl")
include("macros.jl")

end
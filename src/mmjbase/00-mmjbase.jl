module MMJBase

# Exports
## functions.jl
export atol
## macros.jl
export @abstractmethod
export @notimplemented
export @unreachable

# Parts
include("functions.jl")
include("macros.jl")
include("seqintset.jl")

end
module MMJBase

# Exports
## seqintset.jl
export SeqIntSet
## functions.jl
export atol, tomatrix, FromType, ROWS, COLS
## macros.jl
export @abstractmethod
export @notimplemented
export @unreachable

# Parts
include("functions.jl")
include("macros.jl")
include("seqintset.jl")

end
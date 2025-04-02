module MMJBase


# Packages needed
import Symbolics
import SymbolicUtils


# Exports

## macros.jl
export NotImplementedError
export @abstractmethod
export @notimplemented
export @unreachable

## functions.jl
export atol, tomatrix, FromType, ROWS, COLS, pdim, gdim

## seqintset.jl
export SeqIntSet

## symbolics.jl
export rationalize!, integerize, integerize!


# Parts
include("macros.jl")
include("functions.jl")
include("seqintset.jl")
include("symbolics.jl")

end
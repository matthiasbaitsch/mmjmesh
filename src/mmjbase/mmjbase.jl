module MMJBase

export atol
export @abstractmethod
export @notimplemented
export @unreachable

atol(Float64) = 1e-10

include("macros.jl")

end
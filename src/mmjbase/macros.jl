# Macros for error reporting, borrowed from Gridap.jl

macro abstractmethod(message="This function is an abstract method")
  quote
    error($(esc(message)))
  end
end

macro notimplemented(message="This function is not yet implemented")
  quote
    error($(esc(message)))
  end
end

macro unreachable(message="This line of code cannot be reached")
  quote
    error($(esc(message)))
  end
end

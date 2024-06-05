# Macros for error reporting, borrowed from Gridap.jl

struct NotImplementedError <: Exception
  val
  msg::AbstractString
  NotImplementedError(@nospecialize(val)) = (@noinline; new(val, ""))
  NotImplementedError(@nospecialize(val), @nospecialize(msg)) = (@noinline; new(val, msg))
end


macro abstractmethod(message="This function is an abstract method")
  quote
    error($(esc(message)))
  end
end

macro notimplemented(message="This function is not yet implemented")
  quote
    throw(NotImplementedError($(esc(message))))
  end
end

macro unreachable(message="This line of code cannot be reached")
  quote
    error($(esc(message)))
  end
end

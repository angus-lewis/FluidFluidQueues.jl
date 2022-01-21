module MyModule

export f, g, h

f(x) = x^3
g(x) = f(f(x))
h(x) = 1/f(x)

end
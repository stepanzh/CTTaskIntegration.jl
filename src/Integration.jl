module Integration

export Gauss, Kronrod, Midpoint, NewtonCotesClosed, Trapezoid, Simpson
export integrate

abstract type AbstractMethod end

methods() = (
    Gauss,
    Kronrod,
    Midpoint,
    NewtonCotesClosed,
    Trapezoid,
    Simpson,
)

"""
    integrate(f, a, b; method)

Интегрирование `f` на отрезке [`a`, `b`] методом `method`.

Список доступных методов: см. `Integration.methods()`
"""
integrate(f, a, b; method)

"""
    integrate(f, a, b, nnodes; method)

Интегрирование `f` на отрезке [`a`, `b`] составным вариантом метода `method`.
При этом отрезок предварительно разбивается `nnodes` равноотстоящими узлами.

Также см. [`integrate(f, a, b; method)`](@ref)
"""
integrate(f, a, b, nnodes; method)

# Диспетчеризация интерфейса
integrate(args...; method) = __integrate_impl(method, args...)

# Общий метод для составной формулы
function __integrate_impl(method::AbstractMethod, f, a, b, nnodes)
    x = range(a, b; length=nnodes-1)
    int = 0.0
    @views for (x₁, x₂) in zip(x[1:end-1], x[2:end])
        int += __integrate_impl(method, f, x₁, x₂)
    end
    return int
end


"Формула средних прямоугольников"
struct Midpoint <: AbstractMethod end

__integrate_impl(method::Midpoint, f, a, b) = (b-a) * f((b+a)/2)

# составной метод для формулы прямоугольников
function __integrate_impl(method::Midpoint, f, a, b, nnodes)
    h = (b - a) / (nnodes - 1)
    x = range(a + h/2, b; step=h)
    int = h * sum(f, x)
    return int
end

"Формула трапеций"
struct Trapezoid <: AbstractMethod end
__integrate_impl(method::Trapezoid, f, a, b) = error("Не имплементирован")
__integrate_impl(method::Trapezoid, f, a, b, nnodes) = error("Не имплементирован")

"Формула Симпсона"
struct Simpson <: AbstractMethod end
__integrate_impl(method::Simpson, f, a, b) = error("Не имплементирован")
__integrate_impl(method::Simpson, f, a, b, nnodes) = error("Не имплементирован")


# Формулу Ньютона-Котса реализуйте в общем виде
# При этом используйте Лежандровскую матрицу системы на веса
# и равноотстоящие узлы

"Общая формула Ньютона-Котса закрытого типа с `nnodes` узлами"
struct NewtonCotesClosed <: AbstractMethod end
__integrate_impl(method::NewtonCotesClosed, f, a, b, nnodes) = error("Не имплементирован")


# В качестве формулы Гаусса возьмите G_7
# https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula#Example
# Не забудьте пересчитать узлы и веса на отрезок [a, b]

struct Gauss <: AbstractMethod end
"Интгерирование методом Гаусса по 7 точкам."
__integrate_impl(method::Gauss, f, a, b) = error("Не имплементирован")


# В качестве формулы Кронрода возьмите K_15
# https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula#Example
# Не забудьте пересчитать узлы и веса на отрезок [a, b]

struct Kronrod <: AbstractMethod end
__integrate_impl(f, a, b) = error("Не имплементирован")

end # module

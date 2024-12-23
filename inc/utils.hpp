#pragma once

template <typename T>
T abs_val(T const a)
{
    return a > 0 ? a : -a;
}

template <typename T>
T lerp(T const a, T const b, double t)
{
    return a * (1 - t) + b * t;
};
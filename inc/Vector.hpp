#pragma once

#include <array>
#include <vector>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include "utils.hpp"
#include <cmath>

template <typename T, size_t Size>
class Vector
{
private:
    std::array<T, Size> values;

    using own_type = Vector<T, Size>;

public:
    Vector()
    {
        for (size_t i = 0; i < Size; i++)
        {
            values[i] = T(0);
        }
    };

    Vector (T value)
    {
        for (size_t i = 0; i < Size; i++)
        {
            values[i] = value;
        }
    };

    Vector(std::initializer_list<T> init)
    {
        std::copy(init.begin(), init.end(), values.begin());
    }
    Vector(const Vector &) = default;

    Vector operator=(const Vector &other)
    {
        if (this != &other)
        {
            for (size_t i = 0; i < Size; i++)
            {
                values[i] = other[i];
            }
        }
        return *this;
    };

    std::array<T, Size> &get_values()
    {
        return values;
    };

    size_t size() const
    {
        return Size;
    };

    T &operator[](size_t index)
    {
        return values[index];
    };

    const T &operator[](size_t index) const
    {
        return values[index];
    };

    own_type operator+(const own_type &other) const
    {
        Vector<T, Size> result;
        for (size_t i = 0; i < Size; i++)
        {
            result[i] = values[i] + other[i];
        }
        return result;
    };

    own_type operator-(const own_type &other) const
    {
        Vector<T, Size> result;
        for (size_t i = 0; i < Size; i++)
        {
            result[i] = values[i] - other[i];
        }
        return result;
    };

    own_type operator*(const T scalar) const
    {
        Vector<T, Size> result;
        for (size_t i = 0; i < Size; i++)
        {
            result[i] = values[i] * scalar;
        }
        return result;
    };

    bool operator==(const own_type &other) const
    {
        for (size_t i = 0; i < Size; i++)
        {
            if (values[i] != other[i])
            {
                return false;
            }
        }
        return true;
    };

    T dot(own_type const &other) const
    {
        T result = T();
        for (size_t i = 0; i < Size; i++)
        {
            result += (*this)[i] * other[i];
        }
        return result;
    };

    double norm_1() const
    {
        T result = T();
        for (T value : values)
        {
            result += std::max(value, -value);
        }
        return result;
    };

    double norm() const
    {
        T result = T();
        result = this->dot(*this);
        return std::pow(result, 0.5);
    };

    double norm_inf() const
    {
        T result = T();
        for (T value : values)
            result = std::max(result, std::max(value, -value));
        return result;
    }
};

template <typename T, size_t Size>
Vector<T, Size> linear_combination(const std::vector<Vector<T, Size>> &vecs, const std::vector<T> &coefficients)
{
    if (vecs.size() != coefficients.size())
    {
        throw std::invalid_argument("Vector and coefficients have to be of same size.");
    }
    Vector<T, Size> result;
    for (size_t i = 0; i < vecs.size(); i++)
    {
        result = result + vecs[i] * coefficients[i];
    }
    return result;
};

template <typename T, size_t Size>
Vector<T, Size> lerp(Vector<T, Size> const &u, Vector<T, Size> const &v, double t)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Vectors have to be of same size.");
    }
    Vector<T, Size> result;
    result = (u * (1 - t)) + (v * t);
    return result;
};

template <typename T, size_t Size>
T angle_cos(Vector<T, Size> &u, Vector<T, Size> &v)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Vectors have to be of same size.");
    }
    double mag_u = u.norm();
    double mag_v = v.norm();
    if (mag_u == 0 || mag_v == 0)
    {
        throw std::runtime_error("Error: Attempting to devide by 0.");
    }
    Vector<T, Size> normalised_u = u * static_cast<T>(1 / mag_u);
    Vector<T, Size> normalised_v = v * static_cast<T>(1 / mag_v);

    return normalised_u.dot(normalised_v);
};

template <typename T, size_t Size>
Vector<T, Size> cross_product(Vector<T, Size> &u, Vector<T, Size> &v)
{
    if (u.size() != 3 || v.size() != 3)
    {
        std::invalid_argument("Both vectors need to have a length of 3");
    }
    Vector<T, Size> result;
    result[0] = u[1] * v[2] - u[2] * v[1];
    result[1] = u[2] * v[0] - u[0] * v[2];
    result[2] = u[0] * v[1] - u[1] * v[0];
    return result;
};

template <typename T, size_t Size>
std::ostream &operator<<(std::ostream &os, Vector<T, Size> const &vec)
{
    os << "[";
    for (size_t i = 0; i < vec.size(); i++)
    {
        os << vec[i];
        if (i != vec.size() - 1)
        {
            os << ", ";
        }
    }
    os << "]" << std::endl;
    return os;
};
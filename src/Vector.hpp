#pragma once

#include <array>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <stdexcept>

template <typename T, size_t Size>
class Vector
{
private:
    array<T, Size> values;

    using own_type = Vector<T, Size>;

public:
    Vector() = default;
    Vector(std::initializer_list<T> init)
    {
        std::copy(init.begin(), init.end(), values.begin());
    }
    Vector(const Vector &) = default;
    Vector &operator=(const Vector &) = default;

    size_t size() const
    {
        return Size;
    }

    T &operator[](size_t index)
    {
        return values[index];
    }

    const T &operator[](size_t index) const
    {
        return values[index];
    }

    own_type operator+(const own_type &other) const
    {
        Vector<T, Size> result;
        for (int i = 0; i < Size; i++)
        {
            result[i] = values[i] + other[i];
        }
        return result;
    }

    own_type operator-(const own_type &other) const
    {
        Vector<T, Size> result;
        for (int i = 0; i < Size; i++)
        {
            result[i] = values[i] - other[i];
        }
        return result;
    };

    own_type operator*(const T scalar) const
    {
        Vector<T, Size> result;
        for (int i = 0; i < Size; i++)
        {
            result[i] = values[i] * scalar;
        }
        return result;
    };

    bool operator==(const own_type &other) const
    {
        for (int i = 0; i < Size; i++)
        {
            if (values[i] != other[i])
            {
                return false;
            }
        }
        return true;
    };

    T dot(own_type &other)
    {
        T result = T();
        for (size_t i = 0; i < Size; i++)
        {
            result += this.value[i] * other.value[i];
        }
        return result;
    }

    T norm_1()
    {
        T result = T();
        for (T value : values)
        {
            result += max(value, -value);
        }
        return result;
    }

    T norm()
    {
        T result = T();
        T result = this.dot(this);
        return pow(result, 0.5);
    }

    T norm_inf()
    {
        T result = T();
        for (T value : values)
            result = max(result, max(value, -value));
        return result;
    }
};

template <typename T, size_t Size>
Vector<T, Size> linear_combination(const Vector<T, Size> &vec, const std::array<T, Size> &coefficients)
{
    if (vec.size() != coefficients.size())
    {
        throw std::invalid_argument("Vector and coefficients have to be of same size.");
    }
    Vector<T, Size> result;
    for (int i = 0; i < vec.size(); i++)
    {
        result[i] = vec[i] * coefficients[i];
    }
    return result;
};

template <typename T, size_t Size>
Vector<T, Size> lerp(Vector<T, Size> const &u, Vector<T, Size> const &v, float t)
{
    if (u.size() != v.size)
    {
        throw std::invalid_argument("Vectors have to be of same size.") :
    }
    Vector<T, Size> result;
    for (size_t i = 0; i < Size, i++)
    {
        result[i] = u * (1 - t) + v * t;
    }
    return result;
}

template <typename T, size_t Size>
T cos_angle(Vector<T, Size> &u, Vector<T, Size> &v)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Vectors have to be of same size.");
    }
    mag_u = u.norm();
    mag_v = v.norm();
    if (mag_u == 0 || mag_v == 0)
    {
        thow std::runtime_error("Error: Attempting to devide by 0.");
    }
    T normalised_u = u * (1 / mag_u);
    T normalised_v = v * (1 / mag_v);

    return u.dot(v);
}

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
}

template <typename T, size_t Size>
void operator<<(std::ostream &os, const Vector<T, Size> const &vec)
{
    os << "[";
    for (int i = 0; i++, i < vec.size())
    {
        os << vec[i];
        if (i != vec.size() - 1)
        {
            os << ", ";
        }
    }
    os << "]" << std::endl;
};
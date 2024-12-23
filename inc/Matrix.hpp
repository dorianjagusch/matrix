#pragma once

#include <array>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include "Vector.hpp"

#define EPSILON 1e-6

template <typename T, size_t nRow, size_t nCol>
class Matrix
{

private:
    std::array<Vector<T, nRow>, nCol> values;
    size_t rows = nRow;
    size_t cols = nCol;
    bool isDeterminantCached = false;
    T determinant_cache = T();
    using own_type = Matrix<T, nRow, nCol>;

    void swap_rows(size_t row1, size_t row2)
    {
        if (row1 >= nRow || row2 >= nRow)
        {
            throw std::invalid_argument("Row index out of bounds.");
        }
        for (size_t i = 0; i < nCol; i++)
        {
            std::swap((*this)[i][row1], (*this)[i][row2]);
        }
    }

public:
    Matrix()
    {
        for (auto &col : values)
        {
            col = Vector<T, nRow>();
        }
    }

    Matrix(T value)
    {
        for (auto &col : values)
        {
            col = Vector<T, nRow>(value);
        }
    }

    Matrix(const std::array<std::array<T, nRow>, nCol> &init_values)
    {
        for (size_t i = 0; i < nCol; ++i)
        {
            for (size_t j = 0; j < nRow; ++j)
            {
                values[i][j] = init_values[i][j];
            }
        }
    }

    Matrix(own_type const &other)
    {
        for (size_t i = 0; i < nCol; ++i)
        {
            values[i] = other[i];
        }
    }

    own_type &operator=(own_type const &other)
    {
        if (this != &other)
        {
            for (size_t i = 0; i < nCol; ++i)
            {
                values[i] = other[i];
            }
        }
        return *this;
    }

    Matrix(std::initializer_list<std::initializer_list<T>> const init_list)
    {
        if (init_list.size() != nCol)
        {
            throw std::invalid_argument("Initializer list size does not match matrix column size");
        }

        size_t i = 0;
        for (const auto &col : init_list)
        {
            if (col.size() != nRow)
            {
                throw std::invalid_argument("Initializer list size does not match matrix row size");
            }

            std::copy(col.begin(), col.end(), values[i].get_values().begin());
            ++i;
        }
    }

    std::pair<size_t, size_t> size() const
    {
        return std::make_pair(rows, cols);
    }

    bool is_square() const
    {
        return rows == cols;
    }

    own_type operator+(own_type &other) const
    {
        own_type result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i] + other[i];
        }
        return result;
    }

    own_type operator-(own_type &other) const
    {
        own_type result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i] - other[i];
        }
        return result;
    };

    own_type operator*(T scalar) const
    {
        own_type result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i] * scalar;
        }
        return result;
    };

    Vector<T, nCol> operator*(Vector<T, nRow> &vec) const
    {
        Vector<T, nCol> result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i].dot(vec);
        }
        return result;
    };

    bool operator==(own_type &other) const
    {
        for (size_t i = 0; i < nCol; i++)
        {
            for (size_t j = 0; j < nRow; j++)
            {
                if (abs_val(values[i][j] - other[i][j]) > static_cast<T>(EPSILON))
                    return false;
            }
        }
        return true;
    };

    Vector<T, nRow> &operator[](size_t index)
    {
        return values[index];
    }

    const Vector<T, nRow> &operator[](size_t index) const
    {
        return values[index];
    }

    template <size_t nOtherCol>
    Matrix<T, nRow, nOtherCol> operator*(const Matrix<T, nCol, nOtherCol> &other) const
    {
        Matrix<T, nRow, nOtherCol> result;
        Matrix<T, nOtherCol, nCol> other_transposed = other.transpose();
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nOtherCol; j++)
            {
                result[j][i] = values[j].dot(other_transposed[i]);
            }
        }
        return result;
    }

    own_type identity_matrix() const
    {
        own_type result;
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nCol; j++)
            {
                result[i][j] = i == j ? 1 : 0;
            }
        }
        return result;
    };

    T trace() const
    {
        if (!is_square())
        {
            throw std::invalid_argument("Matrix must be square to calculate trace.");
        }
        T result = T();
        for (size_t i = 0; i < nRow; i++)
        {
            result += values[i][i];
        }
        return result;
    }

    Matrix<T, nCol, nRow> transpose() const
    {
        Matrix<T, nCol, nRow> result;
        for (size_t i = 0; i < nCol; i++)
        {
            for (size_t j = 0; j < nRow; j++)
            {
                result[j][i] = values[i][j];
            }
        }
        return result;
    };

    T determinant()
    {
        if (!is_square())
        {
            throw std::invalid_argument("Matrix must be square to calculate determinant.");
        }
        if (nRow == 1)
        {
            return values[0][0];
        }
        if (isDeterminantCached)
        {
            return determinant_cache;
        }
        return row_echelon().second;
    }

    std::pair<Matrix<T, nRow, nCol>, T> row_echelon()
    {
        Matrix<T, nRow, nCol> result = *this;
        size_t currentRow = 0, currentCol = 0;
        T det = T(1);


        while (currentRow < nRow && currentCol < nCol)
        {
            size_t i_max = currentRow;
            for (size_t i = currentRow + 1; i < nRow; i++)
            {
                if (abs_val(result[currentCol][i]) > abs_val(result[currentCol][i_max]))
                {
                    i_max = i;
                }
            }

            if (abs_val(result[currentCol][i_max]) < EPSILON)
            {
                currentCol++;
                det = 0;
                continue;
            }

            if (i_max != currentRow)
            {
                result.swap_rows(currentRow, i_max);
                det *= -1;
            }

            T pivot = result[currentCol][currentRow];
            det *= pivot;
            for (size_t j = 0; j < nCol; j++)
            {
                result[j][currentRow] /= pivot;
            }

            for (size_t i = 0; i < nRow; i++)
            {
                if (i == currentRow)
                {
                    continue;
                }
                T factor = result[currentCol][i];
                for (size_t j = 0; j < nCol; j++)
                {
                    result[j][i] -= factor * result[j][currentRow];
                }
            }

            currentRow++;
            currentCol++;
        }

        return {result, det};
    }

    Matrix<T, nRow, nCol> calculate_inverse()
    {
        Matrix<T, nRow, nCol> result = *this;
        Matrix<T, nRow, nCol> inverse = identity_matrix();
        size_t currentRow = 0, currentCol = 0;

        while (currentRow < nRow && currentCol < nCol)
        {
            size_t i_max = currentRow;
            for (size_t i = currentRow + 1; i < nRow; i++)
            {
                if (abs_val(result[currentCol][i]) > abs_val(result[currentCol][i_max]))
                {
                    i_max = i;
                }
            }

            if (abs_val(result[currentCol][i_max]) < EPSILON)
            {
                currentCol++;
                continue;
            }

            if (i_max != currentRow)
            {
                result.swap_rows(currentRow, i_max);
                inverse.swap_rows(currentRow, i_max);
            }

            T pivot = result[currentCol][currentRow];
            for (size_t j = 0; j < nCol; j++)
            {
                result[j][currentRow] /= pivot;
                inverse[j][currentRow] /= pivot;
            }

            for (size_t i = 0; i < nRow; i++)
            {
                if (i == currentRow)
                {
                    continue;
                }
                T factor = result[currentCol][i];
                for (size_t j = 0; j < nCol; j++)
                {
                    result[j][i] -= factor * result[j][currentRow];
                    inverse[j][i] -= factor * inverse[j][currentRow];
                }
            }

            currentRow++;
            currentCol++;
        }

        return inverse;
    }

    own_type inverse()
    {
        if (!is_square())
        {
            throw std::invalid_argument("Matrix must be square to calculate inverse.");
        }
        T det = determinant();
        if (det == 0)
        {
            throw std::runtime_error("Matrix is singular.");
        }

        own_type rref = row_echelon().first;
        if (rref != identity_matrix())
        {
            throw std::runtime_error("Matrix is not invertible.");
        }
        return calculate_inverse();
    }

    size_t rank()
    {
        Matrix<T, nRow, nCol> rref = row_echelon().first;
        size_t nonZeroRows = 0;
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nCol; j++)
            {
                if (rref[j][i] != 0)
                {
                    nonZeroRows++;
                    break;
                }
            }
        }
        return nonZeroRows;
    }
};

template <typename T, size_t nRow, size_t nCol>
std::ostream &operator<<(std::ostream &os, const Matrix<T, nRow, nCol> &mat)
{
    os << "\n";
    for (size_t col = 0; col < nCol; col++)
    {
        os << mat[col];
    }
    return os;
}

template <typename T, size_t nRow, size_t nCol>
Matrix<T, nRow, nCol> lerp(Matrix<T, nRow, nCol> const &m1, Matrix<T, nRow, nCol> const &m2, double t)
{
    Matrix<T, nRow, nCol> result;
    for (size_t i = 0; i < nCol; i++)
    {
        result[i] = lerp(m1[i], m2[i], t);
    }
    return result;
};

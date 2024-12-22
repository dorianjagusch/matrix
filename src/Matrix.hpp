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
    std::array<Vector<T, nCol>, nRow> values;
    size_t rows = nRow;
    size_t cols = nCol;
    isDeterminantCached = false;
    determinant_cache = T();
    using own_type = Matrix<T, nRow, nCol>;

public:
    Matrix() = default;
    Matrix(const std::array<std::array<T, nCol>, nRow> &init_values)
    {
        for (size_t i = 0; i < nRow; ++i)
        {
            for (size_t j = 0; j < nCol; ++j)
            {
                values[i][j] = init_values[i][j];
            }
        }
    }
    Matrix(const own_type &other)
    {
        for (size_t i = 0; i < nRow; ++i)
        {
            values[i] = other[i];
        }
    }

    own_type &operator=(const own_type &other)
    {
        if (this == &other)
        {
            return *this;
        }
        for (size_t i = 0; i < nRow; ++i)
        {
            values[i] = other[i];
        }
    }

    Matrix(std::initializer_list<std::initializer_list<T>> init_list)
    {
        if (init_list.size() != nRow)
        {
            throw std::invalid_argument("Initializer list size does not match matrix row size");
        }

        size_t i = 0;
        for (const auto &row : init_list)
        {
            if (row.size() != nCol)
            {
                throw std::invalid_argument("Initializer list size does not match matrix column size");
            }

            std::copy(row.begin(), row.end(), values[i].begin());
            ++i;
        }
    }

    std::pair<T, T> size() const
    {
        return std::make_pair(rows, cols);
    }

    bool is_square() const
    {
        return rows == cols;
    }

    own_type operator+(own_type &other)
    {
        Matrix<T, nCol, nRow> result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i] + other[i];
        }
        return result;
    }

    own_type operator-(own_type &other)
    {
        Matrix<T, nCol, nRow> result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i] - other[i];
        }
        return result;
    };

    own_type operator*(T scalar)
    {
        Matrix<T, nCol, nRow> result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i] * scalar;
        }
        return result;
    };

    Vector<T, nCol> operator*(Vector<T, nRow> &vec)
    {
        Vector<T, nCol> result;
        for (size_t i = 0; i < nCol; i++)
        {
            result[i] = values[i].dot(vec);
        }
        return result;
    };

    bool operator==(Matrix<T, nRow, nCol> &other)
    {
        for (size_t i = 0; i < nCol; i++)
        {
            for (size_t j = 0; j < nRow; j++)
            {
                if (abs(values[i][j] - other[i][j]) > EPSILON;)
                    return false;
            }
        }
        return true;
    };

    array<T, nRow> &operator[](size_t index)
    {
        return values[index];
    }

    template <size_t nOtherCol>
    Matrix<T, nRow, nOtherCol> operator*(Matrix<T, nCol, nOtherCol> &other)
    {
        if (cols != other.size().first)
        {
            throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
        }
        Matrix<T, nRow, nOtherCol> result;
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nOtherCol; j++)
            {
                result[i][j] = T();
                for (size_t k = 0; k < nCol; k++)
                {
                    result[i][j] += values[i][k] * other[k][j];
                }
            }
        }
        return result;
    };

    T trace()
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

    own_type transpose()
    {
        Matrix<T, nCol, nRow> result;
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nCol; j++)
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
        if (isDeterminantCached)
        {
            return determinant_cache;
        }
        if (nRow > 4)
        {
            isDeterminantCached = true;
            return row_echelon().second;
        }
        if (nRow == 1)
        {
            return values[0][0];
        }
        T result = T();
        for (size_t i = 0; i < nRow; i++)
        {
            /*
            For small matrices, it doesn't hurt too much to calculate the cofactors, but
            for larger matrices, it's better to use its row echelon form.
            It's technically a little more effecient in algorithmic complexity than RREF for 2 and 3.
            (O(n*n!) vs (On^3))
            */
            result += values[0][i] * cofactor(0, i);
        }
        return result;
    }

    T cofactor(size_t row, size_t col)
    {
        Matrix<T, nRow - 1, nCol - 1> result;
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nCol; j++)
            {
                if (i != row && j != col)
                {
                    result[i < row ? i : i - 1][j < col ? j : j - 1] = values[i][j];
                }
            }
        }
        return result.determinant() * ((row + col) % 2 == 0 ? 1 : -1);
    }


// TODO: Correct row echelon. My matrix is in column major order, so I need to potentially transpose it first and at the end (do research)
    std::pair<own_type, T> row_echelon()
    {
        own_type result = *this;
        size_t currentRow, currentCol = 0;
        T det = T(1);
        while (currentRow < nRow && currentCol < nCol)
        {
            size_t i_max = currentRow;
            for (size_t i = currentRow + 1; i < nRow; i++)
            {
                if (abs(values[i][currentCol]) > abs(values[i_max][currentCol])) // put my  own abs
                {
                    i_max = i;
                }
            }
            if (values[i_max][currentCol] == 0)
            {
                currentCol++;
            }
            else
            {
                std::swap(values[currentRow], values[i_max]);
                det *= -1;
                for (size_t i = currentRow + 1; i < nRow; i++)
                {
                    if (values[currentRow][currentCol] == 0)
                    {
                        if (!is_square())
                        {
                            std::runtime_error("Matrix is (nearly) singular.");
                        }
                        det = 0;
                        isDeterminantCached = true;
                    }
                    T f = values[i][currentCol] / values[currentRow][currentCol];
                    det *= values[currentRow][currentCol];
                    values[i][currentCol] = 0;
                    for (size_t j = currentCol + 1; j < nCol; j++)
                    {
                        values[i][j] = values[i][j] - values[currentRow][j] * f;
                    }
                }
                currentRow++;
                currentCol++;
            }
        }
        determant_cache = det;
        if (nRow == nCol)
        {
            isDeterminantCached = true;
        }
        return {result, det};
    };

    own_type calculate_inverse()
    {
        own_type result = *this;
        own_type inverse = identity_matrix(T, nRow, nCol);
        size_t currentRow, currentCol = 0;
        T det = T(1);
        while (currentRow < nRow && currentCol < nCol)
        {
            size_t i_max = currentRow;
            for (size_t i = currentRow + 1; i < nRow; i++)
            {
                if (abs(result[i][currentCol]) > abs(result[i_max][currentCol])) // put my  own abs
                {
                    i_max = i;
                }
            }
            if (result[i_max][currentCol] == 0)
            {
                currentCol++;
            }
            else
            {
                std::swap(result[currentRow], result[i_max]);
                std::swap(inverse[currentRow], inverse[i_max]);
                for (size_t i = currentRow + 1; i < nRow; i++)
                {
                    T f = result[i][currentCol] / result[currentRow][currentCol];
                    result[i][currentCol] = 0;
                    for (size_t j = currentCol + 1; j < nCol; j++)
                    {
                        result[i][j] = result[i][j] - result[currentRow][j] * f;
                        inverse[i][j] = inverse[i][j] - inverse[currentRow][j] * f;
                    }
                }
                currentRow++;
                currentCol++;
            }
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
        dims = size();
        if (rref != identity_matrix(T, size().first, size().second))
        {
            throw std::runtime_error("Matrix is not invertible.");
        }
        return calculate_inverse();
    }

    size_t rank(){
        size_t nonZeroRows = 0;
        for (size_t i = 0; i < nRow; i++)
        {
            for (size_t j = 0; j < nCol; j++)
            {
                if (values[i][j] != 0)
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
std::ostream &operator<<(std::ostream &os, Matrix<T, nRow, nCol> const &mat)
{
    os << "[ ";
    for (size_t i = 0; i < nRow; ++i)
    {
        os << mat[i];
        if (i != nRow - 1)
        {
            os << ", ";
        }
    }
    os << " ]" << std::endl;
    return os;
};

template <typename T, size_t nRow, size_t nCol>
Matrix<T, nRow, nCol> lerp(Matrix<T, nRow, nCol> &m1, Matrix<T, nRow, nCol> &m2, double t)
{
    Matrix<T, nRow, nCol> result;
    for (int i = 0; i < nRow; i++)
    {
        result[i] = values[i].lerp(other[i], t);
    }
    return result;
};

template <typename T, size_t nRow, size_t nCol>
Matrix<T, nRow, nCol> identity_matrix()
{
    Matrix<T, nRow, nCol> result;
    for (size_t i = 0; i < nRow; i++)
    {
        for (size_t j = 0; j < nCol; j++)
        {
            result[i][j] = i == j ? 1 : 0;
        }
    }
    return result;
};
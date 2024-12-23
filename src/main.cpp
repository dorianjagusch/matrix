#include <iostream>
#include "../inc/Vector.hpp"
#include "../inc/Matrix.hpp"
#include "../inc/utils.hpp"

int main(void)
{
    {
        Vector<double, 2> u = {2., 3.};
        Vector<double, 2> v = {5., 7.};
        std::cout << u + v << std::endl;
        // [7.0]
        // [10.0]
        std::cout << u - v << std::endl;
        // [-3.0]
        // [-4.0]
        std::cout << u * 2 << std::endl;
        // [4.0]
        // [6.0]
    }
    {
        Matrix<double, 2, 2> u = {{1., 2.}, {3., 4.}};
        Matrix<double, 2, 2> v = {{7., 4.}, {-2., 2.}};

        std::cout << u + v << std::endl;
        // [8.0, 6.0]
        // [1.0, 6.0]
        std::cout << u - v << std::endl;
        // [-6.0, -2.0]
        // [5.0, 2.0]
        std::cout << u * 2 << std::endl;
        // [2.0, 4.0]
        // [6.0, 8.0]
    }
    {
        Vector<float, 3> e1({1., 0., 0.});
        Vector<float, 3> e2({0., 1., 0.});
        Vector<float, 3> e3({0., 0., 1.});
        Vector<float, 3> v1({1., 2., 3.});
        Vector<float, 3> v2({0., 10., -100.});
        std::cout << linear_combination(std::vector<Vector<float, 3>>{e1, e2, e3}, std::vector<float>{10., -2., 0.5});
        // [10.]
        // [-2.]
        // [0.5]
        std::cout << linear_combination(std::vector<Vector<float, 3>>{v1, v2}, std::vector<float>{10., -2.});
        // [10.]
        // [0.]
        // [230.]
    }
    {
        std::cout << lerp(0., 1., 0.) << std::endl;
        // 0.0
        std::cout << lerp(0., 1., 1.) << std::endl;
        // 1.0
        std::cout << lerp(0., 1., 0.5) << std::endl;
        // 0.5
        std::cout << lerp(21., 42., 0.3) << std::endl;
        // 27.3
        std::cout << lerp(Vector<double, 2>({2., 1.}), Vector<double, 2>({4., 2.}), 0.3) << std::endl;
        // [2.6]
        // [1.3]
        std::cout << lerp(Matrix<double, 2, 2>({{2., 1.},
                                                {3., 4.}}),
                          Matrix<double, 2, 2>({{20., 10.},
                                                {30., 40.}}),
                          0.5)
                  << std::endl;
        // [[11., 5.5]
        // [16.5, 22.]]
    }
    {

        {
            Vector<float, 2> u = Vector<float, 2>({0., 0.});
            Vector<float, 2> v = Vector<float, 2>({1., 1.});
            std::cout << u.dot(v) << std::endl;
            // 0.0
        }
        {
            Vector<float, 2> u = Vector<float, 2>({1., 1.});
            Vector<float, 2> v = Vector<float, 2>({1., 1.});
            std::cout << u.dot(v) << std::endl;
            // 2.0
        }
        {
            Vector<float, 2> u = Vector<float, 2>({-1., 6.});
            ;
            Vector<float, 2> v = Vector<float, 2>({3., 2.});
            std::cout << u.dot(v) << std::endl;
            // 9.0
        }
    }
    {
        Vector<float, 3> u({0., 0., 0.});
        std::cout << u.norm_1() << "\t" << u.norm() << "\t" << u.norm_inf() << std::endl;
        // 0.0, 0.0, 0.0
        Vector<float, 3> v({1., 2., 3.});
        std::cout << v.norm_1() << "\t" << v.norm() << "\t" << v.norm_inf() << std::endl;
        // 6.0, 3.74165738, 3.0
        Vector<float, 2> w({-1., -2.});
        std::cout << w.norm_1() << "\t" << w.norm() << "\t" << w.norm_inf() << std::endl;
        // 3.0, 2.236067977, 2.0}
    }

    {
        {
            Vector<float, 2> u({1., 0.});
            Vector<float, 2> v({1., 0.});
            std::cout << angle_cos(u, v) << std::endl;
            // 1.0
        }
        {
            Vector<float, 2> u({1., 0.});
            Vector<float, 2> v({0., 1.});
            std::cout << angle_cos(u, v) << std::endl;
            // 0.0
        }
        {
            Vector<float, 2> u({-1., 1.});
            Vector<float, 2> v({1., -1.});
            std::cout << angle_cos(u, v) << std::endl;
            // -1.0
        }
        {
            Vector<float, 2> u({2., 1.});
            Vector<float, 2> v({4., 2.});
            std::cout << angle_cos(u, v) << std::endl;
            // 1.0
        }
        {
            Vector<float, 2> u({1., 2., 3.});
            Vector<float, 2> v({4., 5., 6.});
            std::cout << angle_cos(u, v) << std::endl;
            // 0.974631846
        }
    }
    {
        {
            Vector<double, 3> u({0., 0., 1.});
            Vector<double, 3> v({1., 0., 0.});
            std::cout << cross_product(u, v) << std::endl;
            // [0.]
            // [1.]
            // [0.]
        }
        {
            Vector<double, 3> u({1., 2., 3.});
            Vector<double, 3> v({4., 5., 6.});
            std::cout << cross_product(u, v) << std::endl;
            // [-3.]
            // [6.]
            // [-3.]
        }
        {
            Vector<double, 3> u({4., 2., -3.});
            Vector<double, 3> v({-2., -5., 16.});
            std::cout << cross_product(u, v) << std::endl;
            // [17.]
            // [-58.]
            // [-16.]
        }
    }
    {
        {
            Matrix<float, 2, 2> u({
                {1., 0.},
                {0., 1.},
            });
            Vector<float, 2> v({4., 2.});
            std::cout << u * v << std::endl;
            // [4.]
            // [2.]
        }
        {
            Matrix<float, 2, 2> u({
                {2., 0.},
                {0., 2.},
            });
            Vector<float, 2> v({4., 2.});
            std::cout << u * v << std::endl;
            // [8.]
            // [4.]
        }
        {
            Matrix<float, 2, 2> u({
                {2., -2.},
                {-2., 2.},
            });
            Vector<float, 2> v({4., 2.});
            std::cout << u * v << std::endl;
            // [4.]
            // [-4.]
        }
        {
            Matrix<float, 2, 2> u({
                {1., 0.},
                {0., 1.},
            });
            Matrix<float, 2, 2> v({
                {1., 0.},
                {0., 1.},
            });
            std::cout << u * v << std::endl;
            // [1., 0.]
            // [0., 1.]
        }
        {
            Matrix<float, 2, 2> u({
                {1., 0.},
                {0., 1.},
            });
            Matrix<float, 2, 2> v({
                {2., 1.},
                {4., 2.},
            });
            std::cout << u * v << std::endl;
            // [2., 1.]
            // [4., 2.]
        }
        {
            Matrix<float, 2, 2> u({
                {3., -5.},
                {6., 8.},
            });
            Matrix<float, 2, 2> v({
                {2., 1.},
                {4., 2.},
            });
            std::cout << u * v << std::endl;
            // [-14., -7.]
            // [44., 22.]
        }
    }
    {
        {
            Matrix<float, 2, 2> u({{1., 0.},
                                   {0., 1.}});
            std::cout << u.trace() << std::endl;
            // 2.0
        }
        {
            Matrix<float, 3, 3> u({{2., -5., 0.},
                                   {4., 3., 7.},
                                   {-2., 3., 4.}});
            std::cout << u.trace() << std::endl;
            // 9.0
        }
        {
            Matrix<float, 3, 3> u({{-2., -8., 4.},
                                   {1., -23., 4.},
                                   {0., 6., 4.}});
            std::cout << u.trace() << std::endl;
            // -21.0
        }
    }
    {
        {
            Matrix<float, 2, 2> u({{1., 0.},
                                   {0., 1.}});
            std::cout << u.transpose() << std::endl;
            // [1., 0.]
            // [0., 1.]
        }
        {
            Matrix<float, 2, 2> u({{1., 2.},
                                   {3., 4.}});
            std::cout << u.transpose() << std::endl;
            // [1., 3.]
            // [2., 4.]
        }
        {
            Matrix<float, 2, 3> u({{1., 2.},
                                   {3., 4.},
                                   {5., 6.}});
            std::cout << u.transpose() << std::endl;
            // [1., 3., 5.]
            // [2., 4., 6.]
        }
    }
    {
        {
            Matrix<float, 3, 3> u({{1., 0., 0.},
                                   {0., 1., 0.},
                                   {0., 0., 1.}});
            std::cout << u.row_echelon().first << std::endl;
            // [1.0, 0.0, 0.0]
            // [0.0, 1.0, 0.0]
            // [0.0, 0.0, 1.0]
        }
        {
            Matrix<float, 2, 2> u({
                {1., 2.},
                {3., 4.},
            });
            std::cout << u.row_echelon().first << std::endl;
            // [1.0, 0.0]
            // [0.0, 1.0]
        }
        {
            Matrix<float, 2, 2> u({
                {1., 2.},
                {2., 4.},
            });
            std::cout << u.row_echelon().first << std::endl;
            // [1.0, 0.0]
            // [2.0, 0.0]
        }
        {
            Matrix<float, 5, 3> u({{8., 5., -2., 4., 28.},
                                   {4., 2.5, 20., 4., -4.},
                                   {8., 5., 1., 4., 17.}});
            std::cout << u.row_echelon().first << std::endl;
            // [1., 0., 0., 0., 0.]
            // [0., 1., 0., 0., 0.]
            // [0., 0., 1., 0., 0.]
        }
    }
    {
        {
            Matrix<float, 2, 2> u({
                {1., -1.},
                {-1., 1.},
            });
            std::cout << u.determinant() << std::endl;
            // 0.0
        }
        {
            Matrix<float, 3, 3> u({
                {2., 0., 0.},
                {0., 2., 0.},
                {0., 0., 2.},
            });
            std::cout << u.determinant() << std::endl;
            // 8.0
        }
        {
            Matrix<float, 3, 3> u({
                {8., 5., -2.},
                {4., 7., 20.},
                {7., 6., 1.},
            });
            std::cout << u.determinant() << std::endl;
            // -174.0
        }
        {
            Matrix<float, 4, 4> u({
                {8., 5., -2., 4.},
                {4., 2.5, 20., 4.},
                {8., 5., 1., 4.},
                {28., -4., 17., 1.},
            });
            std::cout << u.determinant() << std::endl;
            // 1032
        }
    }
    {
        {
            Matrix<float, 3, 3> u({
                {1., 0., 0.},
                {0., 1., 0.},
                {0., 0., 1.},
            });
            std::cout << u.inverse() << std::endl;
            // [1.0, 0.0, 0.0]
            // [0.0, 1.0, 0.0]
            // [0.0, 0.0, 1.0]
        }
        {
            Matrix<float, 3, 3> u({
                {2., 0., 0.},
                {0., 2., 0.},
                {0., 0., 2.},
            });
            std::cout << u.inverse() << std::endl;
            // [0.5, 0.0, 0.0]
            // [0.0, 0.5, 0.0]
            // [0.0, 0.0, 0.5]
        }
        {
            Matrix<float, 3, 3> u({
                {8., 5., -2.},
                {4., 7., 20.},
                {7., 6., 1.},
            });
            std::cout << u.inverse() << std::endl;
            // [0.649425287, 0.097701149, -0.655172414]
            // [-0.781609195, -0.126436782, 0.965517241]
            // [0.143678161, 0.074712644, -0.206896552]
        }
    }
    {
        {
            Matrix<float, 3, 3> u({
                {1., 0., 0.},
                {0., 1., 0.},
                {0., 0., 1.},
            });
            std::cout << u.rank() << std::endl;
            // 3
        }
        {
            Matrix<float, 4, 3> u({
                {1., 2., 0., 0.},
                {2., 4., 0., 0.},
                {-1., 2., 1., 1.}
            });
            std::cout << u.rank() << std::endl;
            // 2
        }
        {
            Matrix<float, 3, 4> u({
                {8., 5., -2.},
                {4., 7., 20.},
                {7., 6., 1.},
                {21., 18., 7.},
            });
            std::cout << u.rank() << std::endl;
            // 3
        }
    }
    return 0;
}
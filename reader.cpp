#include <fstream>
#include <vector>
#include <cstddef>
#include <iostream>
#include <algorithm>

#include "array.h"
#include "param.h"

template<typename T>
void put(std::ostream &f, T value) {
    union {
        char buf[sizeof(T)];
        T val;
    } helper;
    helper.val = value;
    std::reverse(helper.buf, helper.buf + sizeof(T));
    f.write(helper.buf, sizeof(T));
}

int main() {
    std::fstream grid("grid", std::ios::in | std::ios::binary);

    array<3, float> x(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> y(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> z(range(-1, n+1), range(-1, n+1), 6);
    array<1, float> r(range(-2, nr+2));
    array<3, float> ox_r(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> oy_r(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> oz_r(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> s_r(n, n, 6);

    grid >> x >> y >> z >> r >> ox_r >> oy_r >> oz_r >> s_r;
    if (!!grid)
        std::cout << "Unread data left in grid" << std::endl;
    grid.close();

    std::fstream dat("120dat", std::ios::in | std::ios::binary);

    array<4, float> ro(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> p (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> s (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> u (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> v (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> w (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> hx(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> hy(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    array<4, float> hz(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);

    dat >> ro >> p >> s >> hx >> hy >> hz >> u >> v >> w;
    if (!!dat)
        std::cout << "Unread data left in dat" << std::endl;
    dat.close();

    for (int slice = 1; slice <= 6; slice ++) {
        std::ofstream vtk("res." + std::to_string(slice) + ".vtk", std::ios::out);
        vtk << "# vtk DataFile Version 3.0\n";
        vtk << "Cubed sphere grid " << slice << "\n";
        vtk << "BINARY\n";
        vtk << "DATASET STRUCTURED_GRID\n";
        vtk << "DIMENSIONS " << n + 1 << " " << n + 1 << " " << nr << "\n";
        vtk << "POINTS " << (n + 1) * (n + 1) * nr << " float\n";
        for (ptrdiff_t k = 1; k <= nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++)
                {
                    float rr = 0.5 * (r(k) + r(k - 1));
                    put<float>(vtk, rr * x(i, j, slice));
                    put<float>(vtk, rr * y(i, j, slice));
                    put<float>(vtk, rr * z(i, j, slice));
                }
        vtk << "POINT_DATA " << (n + 1) * (n + 1) * nr << "\n";
        vtk << "SCALARS rho float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++)
                    put<float>(vtk, ro(k, i, j, slice));
        vtk << "SCALARS p float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++)
                    put<float>(vtk, p(k, i, j, slice));
        vtk << "SCALARS s float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++)
                    put<float>(vtk, s(k, i, j, slice));
        vtk << "VECTORS v float\n";
        for (ptrdiff_t k = 1; k <= nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++) {
                    put<float>(vtk, u(k, i, j, slice));
                    put<float>(vtk, v(k, i, j, slice));
                    put<float>(vtk, w(k, i, j, slice));
                }
        vtk << "VECTORS H float\n";
        for (ptrdiff_t k = 1; k <= nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++) {
                    put<float>(vtk, hx(k, i, j, slice));
                    put<float>(vtk, hy(k, i, j, slice));
                    put<float>(vtk, hz(k, i, j, slice));
                }
        vtk.close();
    }
    return 0;
}

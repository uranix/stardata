#include <fstream>
#include <vector>
#include <cstddef>
#include <iostream>
#include <algorithm>

#include <map>

#include "array.h"
#include "param.h"

using namespace std;
using fort::range;

template<typename T>
void put(ostream &f, T value) {
    union {
        char buf[sizeof(T)];
        T val;
    } helper;
    helper.val = value;
    reverse(helper.buf, helper.buf + sizeof(T));
    f.write(helper.buf, sizeof(T));
}

int main() {
    fstream grid("grid", ios::in | ios::binary);

    fort::array<3, float> x(range(-1, n+1), range(-1, n+1), 6);
    fort::array<3, float> y(range(-1, n+1), range(-1, n+1), 6);
    fort::array<3, float> z(range(-1, n+1), range(-1, n+1), 6);
    fort::array<1, float> r(range(-2, nr+2));

    grid >> x >> y >> z >> r;

    cout << "Rmin = " << r(0) << endl;
    cout << "Rmax = " << r(nr) << endl;

    if (!!grid)
        cout << "Unread data left in grid" << endl;
    grid.close();

    fstream dat("120dat", ios::in | ios::binary);

    fort::array<4, float> ro(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> p (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> s (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> u (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> v (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> w (range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> hx(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> hy(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);
    fort::array<4, float> hz(range(0, nr+1), range(-1, n+2), range(-1, n+2), 6);

    dat >> ro >> p >> s >> hx >> hy >> hz >> u >> v >> w;
    if (!!dat)
        cout << "Unread data left in dat" << endl;
    dat.close();

    const int Nr = 1;

    for (int side = 1; side <= 6; side++) {
        ofstream vtk("side." + to_string(side) + ".vtk");

        vtk << "# vtk DataFile Version 3.0\n";
        vtk << "Cubed sphere grid, side " << side << "\n";
        vtk << "BINARY\n";
        vtk << "DATASET STRUCTURED_GRID\n";
        vtk << "DIMENSIONS " << n + 1 << " " << n + 1 << " " << Nr + 1 << "\n";
        vtk << "POINTS " <<  (Nr + 1) * (n + 1) * (n + 1) << " float\n";
        for (ptrdiff_t k = 0; k <= Nr; k++)
            for (ptrdiff_t j = 0; j <= n; j++)
                for (ptrdiff_t i = 0; i <= n; i++) {
                    float rr = r(k);
                    put<float>(vtk, rr * x(i, j, side));
                    put<float>(vtk, rr * y(i, j, side));
                    put<float>(vtk, rr * z(i, j, side));
                }

        const int ncells = n * n * Nr;
        vtk << "CELL_DATA " << ncells << endl;

        vtk << "SCALARS rho float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= n; j++)
                for (ptrdiff_t i = 1; i <= n; i++)
                    put<float>(vtk, ro(k, i, j, side));
        vtk << "VECTORS v float\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= n; j++)
                for (ptrdiff_t i = 1; i <= n; i++) {
                    put<float>(vtk, u(k, i, j, side));
                    put<float>(vtk, v(k, i, j, side));
                    put<float>(vtk, w(k, i, j, side));
                }

        vtk.close();
    }

    return 0;
}

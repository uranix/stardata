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

struct ind {
    int s;
    ptrdiff_t i, j, k;
};

ind probe(float x, float y, float z, const fort::array<1, float> &r, const int nr, const int n) {
    float yz = std::max(std::abs(y), std::abs(z));
    float xz = std::max(std::abs(x), std::abs(z));
    float xy = std::max(std::abs(x), std::abs(y));

    float rr = sqrt(x * x + y * y + z * z);

    int k = std::lower_bound(&r(0), &r(nr) + 1, rr) - &r(0);
    if (k < 1)
        k = 1;
    if (k > nr)
        k = nr;

    int s;
    float phi, psi;
    if (-x >= yz) {
        s = 1;
        phi = std::atan2(-y, -x);
        psi = std::atan2(+z, -x);
    } else if (+z >= xy) {
        s = 2;
        phi = std::atan2(-y, +z);
        psi = std::atan2(+x, +z);
    } else if (+y >= xz) {
        s = 3;
        phi = std::atan2(-x, +y);
        psi = std::atan2(+z, +y);
    } else if (-y >= xz) {
        s = 4;
        phi = std::atan2(+x, -y);
        psi = std::atan2(+z, -y);
    } else if (+x >= yz) {
        s = 5;
        phi = std::atan2(+y, +x);
        psi = std::atan2(+z, +x);
    } else { /* -z >= xy */
        s = 6;
        phi = std::atan2(-y, -z);
        psi = std::atan2(-x, -z);
    }

    float pi2 = 2 * std::atan(1.);
    phi /= pi2;
    psi /= pi2;
    phi += .5;
    psi += .5;
    phi *= n;
    psi *= n;

    int i, j;
    i = phi + 1;
    j = psi + 1;
    if (i < 1)
        i = 1;
    if (i > n)
        i = n;
    if (j < 1)
        j = 1;
    if (j > n)
        j = n;

    return ind{s, i, j, k};
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

    const int Nr = nr;

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

        vtk << "SCALARS del float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= n; j++)
                for (ptrdiff_t i = 1; i <= n; i++) {
                    float rr = 0.5 * (r(k-1) + r(k));
                    float xx = 0.5 * (x(i-1, j-1, side) + x(i, j, side));
                    float yy = 0.5 * (y(i-1, j-1, side) + y(i, j, side));
                    float zz = 0.5 * (z(i-1, j-1, side) + z(i, j, side));
                    const auto &ijk = probe(rr * xx, rr * yy, rr * zz, r, nr, n);
                    float del =
                        std::abs(ijk.s - side) +
                        std::abs(ijk.i - i) +
                        std::abs(ijk.j - j) +
                        std::abs(ijk.k - k);
                    put<float>(vtk, del);
                }
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

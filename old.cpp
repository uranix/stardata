#include <fstream>
#include <vector>
#include <cstddef>
#include <iostream>
#include <algorithm>

#include <map>

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

struct idx {
    int side;
    int i, j;
    idx(int side, int i, int j) : side(side), i(i), j(j) { }
    bool operator<(const idx &o) {
        if (side < o.side)
            return true;
        if (side > o.side)
            return false;
        if (i < o.i)
            return true;
        if (i > o.i)
            return false;
        return j < o.j;
    }
};

struct face {
    int a, b, c, d;
};

std::vector<face> mapping{
    face{-1, -1, -1, -1},
    face{0, 1, 2, 3},
    face{3, 2, 7, 6},
    face{5, 0, 3, 6},
    face{1, 4, 7, 2},
    face{4, 5, 6, 7},
    face{5, 4, 1, 0}};

size_t edge(int a, int b, int s, int n) {
    static std::map<std::pair<int, int>, int> v2e;
    if (v2e.empty()) {
        v2e[std::make_pair(0, 1)] = 0;
        v2e[std::make_pair(1, 2)] = 1;
        v2e[std::make_pair(2, 3)] = 2;
        v2e[std::make_pair(0, 3)] = 3;

        v2e[std::make_pair(4, 5)] = 4;
        v2e[std::make_pair(5, 6)] = 5;
        v2e[std::make_pair(6, 7)] = 6;
        v2e[std::make_pair(4, 7)] = 7;

        v2e[std::make_pair(2, 7)] = 8;
        v2e[std::make_pair(1, 4)] = 9;
        v2e[std::make_pair(0, 5)] = 10;
        v2e[std::make_pair(3, 6)] = 11;
    }

    if (a > b)
        return edge(b, a, n - 2 - s, n);

    int eno = v2e[std::make_pair(a, b)];

    return eno * (n - 1) + s;
}

size_t translate(const idx &p, const int n) {
    const face &f = mapping[p.side];
    bool imid = p.i > 0 && p.i < n;
    bool jmid = p.j > 0 && p.j < n;
    if (!imid && !jmid) {
        if (p.i == 0 && p.j == 0)
            return f.a;
        if (p.i == 0 && p.j == n)
            return f.b;
        if (p.i == n && p.j == n)
            return f.c;
        if (p.i == n && p.j == 0)
            return f.d;
    }
    if (p.i == 0)
        return 8 + edge(f.d, f.a, n - 1 - p.j, n);
    if (p.j == 0)
        return 8 + edge(f.a, f.b, p.i - 1, n);
    if (p.i == n)
        return 8 + edge(f.b, f.c, p.j - 1, n);
    if (p.j == n)
        return 8 + edge(f.c, f.d, n - 1 - p.i, n);

    if (imid && jmid)
        return
            8 + 12 * (n - 1) +
            (p.side - 1) * (n - 1) * (n - 1) +
            + (p.j - 1) * (n - 1) + (p.i - 1);
}

int main() {
    std::fstream grid("grid", std::ios::in | std::ios::binary);

    array<3, float> x(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> y(range(-1, n+1), range(-1, n+1), 6);
    array<3, float> z(range(-1, n+1), range(-1, n+1), 6);
    array<1, float> r(range(-2, nr+2));

    grid >> x >> y >> z >> r;
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

    const int N = n;
    const int Nr = nr;

    for (int side = 1; side <= 6; side ++) {
        std::ofstream vtk("res." + std::to_string(side) + ".vtk", std::ios::out);
        vtk << "# vtk DataFile Version 3.0\n";
        vtk << "Cubed sphere grid " << side << "\n";
        vtk << "BINARY\n";
        vtk << "DATASET STRUCTURED_GRID\n";
        vtk << "DIMENSIONS " << N + 1 << " " << N + 1 << " " << Nr + 1 << "\n";
        vtk << "POINTS " << (N + 1) * (N + 1) * (Nr + 1) << " float\n";
        for (ptrdiff_t k = 0; k <= Nr; k++)
            for (ptrdiff_t j = 0; j <= N; j++)
                for (ptrdiff_t i = 0; i <= N; i++)
                {
                    float rr = r(k);
                    put<float>(vtk, rr * x(i, j, side));
                    put<float>(vtk, rr * y(i, j, side));
                    put<float>(vtk, rr * z(i, j, side));
                }
        vtk << "CELL_DATA " << N * N * Nr << "\n";
        vtk << "SCALARS ridx float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++)
                    put<float>(vtk, k);
        vtk << "SCALARS iidx float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++)
                    put<float>(vtk, i);
        vtk << "SCALARS jidx float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++)
                    put<float>(vtk, j);
        vtk << "SCALARS rho float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++)
                    put<float>(vtk, ro(k, i, j, side));
        vtk << "SCALARS p float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++)
                    put<float>(vtk, p(k, i, j, side));
        vtk << "SCALARS s float\nLOOKUP_TABLE default\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++)
                    put<float>(vtk, s(k, i, j, side));
        vtk << "VECTORS v float\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++) {
                    put<float>(vtk, u(k, i, j, side));
                    put<float>(vtk, v(k, i, j, side));
                    put<float>(vtk, w(k, i, j, side));
                }
        vtk << "VECTORS H float\n";
        for (ptrdiff_t k = 1; k <= Nr; k++)
            for (ptrdiff_t j = 1; j <= N; j++)
                for (ptrdiff_t i = 1; i <= N; i++) {
                    put<float>(vtk, hx(k, i, j, side));
                    put<float>(vtk, hy(k, i, j, side));
                    put<float>(vtk, hz(k, i, j, side));
                }
        vtk.close();
    }
    return 0;
}

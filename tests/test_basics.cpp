/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "../src/util/elementaries.h"


#include "effedup.hpp"

#include <utility>
#include <vector>

static
auto testSignTemplated()
{
    std::vector<std::pair<int, int>> values{
        { -2, -1 },
        { -1, -1 },
        {  0,  0 },
        {  1,  1 },
        {  2,  1 }
    };

    for (auto [v, s] : values) {
        auto as = LLKAInternal::sign(v);
        EFF_expect(as, s, "Templated sign() function returned incorrect value");
    }
}


static
auto testSignDouble()
{
    std::vector<std::pair<double, double>> values{
        { -2.3, -1.0 },
        { -1.4, -1.0 },
        { -0.4, -1.0 },
        {  -0,  0 },
        {   0,  0 },
        { 0.7,  1 },
        { 1.6,  1 },
        { 2.0,  1 }
    };

    for (auto [v, s] : values) {
        auto as = LLKAInternal::sign(v);
        EFF_expect(as, s, "Double sign() function returned incorrect value");
    }
}

auto main() -> int
{
    testSignTemplated();
    testSignDouble();
}


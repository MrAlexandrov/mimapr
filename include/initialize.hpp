#pragma once
#include "options.hpp"
#include <vector>
#include <algorithm>

namespace NOptions {

constexpr long double A = 1, B = 0, C = -1, D = -10;

const std::vector<TRestriction> restrictions = {
    {ERestrictionGrade::SECOND, 2, 10},
    // {ERestrictionGrade::FIRST, 8, 5},
    {ERestrictionGrade::THIRD, 8, 1}
};

const long double minimum = std::min(static_cast<long double>(2), std::min_element(restrictions.begin(), restrictions.end())->Position);
const long double maximum = std::min(static_cast<long double>(8), std::max_element(restrictions.begin(), restrictions.end())->Position);

} // NOptions
#pragma once
#include "options.hpp"

namespace NOptions {

constexpr long double A = 1, B = 0, C = -1, D = -10;
constexpr TRestriction LOWER = {ERestrictionGrade::SECOND, 2, 10};
constexpr TRestriction UPPER = {ERestrictionGrade::FIRST, 8, 5};

} // NOptions
#include "Demography.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

TEST_CASE("Demography constructor exceptions") {
  SECTION("Throws if times and sizes vector lengths differ") {
    CHECK_THROWS_WITH(Demography({0, 0, 0}, {0}), Catch::Matchers::ContainsSubstring("Demography times and sizes must have equal length"));
    CHECK_THROWS_WITH(Demography({0}, {0, 0, 0}), Catch::Matchers::ContainsSubstring("Demography times and sizes must have equal length"));
  }

  SECTION("Throws if times not non-negative") {
    CHECK_THROWS_WITH(Demography({1}, {-1}), Catch::Matchers::ContainsSubstring("Demography expects non-negative times"));
  }

  SECTION("Throws if sizes not positive") {
    CHECK_THROWS_WITH(Demography({0}, {0}), Catch::Matchers::ContainsSubstring("strictly positive sizes"));
    CHECK_THROWS_WITH(Demography({-1}, {0}), Catch::Matchers::ContainsSubstring("strictly positive sizes"));
  }

  SECTION("Throws if times not increasing") {
    CHECK_THROWS_WITH(Demography({1, 1, 1}, {3, 2, 1}), Catch::Matchers::ContainsSubstring("Demography times must be strictly increasing. Found 2 after 3 at index 1"));
  }

  SECTION("Throws if times does not start at 0") {
    CHECK_THROWS_WITH(Demography({1}, {5}), Catch::Matchers::ContainsSubstring("Demography must start at time 0.0"));
  }
}

TEST_CASE("Demography std_to_gen exceptions") {
  SECTION("std_to_gen value must be non-negative") {
    Demography d{{1}, {0}};
    CHECK_THROWS_WITH(d.std_to_gen(-1), Catch::Matchers::ContainsSubstring("Demography can only convert non-negative times to std"));
  }
}

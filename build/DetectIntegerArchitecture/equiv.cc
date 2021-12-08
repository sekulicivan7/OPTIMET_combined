#include <cstdint>
#include <type_traits>
static_assert(
  std::is_same<unsigned long, std::int8_t>::value
  or std::is_same<unsigned long, std::int16_t>::value
  or std::is_same<unsigned long, std::int32_t>::value
  or std::is_same<unsigned long, std::int64_t>::value
  or std::is_same<unsigned long, std::uint8_t>::value
  or std::is_same<unsigned long, std::uint16_t>::value
  or std::is_same<unsigned long, std::uint32_t>::value
  or std::is_same<unsigned long, std::uint64_t>::value
  or std::is_same<unsigned long, unsigned int>::value
,
  "unsigned long vs std::int8_t;std::int16_t;std::int32_t;std::int64_t;std::uint8_t;std::uint16_t;std::uint32_t;std::uint64_t;unsigned int");
int main(int nargs, char **argv) {
  return 0;
}

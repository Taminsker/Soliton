#ifndef SRC_ENUMS_COMMON_HEAD_HPP
#define SRC_ENUMS_COMMON_HEAD_HPP

#include <iostream>
#include <string>
#include "../solitonheader.hpp"

template <typename T>
T from_string (const std::string & s);

template <typename T>
std::string to_string (const T & type);

template <typename T, typename U>
U Convert (const T & type);

#endif /* SRC_ENUMS_COMMON_HEAD_HPP */

#ifndef COMMON_HEAD_H
#define COMMON_HEAD_H

#include <iostream>
#include <string>

template<typename T>
T from_string (const std::string& s);

template<typename T>
std::string to_string (const T& type);

template <typename T, typename U>
U Convert(const T& type);

#endif // COMMON_HEAD_H

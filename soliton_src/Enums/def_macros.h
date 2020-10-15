#ifndef DEF_MACROS_H
#define DEF_MACROS_H

#define TO_STRING_SPECIALIZATION(X)         \
    template<>                              \
    std::string to_string (const X& type)   \
{                                           \
    typedef X _ThisType;                    \
    switch (type)

#define CASE(X)                             \
    case _ThisType::X:                      \
        return #X

#define END_CASE_DEFAULT(X)                 \
    default:                                \
        return #X;                          \
    } (void) type

#define FROM_STRING_SPECIALIZATION(X)       \
    template<>                              \
    X from_string (const std::string& s)    \
{                                           \
    typedef X _ThisType;

#define SIMPLE_IF(X)                        \
    if (s == #X)                            \
        return _ThisType::X

#define ELSEIF(X)                           \
    else if (s == #X)                       \
        return _ThisType::X

#define ELSE(X)                             \
    else                                    \
        return _ThisType::X;                \
    } (void)s

#define CONVERT_SPECIALIZATION(X, Y)        \
    template<>                              \
    Y Convert (const X& type)               \
{                                           \
    typedef X _def_Type;                    \
    typedef Y _after_Type;                  \
    switch (type)


#define CASE_CONVERT(X, Y)                  \
    case _def_Type::X:                      \
        return _after_Type::Y

#define END_CASE_DEFAULT_CONVERT(X)         \
    default:                                \
        return _after_Type::X;              \
    } (void) type

#endif // DEF_MACROS_H

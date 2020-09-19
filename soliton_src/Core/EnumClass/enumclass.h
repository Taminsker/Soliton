#ifndef ENUMCLASS_H
#define ENUMCLASS_H

/**
 * from https://stackoverflow.com/a/8498694
 */
template< typename T >
class EnumClass
{
public:
    class Iterator
    {
    public:
        Iterator (int value) :
            m_value( value )
        {}

        T operator*() const
        {
            return static_cast<T>(m_value);
        }

        void operator++()
        {
            ++m_value;
            return;
        }

        bool operator!=(Iterator it)
        {
            return m_value != it.m_value;
        }

    private:
        int m_value;
    };

};

template <typename T>
typename EnumClass <T>::Iterator begin (EnumClass<T>)
{
    return typename EnumClass <T>::Iterator (static_cast<int>(T::FIRST));
}

template <typename T>
typename EnumClass <T>::Iterator end (EnumClass<T>)
{
    return typename EnumClass <T>::Iterator (static_cast<int>(T::LAST) + 1);
}


#endif // ENUMCLASS_H

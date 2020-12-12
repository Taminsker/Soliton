#ifndef SRC_ALGORITHMS_MATH_ITERATOR_ITERATORS_HPP
#define SRC_ALGORITHMS_MATH_ITERATOR_ITERATORS_HPP

#include <iterator>

#include "../../../solitonheader.hpp"

namespace internal
{
// inherit from this
class IterativeObject
{
public:
    typedef IterativeObject SuperClass;

    IterativeObject ();
    IterativeObject (const IterativeObject & tocopy);
    virtual ~IterativeObject ();

    virtual IterativeObject * operator+ (ul_t i) = 0;
    virtual void              operator++ ()      = 0;
};

// inherit from this
template <typename T>
class IterativeContainer
{
public:
    class iterator
    {
    public:
        typedef iterator                  self_type;
        typedef T                         value_type;
        typedef T &                       reference;
        typedef T *                       pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef ul_t                      difference_type;

        iterator (pointer ptr) : m_ptr (ptr)
        {
        }

        self_type
        operator++ ()
        {
            self_type i = *this;
            m_ptr->   operator++ ();
            return i;
        }

        self_type
        operator++ (int)
        {
            m_ptr->operator++ ();
            return *this;
        }

        reference
        operator* ()
        {
            return *m_ptr;
        }

        pointer
        operator-> ()
        {
            return m_ptr;
        }

        bool
        operator== (const self_type & rhs)
        {
            return m_ptr->operator== (*rhs.m_ptr);
        }

        bool
        operator!= (const self_type & rhs)
        {
            return m_ptr->operator!= (*rhs.m_ptr);
        }

    private:
        pointer m_ptr;
    };

    typedef IterativeContainer<T> SuperClass;

    IterativeContainer ()
    {
    }

    IterativeContainer (const IterativeContainer & tocopy) : m_data (new T (*tocopy.m_data)),
                                                             m_size (tocopy.m_size)
    {
    }

    virtual ~IterativeContainer ()
    {
        delete m_data;
    }

    ul_t
    size () const
    {
        return m_size;
    }

    iterator
    begin ()
    {
        return iterator (m_data);
    }

    iterator
    end ()
    {
        T copy = *m_data;
        return iterator (copy + (m_size));
    }

protected:
    T *  m_data;
    ul_t m_size;
};
}  // namespace internal

#endif /* SRC_ALGORITHMS_MATH_ITERATOR_ITERATORS_HPP */

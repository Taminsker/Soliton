#ifndef HETCONTAINER_H
#define HETCONTAINER_H

#include "../Point/point.h"
#include <ProgDef/proddef.h>

#include <vector>
#include <string>

#define SPACEARRAY std::left << std::setw (25)

class Cell;
class Edge;
class Point;

template <typename T>
class HetContainer
{
public:
    typedef struct
    {
        std::string name;
        std::vector <T> vec;
    } Array;

    typedef std::vector <Array*> Arrays;

    HetContainer () : m_loa(new Arrays()) {}
    ~HetContainer ()
    {
        for (Array* a : *m_loa)
            delete a;

        delete m_loa;
    }

    Arrays* All ()
    {
        return m_loa;
    }

    Array* Get (int num)
    {
        if (num >= 0 && num < static_cast<int>(m_loa->size ()))
            return m_loa->at (ul_t (num));
        return nullptr;
    }

    Array* Get (std::string name)
    {
        for (Array* arr : *m_loa)
            if (arr->name == name)
                return arr;

        return nullptr;
    }

    void Add (std::string name, std::vector <T>& vec)
    {
        Array* arr = new Array();
        arr->name = name;
        arr->vec = vec;

        m_loa->push_back (arr);

        return;
    }

    void Remove (int num)
    {

        delete m_loa->at (ul_t (num));
        m_loa->erase (m_loa->begin () + num);
    }

    void Remove (std::string name)
    {
        for (ul_t i = 0; i < m_loa->size (); ++i)
        {
            if (m_loa->at (i)->name == name)
            {
                delete m_loa->at (i);
                m_loa->erase (m_loa->begin () + i);
                return;
            }
        }

        return;
    }

    int GetSize ()
    {
        return static_cast<int>(m_loa->size ());
    }

    void Print (std::ostream& out = std::cout)
    {
        for (Array* a : *m_loa)
        {
            out << " * " << COLOR_YELLOW << "\t\"" << a->name << "\"" << COLOR_DEFAULT << ", " << std::flush;
            if (std::is_same<T, std::string>::value)      out << "std::string";
            if (std::is_same<T, bool>::value)          out << "bool";
            if (std::is_same<T, int>::value)          out << "int";
            if (std::is_same<T, real_t>::value)         out << "real_t";

            out << ", " << a->vec.size() << ENDLINE;
        }
    }

    void Delete (int i)
    {
        for (Array* array : *m_loa)
        {
            array->vec.erase (array->vec.begin () + i);
        }
    }

    void Add (int i)
    {
        for (Array* array : *m_loa)
        {
            array->vec.insert (ul_t (i)) = T();
        }
    }

private:
    Arrays* m_loa;
};

template <typename T>
class HetContainer<T*>
{
public:
    typedef struct
    {
        std::string name;
        std::vector <T*> vec;
    } Array;

    typedef std::vector <Array*> Arrays;

    HetContainer () : m_loa(new Arrays()) {}
    ~HetContainer ()
    {
        for (Array* a : *m_loa)
        {
            for (Point* p : a->vec)
                delete p;
            delete a;
        }
        delete m_loa;
    }

    Arrays* All ()
    {
        return m_loa;
    }

    Array* Get (int num)
    {
        if (num >= 0 && num < static_cast<int>(m_loa->size ()))
            return m_loa->at (ul_t (num));
        return nullptr;
    }

    Array* Get (std::string name)
    {
        for (Array* arr : *m_loa)
            if (arr->name == name)
                return arr;

        return nullptr;
    }

    void Add (std::string name, std::vector <T*>& vec)
    {
        Array* arr = new Array();
        arr->name = name;
        arr->vec = vec;
        m_loa->push_back (arr);

        return;
    }

    void Remove (int num)
    {
        for (T* p : m_loa->at (ul_t (num))->vec)
            delete p;

        delete m_loa->at (ul_t (num));
        m_loa->erase (m_loa->begin () + num);
    }

    void Remove (std::string name)
    {
        for (ul_t i = 0; i < m_loa->size (); ++i)
        {
            if (m_loa->at (i)->name == name)
            {
                for (Point* p : m_loa->at (i)->vec)
                    delete p;

                delete m_loa->at (i);
                m_loa->erase (m_loa->begin () + static_cast<int>(i));

                return;
            }
        }

        return;
    }

    int GetSize ()
    {
        return static_cast<int>(m_loa->size ());
    }

    void Print (std::ostream& out = std::cout)
    {
        for (Array* a : *m_loa)
        {
            out << " * " << COLOR_YELLOW << "\t\"" << a->name << "\"" << COLOR_DEFAULT << ", " << std::flush;
            if (std::is_same<T, Point>::value)
                out << "Point*";
            out << ", " << a->vec.size() << ENDLINE;
        }
    }

    void Delete (int i)
    {
        for (Array* array : *m_loa)
        {
            delete array->vec.at (ul_t (i));
            array->vec.erase (array->vec.begin () + i);
        }

        return;
    }

    void Add (int i)
    {
        for (Array* array : *m_loa)
            array->vec.insert (ul_t (i)) = new T ();

        return;
    }

private:
    Arrays* m_loa;
};


template<typename T>
struct identity {typedef T type;};

template <typename T>
class DataContainer
{
public:
    DataContainer (std::vector <T>* vector) :
        m_link   (vector),
        m_bool   (new HetContainer<bool> ()),
        m_real_t  (new HetContainer<real_t> ()),
        m_point   (new HetContainer<Point*> ()),
        m_int    (new HetContainer<int> ()),
        m_string  (new HetContainer<std::string> ())
    {}

    ~DataContainer ()
    {
        delete m_bool;
        delete m_real_t;
        delete m_int;
        delete m_string;
        delete m_point;
    }

    template<typename U>
    HetContainer<U>* GetAll ()
    {
        return GetByIdentity (identity<U>());
    }

    template<typename U>
    typename HetContainer<U>::Array* Get (int num)
    {
        return GetAll<U> ()->Get(num);
    }

    template<typename U>
    typename HetContainer<U>::Array* Get (std::string name)
    {
        return GetAll<U> ()->Get(name);
    }

    template<typename U>
    void Add (std::string name, std::vector <U>& vec)
    {
        return GetAll<U> ()->Add (name, vec);
    }

    template<typename U>
    void Remove (int num)
    {
        return GetAll<U> ()->Remove (num);
    }

    template<typename U>
    void Remove (std::string name)
    {
        return GetAll<U> ()->Remove(name);
    }

    int GetNumberOfArrays ()
    {
        int count = 0;
        count += m_bool->GetSize ();
        count += m_real_t->GetSize ();
        count += m_point->GetSize ();
        count += m_int->GetSize ();
        count += m_string->GetSize ();

        return count;
    }

    void Print (std::ostream& out = std::cout)
    {
        INFOS << "Type of Data : \t" << std::flush;
        if (std::is_same<T, Cell*>::value)
            std::cout << "Cell*" << ENDLINE;
        else if (std::is_same<T, Point*>::value)
            std::cout << "Point*" << ENDLINE;
        else if (std::is_same<T, Edge*>::value)
            std::cout << "Edge*" << ENDLINE;
        else
            std::cout << COLOR_RED << "unknow-type" << ENDLINE;

        INFOS << "Number of arrays : \t" << GetNumberOfArrays () << " " << ENDLINE;

        m_bool->Print (out);
        m_real_t->Print (out);
        m_point->Print (out);
        m_int->Print (out);
        m_string->Print (out);
    }

    void Delete (int i)
    {
        m_bool->Delete (i);
        m_real_t->Delete (i);
        m_point->Delete (i);
        m_int->Delete (i);
        m_string->Delete (i);

        return;
    }

    void Add (int i)
    {
        m_bool->Add (i);
        m_real_t->Add (i);
        m_point->Add (i);
        m_int->Add (i);
        m_string->Add (i);

        return;
    }

private:
    std::vector <T>*      m_link;
    HetContainer<bool>*     m_bool;
    HetContainer<real_t>*    m_real_t;
    HetContainer<Point*>*    m_point;
    HetContainer<int>*     m_int;
    HetContainer<std::string>* m_string;

    template<typename U>
    HetContainer<U>* GetByIdentity (identity<U>)
    {
        return nullptr;
    }

    HetContainer<bool>* GetByIdentity (identity<bool>)
    {
        return m_bool;
    }

    HetContainer<real_t>* GetByIdentity (identity<real_t>)
    {
        return m_real_t;
    }

    HetContainer<Point*>* GetByIdentity (identity<Point*>)
    {
        return m_point;
    }

    HetContainer<int>* GetByIdentity (identity<int>)
    {
        return m_int;
    }

    HetContainer<std::string>* GetByIdentity (identity<std::string>)
    {
        return m_string;
    }

};

template <typename T>
std::ostream & operator<< (std::ostream &out, const std::vector<T> vec)
{
    for (ul_t i = 0; i < vec.size (); ++i)
        out << vec.at (i) << std::endl;
    return out;
}


#endif // HETCONTAINER_H

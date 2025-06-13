#ifndef POLY_H
#define POLY_H

#include <memory>
#include <ranges>
#include <iostream>
#include <complex>

template<class __T>
class __iterator{
public:
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef __T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef __T* pointer;
    typedef __T& reference;
private:
    pointer __ptr;
public:
    constexpr __iterator(): __ptr(nullptr){}
    constexpr __iterator(pointer _ptr): __ptr(_ptr){}
    constexpr __iterator(const __iterator& other): __ptr(other.__ptr){}
    constexpr __iterator& operator=(const __iterator& other){if(this != &other) __ptr = other.__ptr; return *this;}
    constexpr reference operator*() const{return *__ptr;}
    constexpr __iterator& operator++(){++__ptr; return *this;}
    constexpr __iterator operator++(int){auto __tmp = *this; ++__ptr; return __tmp;}
    constexpr __iterator& operator--(){--__ptr; return *this;}
    constexpr __iterator operator--(int){auto __tmp = *this; --__ptr; return __tmp;}
    constexpr bool operator==(const __iterator& other) const{return __ptr == other.__ptr;}
};

template<class __T>
class __riterator{
public:
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef __T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef __T* pointer;
    typedef __T& reference;
private:
    pointer __ptr;
public:
    constexpr __riterator(): __ptr(nullptr){}
    constexpr __riterator(pointer _ptr): __ptr(_ptr){}
    constexpr __riterator(const __riterator& other): __ptr(other.__ptr){}
    constexpr __riterator& operator=(const __riterator& other){if(this != &other) __ptr = other.__ptr; return *this;}
    constexpr reference operator*() const{return *__ptr;}
    constexpr __riterator& operator++(){--__ptr; return *this;}
    constexpr __riterator operator++(int){auto __tmp = *this; --__ptr; return __tmp;}
    constexpr __riterator& operator--(){++__ptr; return *this;}
    constexpr __riterator operator--(int){auto __tmp = *this; ++__ptr; return __tmp;}
    constexpr bool operator==(const __riterator& other) const{return __ptr == other.__ptr;}
};



/**
 * @brief A polynomial class with math operators. Features specializations for constants and conversions between complex types.
 * It is always constructible via aggregates.
 * @tparam __T Type of element. Must be a complete type.
 * @tparam __N The degree of the polynomial.
 */
template<class __T, std::size_t __N>
class poly{
public:
    typedef __T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef __iterator<__T> iterator;
    typedef __iterator<__T> const_iterator;
    typedef __riterator<__T> reverse_iterator;
    typedef __riterator<__T> const_reverse_iterator;

    __T __c[__N + 1];

    constexpr void
    fill(const value_type& __u)
    { std::fill_n(begin(), size(), __u); }

    constexpr void
    swap(poly& __other)
    noexcept(std::is_nothrow_swappable_v<__T>)
    { std::swap_ranges(begin(), end(), __other.begin()); }

    [[nodiscard, gnu::const]]
    constexpr iterator
    begin() noexcept
    { return iterator(data()); }

    [[nodiscard]]
    constexpr const_iterator
    begin() const noexcept
    { return const_iterator(data()); }

    [[nodiscard, gnu::const]]
    constexpr iterator
    end() noexcept
    { return iterator(data() + __N + 1); }

    [[nodiscard]]
    constexpr const_iterator
    end() const noexcept
    { return const_iterator(data() + __N + 1); }

    [[nodiscard, gnu::const]]
    constexpr reverse_iterator
    rbegin() noexcept
    { return reverse_iterator(data() + __N); }

    [[nodiscard]]
    constexpr const_reverse_iterator
    rbegin() const noexcept
    { return const_reverse_iterator(data() + __N); }

    [[nodiscard, gnu::const]]
    constexpr reverse_iterator
    rend() noexcept
    { return reverse_iterator(data() - 1); }

    [[nodiscard]]
    constexpr const_reverse_iterator
    rend() const noexcept
    { return const_reverse_iterator(data() - 1); }

    [[nodiscard]]
    constexpr const_iterator
    cbegin() const noexcept
    { return const_iterator(data()); }

    [[nodiscard]]
    constexpr const_iterator
    cend() const noexcept
    { return const_iterator(data() + __N + 1); }

    [[nodiscard]]
    constexpr const_reverse_iterator
    crbegin() const noexcept
    { return const_reverse_iterator(data() + __N); }

    [[nodiscard]]
    constexpr const_reverse_iterator
    crend() const noexcept
    { return const_reverse_iterator(data() - 1); }

    [[nodiscard, gnu::const]]
    constexpr size_type
    size() const noexcept { return __N + 1; }

    [[nodiscard, gnu::const]]
    constexpr size_type
    max_size() const noexcept { return __N + 1; }

    [[nodiscard, gnu::const]]
    constexpr bool
    empty() const noexcept { return false; }

    [[nodiscard]]
    constexpr reference
    operator[](size_type __n) noexcept
    { return __c[__n]; }

    [[nodiscard]]
    constexpr const_reference 
    operator[](size_type __n) const noexcept
    { return __c[__n]; }

    constexpr reference
    at(size_type __n)
    { if(__n >= __N + 1) 
        std::__throw_out_of_range_fmt("poly::at: __n (which is %zu) >= __N + 1 (which is %zu)", __n, __N + 1);
    return __c[__n]; }

    constexpr const_reference
    at(size_type __n) const
    { return __n < __N + 1 ? __c[__n]
    : ( std::__throw_out_of_range_fmt("poly::at: __n (which is %zu) >= __N + 1 (which is %zu)", __n, __N + 1),
        __c[__n]); }

    [[nodiscard]]
    constexpr reference
    front() noexcept
    { return __c[0]; }

    [[nodiscard]]
    constexpr const_reference
    front() const noexcept
    { return __c[0]; }
    
    [[nodiscard]]
    constexpr reference
    back() noexcept
    { return __c[__N]; }

    [[nodiscard]]
    constexpr const_reference
    back() const noexcept
    { return __c[__N]; }

    [[nodiscard, gnu::const]]
    constexpr pointer
    data() noexcept
    { return static_cast<pointer>(__c); }

    [[nodiscard]]
    constexpr const_pointer
    data() const noexcept
    { return static_cast<const_pointer>(__c); }

    [[nodiscard, gnu::const]]
    constexpr size_type
    degree() const noexcept
    { return __N; }

    [[nodiscard]]
    constexpr bool
    operator==(const poly& __other) const noexcept
    { return std::equal(begin(), end(), __other.begin()); }

    [[nodiscard]]
    constexpr std::compare_three_way_result_t<__T>
    operator<=>(const poly& __other) const noexcept
    { return std::lexicographical_compare_three_way(begin(), end(), __other.begin(), __other.end()); }
};

#endif
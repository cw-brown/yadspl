#ifndef POLY_H
#define POLY_H

#include <memory>
#include <ranges>
#include <iostream>
#include <complex>

template<class __T>
class __iterator{
public:
    typedef std::contiguous_iterator_tag iterator_category;
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
    constexpr pointer operator->() const{return __ptr;}
    constexpr reference operator[](const difference_type& _n) const{return __ptr[_n];}
    constexpr __iterator& operator++(){++__ptr; return *this;}
    constexpr __iterator operator++(int){auto __tmp = *this; ++__ptr; return __tmp;}
    constexpr __iterator& operator--(){--__ptr; return *this;}
    constexpr __iterator operator--(int){auto __tmp = *this; --__ptr; return __tmp;}
    constexpr __iterator& operator+=(const difference_type& n){__ptr += n; return *this;}
    constexpr __iterator& operator-=(const difference_type& n){__ptr -= n; return *this;}
    constexpr friend __iterator operator+(const __iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr += n; return __tmp;}
    constexpr friend __iterator operator+(const difference_type& n, const __iterator& __this){auto __tmp = __this; __tmp.__ptr += n; return __tmp;}
    constexpr friend __iterator operator-(const __iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr -= n; return __tmp;}
    constexpr friend difference_type operator-(const __iterator& __this, const __iterator& other){return __this.__ptr - other.__ptr;}
    constexpr bool operator==(const __iterator& other) const noexcept{return __ptr == other.__ptr;}
    constexpr bool operator!=(const __iterator& other) const noexcept{return __ptr != other.__ptr;}
    constexpr bool operator< (const __iterator& other) const noexcept{return __ptr <  other.__ptr;}
    constexpr bool operator> (const __iterator& other) const noexcept{return __ptr >  other.__ptr;}
    constexpr bool operator<=(const __iterator& other) const noexcept{return __ptr <= other.__ptr;}
    constexpr bool operator>=(const __iterator& other) const noexcept{return __ptr >= other.__ptr;}
};

template<class __T>
class __const_iterator{
public:
    typedef std::contiguous_iterator_tag iterator_category;
    typedef __T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef const __T* pointer;
    typedef const __T& reference;
private:
    pointer __ptr;
public:
    constexpr __const_iterator(): __ptr(nullptr){}
    constexpr __const_iterator(pointer _ptr): __ptr(_ptr){}
    constexpr __const_iterator(const __const_iterator& other): __ptr(other.__ptr){}
    constexpr __const_iterator& operator=(const __const_iterator& other){if(this != &other) __ptr = other.__ptr; return *this;}
    constexpr reference operator*() const{return *__ptr;}
    constexpr pointer operator->() const{return __ptr;}
    constexpr reference operator[](const difference_type& _n) const{return __ptr[_n];}
    constexpr __const_iterator& operator++(){++__ptr; return *this;}
    constexpr __const_iterator operator++(int){auto __tmp = *this; ++__ptr; return __tmp;}
    constexpr __const_iterator& operator--(){--__ptr; return *this;}
    constexpr __const_iterator operator--(int){auto __tmp = *this; --__ptr; return __tmp;}
    constexpr __const_iterator& operator+=(const difference_type& n){__ptr += n; return *this;}
    constexpr __const_iterator& operator-=(const difference_type& n){__ptr -= n; return *this;}
    constexpr friend __const_iterator operator+(const __const_iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr += n; return __tmp;}
    constexpr friend __const_iterator operator+(const difference_type& n, const __const_iterator& __this){auto __tmp = __this; __tmp.__ptr += n; return __tmp;}
    constexpr friend __const_iterator operator-(const __const_iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr -= n; return __tmp;}
    constexpr friend difference_type operator-(const __const_iterator& __this, const __const_iterator& other){return __this.__ptr - other.__ptr;}
    constexpr bool operator==(const __const_iterator& other) const noexcept{return __ptr == other.__ptr;}
    constexpr bool operator!=(const __const_iterator& other) const noexcept{return __ptr != other.__ptr;}
    constexpr bool operator< (const __const_iterator& other) const noexcept{return __ptr <  other.__ptr;}
    constexpr bool operator> (const __const_iterator& other) const noexcept{return __ptr >  other.__ptr;}
    constexpr bool operator<=(const __const_iterator& other) const noexcept{return __ptr <= other.__ptr;}
    constexpr bool operator>=(const __const_iterator& other) const noexcept{return __ptr >= other.__ptr;}
};

template<class __T>
class __reverse_iterator{
public:
    typedef std::contiguous_iterator_tag iterator_category;
    typedef __T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef __T* pointer;
    typedef __T& reference;
private:
    pointer __ptr;
public:
    constexpr __reverse_iterator(): __ptr(nullptr){}
    constexpr __reverse_iterator(pointer _ptr): __ptr(_ptr){}
    constexpr __reverse_iterator(const __reverse_iterator& other): __ptr(other.__ptr){}
    constexpr __reverse_iterator& operator=(const __reverse_iterator& other){if(this != &other) __ptr = other.__ptr; return *this;}
    constexpr reference operator*() const{return *__ptr;}
    constexpr pointer operator->() const{return __ptr;}
    constexpr reference operator[](const difference_type& _n) const{return __ptr[_n];}
    constexpr __reverse_iterator& operator++(){--__ptr; return *this;}
    constexpr __reverse_iterator operator++(int){auto __tmp = *this; --__ptr; return __tmp;}
    constexpr __reverse_iterator& operator--(){++__ptr; return *this;}
    constexpr __reverse_iterator operator--(int){auto __tmp = *this; ++__ptr; return __tmp;}
    constexpr __reverse_iterator& operator+=(const difference_type& n){__ptr -= n; return *this;}
    constexpr __reverse_iterator& operator-=(const difference_type& n){__ptr += n; return *this;}
    constexpr friend __reverse_iterator operator+(const __reverse_iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr -= n; return __tmp;}
    constexpr friend __reverse_iterator operator+(const difference_type& n, const __reverse_iterator& __this){auto __tmp = __this; __tmp.__ptr -= n; return __tmp;}
    constexpr friend __reverse_iterator operator-(const __reverse_iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr += n; return __tmp;}
    constexpr friend difference_type operator-(const __reverse_iterator& __this, const __reverse_iterator& other){return __this.__ptr - other.__ptr;}
    constexpr bool operator==(const __reverse_iterator& other) const noexcept{return __ptr == other.__ptr;}
    constexpr bool operator!=(const __reverse_iterator& other) const noexcept{return __ptr != other.__ptr;}
    constexpr bool operator< (const __reverse_iterator& other) const noexcept{return __ptr <  other.__ptr;}
    constexpr bool operator> (const __reverse_iterator& other) const noexcept{return __ptr >  other.__ptr;}
    constexpr bool operator<=(const __reverse_iterator& other) const noexcept{return __ptr <= other.__ptr;}
    constexpr bool operator>=(const __reverse_iterator& other) const noexcept{return __ptr >= other.__ptr;}
};

template<class __T>
class __const_reverse_iterator{
public:
    typedef std::contiguous_iterator_tag iterator_category;
    typedef __T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef const __T* pointer;
    typedef const __T& reference;
private:
    pointer __ptr;
public:
    constexpr __const_reverse_iterator(): __ptr(nullptr){}
    constexpr __const_reverse_iterator(pointer _ptr): __ptr(_ptr){}
    constexpr __const_reverse_iterator(const __const_reverse_iterator& other): __ptr(other.__ptr){}
    constexpr __const_reverse_iterator& operator=(const __const_reverse_iterator& other){if(this != &other) __ptr = other.__ptr; return *this;}
    constexpr reference operator*() const{return *__ptr;}
    constexpr pointer operator->() const{return __ptr;}
    constexpr reference operator[](const difference_type& _n) const{return __ptr[_n];}
    constexpr __const_reverse_iterator& operator++(){--__ptr; return *this;}
    constexpr __const_reverse_iterator operator++(int){auto __tmp = *this; --__ptr; return __tmp;}
    constexpr __const_reverse_iterator& operator--(){++__ptr; return *this;}
    constexpr __const_reverse_iterator operator--(int){auto __tmp = *this; ++__ptr; return __tmp;}
    constexpr __const_reverse_iterator& operator+=(const difference_type& n){__ptr -= n; return *this;}
    constexpr __const_reverse_iterator& operator-=(const difference_type& n){__ptr += n; return *this;}
    constexpr friend __const_reverse_iterator operator+(const __const_reverse_iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr -= n; return __tmp;}
    constexpr friend __const_reverse_iterator operator+(const difference_type& n, const __const_reverse_iterator& __this){auto __tmp = __this; __tmp.__ptr -= n; return __tmp;}
    constexpr friend __const_reverse_iterator operator-(const __const_reverse_iterator& __this, const difference_type& n){auto __tmp = __this; __tmp.__ptr += n; return __tmp;}
    constexpr friend difference_type operator-(const __const_reverse_iterator& __this, const __const_reverse_iterator& other){return __this.__ptr - other.__ptr;}
    constexpr bool operator==(const __const_reverse_iterator& other) const noexcept{return __ptr == other.__ptr;}
    constexpr bool operator!=(const __const_reverse_iterator& other) const noexcept{return __ptr != other.__ptr;}
    constexpr bool operator< (const __const_reverse_iterator& other) const noexcept{return __ptr <  other.__ptr;}
    constexpr bool operator> (const __const_reverse_iterator& other) const noexcept{return __ptr >  other.__ptr;}
    constexpr bool operator<=(const __const_reverse_iterator& other) const noexcept{return __ptr <= other.__ptr;}
    constexpr bool operator>=(const __const_reverse_iterator& other) const noexcept{return __ptr >= other.__ptr;}
};

/**
 * @brief A polynomial class with math operators. Features specializations for constants and conversions between complex types.
 * It is always constructible via aggregates. Represents a polynomial as p[0] + p[1]*x + ... p[N]*x^N.
 * @tparam __T Type of element. Must be a complete type.
 * @tparam __N The degree of the polynomial.
 */
template<class __T, std::size_t __N>
requires(std::is_move_constructible_v<__T> && std::is_move_assignable_v<__T> && std::is_arithmetic_v<__T>)
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
    typedef __const_iterator<__T> const_iterator;
    typedef __reverse_iterator<__T> reverse_iterator;
    typedef __const_reverse_iterator<__T> const_reverse_iterator;

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
    { if(__n >= __N + 1) 
        std::__throw_out_of_range_fmt("poly::at: __n (which is %zu) >= __N + 1 (which is %zu)", __n, __N + 1);
    return __c[__n]; }

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

    template<size_type __M>
    requires(__N >= __M)
    constexpr poly&
    operator+=(const poly<__T, __M>& __other) noexcept
    {
        for(size_type k = 0; k < __M + 1; ++k)
            __c[k] += __other[k];
        return *this;
    }

    template<size_type __M>
    requires(__N >= __M)
    constexpr poly&
    operator-=(const poly<__T, __M>& __other) noexcept
    {
        for(size_type k = 0; k < __M + 1; ++k)
            __c[k] -= __other[k];
        return *this;
    }

    template<size_type __M>
    requires(__N >= __M)
    [[nodiscard]]
    constexpr poly
    operator+(const poly<__T, __M>& __other) noexcept
    {
        poly<__T, __N> _out = *this;
        for(size_type k = 0; k < __M + 1; ++k)
            _out[k] += __other[k];
        return _out;
    }

    template<size_type __M>
    requires(__M > __N)
    [[nodiscard]]
    constexpr poly<__T, __M>
    operator+(const poly<__T, __M>& __other) noexcept
    {
        poly<__T, __M> _out = __other;
        for(size_type k = 0; k < __N + 1; ++k)
            _out[k] += __c[k];
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator+(const poly& __obj, const __T& __constant) noexcept
    {
        poly _out = __obj;
        _out[0] += __constant;
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator+(const __T& __constant, const poly& __obj) noexcept
    {
        poly _out = __obj;
        _out[0] += __constant;
        return _out;
    }

    template<size_type __M>
    requires(__N >= __M)
    [[nodiscard]]
    constexpr poly
    operator-(const poly<__T, __M>& __other) noexcept
    {
        poly<__T, __N> _out = *this;
        for(size_type k = 0; k < __M + 1; ++k)
            _out[k] -= __other[k];
        return _out;
    }

    template<size_type __M>
    requires(__M > __N)
    [[nodiscard]]
    constexpr poly<__T, __M>
    operator-(const poly<__T, __M>& __other) noexcept
    {
        poly<__T, __M> _out = __T() - __other;
        for(size_type k = 0; k < __N + 1; ++k)
            _out[k] += __c[k];
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator-(const poly& __obj, const __T& __constant) noexcept
    {
        poly _out = __obj;
        _out[0] -= __constant;
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator-(const __T& __constant, const poly& __obj) noexcept
    {
        poly _out = __obj;
        for(size_type k = 0; k < __N + 1; ++k)
            if(k == 0) _out[k] = __constant - __obj[k]; 
            else _out[k] = __T() - __obj[k];
        return _out;
    }

    template<size_type __M>
    requires(__N >= __M)
    [[nodiscard]]
    constexpr poly<__T, __N + __M>
    operator*(const poly<__T, __M>& __other) noexcept
    { // Developed from reference https://www.mathworks.com/help/matlab/ref/conv.html
        poly<__T, __N + __M> _out;
        _out.fill(__T());
        for(size_type k = 0; k <= __N + __M; ++k){
            for(size_type j = std::max(1ULL, k - __M + 1) - 1; j <= std::min(k, __N); ++j)
                _out[k] += __c[j] * __other[k - j];
        }
        return _out;
    }

    template<size_type __M>
    requires(__M > __N)
    [[nodiscard]]
    constexpr poly<__T, __N + __M>
    operator*(const poly<__T, __M>& __other) noexcept
    { // Developed from reference https://www.mathworks.com/help/matlab/ref/conv.html
        poly<__T, __N + __M> _out;
        _out.fill(__T());
        for(size_type k = 0; k <= __N + __M; ++k){
            for(size_type j = std::max(1ULL, k - __N + 1) - 1; j <= std::min(k, __M); ++j)
                _out[k] += __other[j] * __c[k - j];
        }
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator*(const poly& __obj, const __T& __constant) noexcept
    {
        poly _out = __obj;
        for(size_type k = 0; k < __N + 1; ++k)
            _out[k] *= __constant;
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator*(const __T& __constant,const poly& __obj) noexcept
    {
        poly _out = __obj;
        for(size_type k = 0; k < __N + 1; ++k)
            _out[k] *= __constant;
        return _out;
    }

    [[nodiscard]]
    constexpr friend poly
    operator/(const poly& __obj, const __T& __constant) noexcept
    {
        poly _out = __obj;
        for(size_type k = 0; k < __N + 1; ++k)
            _out[k] /= __constant;
        return _out;
    }

    [[nodiscard]]
    constexpr poly
    operator/(const poly& __other) noexcept
    { // Element-wise division operation, because I'm not implementing actual polynomial division.
        poly _out = *this;
        for(size_type k = 0; k < __N + 1; ++k)
            _out[k] /= __other[k];
        return _out;
    }
};
#endif
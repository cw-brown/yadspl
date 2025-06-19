/**
 * @file polynomial.hpp
 * @author cw-brown (https://github.com/cw-brown)
 * @brief Generic polynomial container with math functions.
 * @version 0.1
 * @date 2025-06-18
 */
#ifndef POLY_H
#define POLY_H

#include <memory>
#include <ranges>
#include <iostream>
#include <complex>

#include "basic_iterators.hpp"

/**
 * @brief Polynomial implements a container for real polynomial coefficients. It has built in math functions for addition, subtraction, and multiplication.
 * @tparam _Tp An arithmetic type, either floating or integral.
 */
template<class _Tp>
requires std::is_move_constructible_v<_Tp> && std::is_move_assignable_v<_Tp> && std::is_arithmetic_v<_Tp>
class polynomial{
public:
    typedef _Tp value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef __iterator<value_type> iterator;
    typedef __const_iterator<value_type> const_iterator;
    typedef __reverse_iterator<value_type> reverse_iterator;
    typedef __const_reverse_iterator<value_type> const_reverse_iterator;
private:
    size_type _n; // polynomial order
    pointer _c; // internal array
public:
    constexpr polynomial(): _n(0), _c(new value_type[1]){_c[0] = 0.0;}
    constexpr polynomial(size_type order): _n(order), _c(new value_type[order + 1]){std::uninitialized_fill(begin(), end(), value_type{});}
    constexpr polynomial(std::initializer_list<value_type> init): _n(init.size() - 1), _c(new value_type[init.size()]){std::uninitialized_copy(init.begin(), init.end(), begin());}
    constexpr polynomial(const polynomial& other): _n(other._n), _c(new value_type[other._n + 1]){std::uninitialized_copy(other.begin(), other.end(), begin());}
    constexpr polynomial(polynomial&& other) noexcept: _n(other._n), _c(std::move(other._c)){other._c = nullptr;}

    constexpr ~polynomial(){std::destroy(begin(), end());}

    constexpr polynomial& operator=(const polynomial& other){
        if(this != &other){
            _n = other._n;
            _c = new value_type[other._n + 1];
            std::uninitialized_copy(other.begin(), other.end(), begin());
        }
        return *this;
    }
    constexpr polynomial& operator=(polynomial&& other){
        if(this != &other){
            _n = other._n;
            _c = other._c;
            other._c = nullptr;
        }
        return *this;
    }
    constexpr polynomial& operator=(std::initializer_list<_Tp> init){
        _n = init.size() - 1;
        _c = new value_type[init.size()];
        std::uninitialized_copy(init.begin(), init.end(), begin());
    }
    
    constexpr void fill(const value_type& constant){std::fill(begin(), end(), constant);}
    constexpr void swap(polynomial& other) noexcept(std::is_nothrow_swappable_v<value_type>){std::swap_ranges(begin(), end(), other.begin());}

    constexpr reference operator[](size_type n) noexcept{return _c[n];}
    constexpr const_reference operator[](size_type n) const noexcept{return _c[n];}
    constexpr reference at(size_type n) noexcept{return _c[n];}
    constexpr const_reference at(size_type n) const noexcept{return _c[n];}
    constexpr reference front() noexcept{return *begin();}
    constexpr const_reference front() const noexcept{return *begin();}
    constexpr reference back() noexcept{return *rbegin();}
    constexpr const_reference back() const noexcept{return *rbegin();}
    constexpr pointer data() noexcept{return _c;}
    constexpr const_pointer data() const noexcept{return _c;}

    constexpr size_type order() const noexcept{return _n;}

    constexpr bool operator==(const polynomial& other) const noexcept{return std::equal(begin(), end(), other.begin());}
    constexpr std::compare_three_way_result_t<value_type> operator<=>(const polynomial& other) const noexcept{return std::lexicographical_compare_three_way(begin(), end(), other.begin(), other.end());}

    constexpr polynomial operator+(const polynomial& other) const noexcept{
        polynomial out = _n >= other._n ? *this : other;
        for(size_type k = 0; k < std::min(_n, other._n) + 1; ++k){
            out[k] += _n >= other._n ? other[k] : _c[k];
        }
        return out;
    }
    constexpr polynomial operator+(const value_type& constant) const noexcept{
        polynomial out = *this;
        out[0] += constant;
        return out;
    }
    constexpr friend polynomial operator+(const value_type& constant, const polynomial& obj) noexcept{
        polynomial out = obj;
        out[0] += constant;
        return out;
    }
    
    constexpr polynomial operator-(const polynomial& other) const noexcept{
        polynomial out = _n >= other._n ? *this : other * -1.0;
        for(size_type k = 0; k < std::min(_n, other._n) + 1; ++k){
            out[k] -= _n >= other._n ? other[k] : -1.0 * _c[k];
        }
        return out;
    }
    constexpr polynomial operator-(const value_type& constant) const noexcept{
        polynomial out = *this;
        out[0] -= constant;
        return out;
    }
    constexpr friend polynomial operator-(const value_type& constant, const polynomial& obj) noexcept{
        polynomial out = -1.0 * obj;
        out[0] += constant;
        return out;
    }

    constexpr polynomial operator*(const polynomial& other) const noexcept{
        // Developed from reference https://www.mathworks.com/help/matlab/ref/conv.html
        polynomial out(_n + other._n);
        for(size_type k = 0; k <= _n + other._n; ++k){
            for(size_type j = static_cast<size_type>(std::max(1, (int)k - (int)_n + 1)) - 1; j <= std::min(k, other._n); ++j){
                out._c[k] += other._c[j] * _c[k - j];
            }
        }
        return out;
    }
    constexpr polynomial operator*(const value_type& constant) const noexcept{
        polynomial out = *this;
        for(size_type k = 0; k < _n + 1; ++k) out[k] *= constant;
        return out;
    }
    constexpr friend polynomial operator*(const value_type& constant, const polynomial& obj) noexcept{
        polynomial out = obj;
        for(size_type k = 0; k < obj._n + 1; ++k) out[k] *= constant;
        return out;
    }

    constexpr iterator begin() noexcept{return iterator(&_c[0]);}
    constexpr const_iterator begin() const noexcept{return const_iterator(&_c[0]);}
    constexpr const_iterator cbegin() const noexcept{return const_iterator(&_c[0]);}

    constexpr iterator end() noexcept{return iterator(&_c[_n + 1]);}
    constexpr const_iterator end() const noexcept{return const_iterator(&_c[_n + 1]);}
    constexpr const_iterator cend() const noexcept{return const_iterator(&_c[_n + 1]);}

    constexpr reverse_iterator rbegin() noexcept{return reverse_iterator(&_c[_n]);}
    constexpr const_reverse_iterator rbegin() const noexcept{return const_reverse_iterator(&_c[_n]);}
    constexpr const_reverse_iterator crbegin() const noexcept{return const_reverse_iterator(&_c[_n]);}

    constexpr reverse_iterator rend() noexcept{return reverse_iterator(&_c[-1]);}
    constexpr const_reverse_iterator rend() const noexcept{return const_reverse_iterator(&_c[-1]);}
    constexpr const_reverse_iterator crend() const noexcept{return const_reverse_iterator(&_c[-1]);}
};

/**
 * @brief Complex_Polynomial implements a container for complex polynomial coefficients. It has built in math functions for addition, subtraction, and multiplication.
 * @tparam _Tp An arithmetic type, either floating or integral.
 */
template<class _Tp>
requires std::is_move_constructible_v<_Tp> && std::is_move_assignable_v<_Tp> && std::is_arithmetic_v<_Tp>
class complex_polynomial{
public:
    typedef std::complex<_Tp> value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef __iterator<value_type> iterator;
    typedef __const_iterator<value_type> const_iterator;
    typedef __reverse_iterator<value_type> reverse_iterator;
    typedef __const_reverse_iterator<value_type> const_reverse_iterator;
private:
    size_type _n; // polynomial order
    pointer _c; // internal array
public:
    constexpr complex_polynomial(): _n(0), _c(new value_type[1]){_c[0] = value_type{};}
    constexpr complex_polynomial(size_type order): _n(order), _c(new value_type[order + 1]){std::uninitialized_fill(begin(), end(), value_type{});}
    constexpr complex_polynomial(std::initializer_list<value_type> init): _n(init.size() - 1), _c(new value_type[init.size()]){std::uninitialized_copy(init.begin(), init.end(), begin());}
    constexpr complex_polynomial(const complex_polynomial& other): _n(other._n), _c(new value_type[other._n + 1]){std::uninitialized_copy(other.begin(), other.end(), begin());}
    constexpr complex_polynomial(complex_polynomial&& other): _n(other._n), _c(std::move(other._c)){other._c = nullptr;}
    
    constexpr ~complex_polynomial(){std::destroy(begin(), end());}

    constexpr complex_polynomial& operator=(const complex_polynomial& other){
        if(this != &other){
            _n = other._n;
            _c = new value_type(_n + 1);
            std::uninitialized_copy(other.begin(), other.end(), begin());
        }   
        return *this;
    }
    constexpr complex_polynomial& operator=(complex_polynomial&& other){
        if(this != &other){
            _n = other._n;
            _c = std::move(other._c);
            other._c = nullptr;
        }   
        return *this;
    }
    constexpr complex_polynomial& operator=(std::initializer_list<value_type> init){
        _n = init.size() - 1;
        _c = new value_type[init.size()];
        std::uninitialized_copy(init.begin(), init.end(), begin());
        return *this;
    }
    
    constexpr void fill(const value_type& constant){std::fill(begin(), end(), constant);}
    constexpr void swap(complex_polynomial& other) noexcept(std::is_nothrow_swappable_v<value_type>){std::swap_ranges(begin(), end(), other.begin());}
   
    constexpr reference operator[](size_type n) noexcept{return _c[n];}
    constexpr const_reference operator[](size_type n) const noexcept{return _c[n];}
    constexpr reference at(size_type n) noexcept{return _c[n];}
    constexpr const_reference at(size_type n) const noexcept{return _c[n];}
    constexpr reference front() noexcept{return *begin();}
    constexpr const_reference front() const noexcept{return *begin();}
    constexpr reference back() noexcept{return *rbegin();}
    constexpr const_reference back() const noexcept{return *rbegin();}
    constexpr pointer data() noexcept{return _c;}
    constexpr const_pointer data() const noexcept{return _c;}

    constexpr size_type order() const noexcept{return _n;}

    constexpr iterator begin() noexcept{return iterator(&_c[0]);}
    constexpr const_iterator begin() const noexcept{return const_iterator(&_c[0]);}
    constexpr const_iterator cbegin() const noexcept{return const_iterator(&_c[0]);}

    constexpr iterator end() noexcept{return iterator(&_c[_n + 1]);}
    constexpr const_iterator end() const noexcept{return const_iterator(&_c[_n + 1]);}
    constexpr const_iterator cend() const noexcept{return const_iterator(&_c[_n + 1]);}

    constexpr reverse_iterator rbegin() noexcept{return reverse_iterator(&_c[_n]);}
    constexpr const_reverse_iterator rbegin() const noexcept{return const_reverse_iterator(&_c[_n]);}
    constexpr const_reverse_iterator crbegin() const noexcept{return const_reverse_iterator(&_c[_n]);}

    constexpr reverse_iterator rend() noexcept{return reverse_iterator(&_c[-1]);}
    constexpr const_reverse_iterator rend() const noexcept{return const_reverse_iterator(&_c[-1]);}
    constexpr const_reverse_iterator crend() const noexcept{return const_reverse_iterator(&_c[-1]);}
};

#endif
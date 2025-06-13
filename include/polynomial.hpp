#ifndef POLY_H
#define POLY_H

#include <memory>
#include <ranges>
#include <iostream>
#include <complex>

/**
 * @brief A generic polynomial class. Stores the coefficients as p[0] + p[1]*x + p[2]*x^2 + ... p[n-1]*x^(n-1). Features specializations for conversions between complex types.
 * @tparam T 
 * @tparam Allocator 
 */
template<class T, class Allocator = std::allocator<T>>
class polynomial{
private:
    class it{
    public:
        typedef std::contiguous_iterator_tag iterator_category;
        typedef T value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef T* pointer;
        typedef T& reference;
    private:
        pointer _ptr;
    public:
        constexpr it(): _ptr(nullptr){}
        constexpr it(pointer ptr): _ptr(ptr){}
        constexpr it(const it& other): _ptr(other._ptr){}
        constexpr it(it&& other) noexcept(true): _ptr(std::move(other._ptr)){}
        constexpr it& operator=(const it& other){if(this != &other) _ptr = other._ptr; return this;}
        constexpr it& operator=(it&& other) noexcept{if(this != &other) _ptr = std::move(other._ptr); return this;}
        constexpr reference operator*() const{return *_ptr;}
        constexpr pointer operator->() noexcept{return _ptr;}
        constexpr pointer operator->() const noexcept{return _ptr;}
        constexpr reference operator[](const difference_type& n) noexcept{return *(_ptr + n);}
        constexpr const reference operator[](const difference_type& n) const noexcept{return *(this->_ptr + n);}
        constexpr bool operator==(const it& other) const noexcept{return _ptr == other._ptr;}
        constexpr bool operator!=(const it& other) const noexcept{return _ptr != other._ptr;}
        constexpr bool operator>=(const it& other) const noexcept{return _ptr >= other._ptr;}
        constexpr bool operator<=(const it& other) const noexcept{return _ptr <= other._ptr;}
        constexpr bool operator> (const it& other) const noexcept{return _ptr >  other._ptr;}
        constexpr bool operator< (const it& other) const noexcept{return _ptr <  other._ptr;}
        constexpr it& operator++() noexcept{++_ptr; return *this;}
        constexpr it operator++(int) noexcept{it temp(*this); ++_ptr; return temp;}
        constexpr it& operator--() noexcept{--_ptr; return *this;}
        constexpr it operator--(int) noexcept{it temp(*this); --_ptr; return temp;}
        constexpr it& operator+=(const difference_type& n) noexcept{_ptr += n; return *this;}
        constexpr it& operator-=(const difference_type& n) noexcept{_ptr -= n; return *this;}
        constexpr friend it operator+(const it& pred, const difference_type& n) noexcept{it temp(pred); temp._ptr += n; return temp;}
        constexpr friend it operator+(const difference_type& n, const it& pred) noexcept{it temp(pred); temp._ptr += n; return temp;}
        constexpr friend it operator-(const it& pred, const difference_type& n) noexcept{it temp(pred); temp._ptr -= n; return temp;}
        constexpr friend difference_type operator-(const it& lhs, const it& rhs) noexcept{return lhs._ptr - rhs._ptr;}
    };

    class rit{
    public:
        typedef std::contiguous_iterator_tag iterator_category;
        typedef T value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef T* pointer;
        typedef T& reference;
    private:
        pointer _ptr;
    public:
        constexpr rit(): _ptr(nullptr){}
        constexpr rit(pointer ptr): _ptr(ptr){}
        constexpr rit(const rit& other): _ptr(other._ptr){}
        constexpr rit(rit&& other) noexcept(true): _ptr(std::move(other._ptr)){}
        constexpr rit& operator=(const rit& other){if(this != &other) _ptr = other._ptr; return this;}
        constexpr rit& operator=(rit&& other) noexcept{if(this != &other) _ptr = std::move(other._ptr); return this;}
        constexpr reference operator*() const{return *_ptr;}
        constexpr pointer operator->() noexcept{return _ptr;}
        constexpr pointer operator->() const noexcept{return _ptr;}
        constexpr reference operator[](const difference_type& n) noexcept{return *(_ptr + n);}
        constexpr const reference operator[](const difference_type& n) const noexcept{return *(this->_ptr + n);}
        constexpr bool operator==(const rit& other) const noexcept{return _ptr == other._ptr;}
        constexpr bool operator!=(const rit& other) const noexcept{return _ptr != other._ptr;}
        constexpr bool operator>=(const rit& other) const noexcept{return _ptr >= other._ptr;}
        constexpr bool operator<=(const rit& other) const noexcept{return _ptr <= other._ptr;}
        constexpr bool operator> (const rit& other) const noexcept{return _ptr >  other._ptr;}
        constexpr bool operator< (const rit& other) const noexcept{return _ptr <  other._ptr;}
        constexpr rit& operator--() noexcept{++_ptr; return *this;}
        constexpr rit operator--(int) noexcept{auto temp(*this); ++_ptr; return temp;}
        constexpr rit& operator++() noexcept{--_ptr; return *this;}
        constexpr rit operator++(int) noexcept{auto temp(*this); --_ptr; return temp;}
        constexpr rit& operator-=(const difference_type& n) noexcept{_ptr += n; return *this;}
        constexpr rit& operator+=(const difference_type& n) noexcept{_ptr -= n; return *this;}
        constexpr friend rit operator+(const rit& pred, const difference_type& n) noexcept{auto temp(pred); temp._ptr -= n; return temp;}
        constexpr friend rit operator+(const difference_type& n, const rit& pred) noexcept{auto temp(pred); temp._ptr -= n; return temp;}
        constexpr friend rit operator-(const rit& pred, const difference_type& n) noexcept{auto temp(pred); temp._ptr += n; return temp;}
        constexpr friend difference_type operator-(const rit& lhs, const rit& rhs) noexcept{return rhs._ptr - lhs._ptr;}
    };
public:
    typedef T value_type;
    typedef Allocator allocator_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef it iterator;
    typedef it const_iterator;
    typedef rit reverse_iterator;
    typedef rit const_reverse_iterator;
private:
    Allocator _alloc;
    pointer _coef;
    size_type _n;
public:
    /**
     * @brief Default constructs a polynomial of degree 0 as a constant coefficient of 0.
     */
    constexpr polynomial()
        : polynomial(Allocator()){}
    
    /**
     * @brief Default constructs a polynomial of degree 0 as a constant coefficient of 0.
     * @param alloc The allocator.
     */
    constexpr explicit polynomial(const Allocator& alloc) noexcept(true)
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, 1)), _n(0){std::uninitialized_value_construct(begin(), end());}
    
    /**
     * @brief Construct a polynomial of degree n.
     * @param n The polynomial degree.
     * @param alloc The allocator.
     */
    constexpr explicit polynomial(const size_type& n, const Allocator& alloc = Allocator())
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, n + 1)), _n(n){std::uninitialized_value_construct(begin(), end());}
    
    /**
     * @brief Construct a polynomial of degree n with all coefficients initialized to a value.
     * @param n The polynomial degree.
     * @param value The initialized value.
     * @param alloc The allocator.
     */
    constexpr explicit polynomial(const size_type& n, const_reference value, const Allocator& alloc = Allocator())
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, n + 1)), _n(n){std::uninitialized_fill(begin(), end(), value);}
    
    /**
     * @brief Construct a polynomial from the range [first, last).
     * @param first An input iterator.
     * @param last An input iterator.
     * @param alloc The allocator.
     */
    template<class InputIt> 
    constexpr polynomial(InputIt first, InputIt last, const Allocator& alloc = Allocator()) requires(std::input_iterator<InputIt>)
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, std::distance(first, last))), _n(std::distance(first, last) - 1){
            for(size_type i = 0; first != last; ++first, ++i){
                std::allocator_traits<Allocator>::construct(_alloc, _coef + i, *first);
            }
        }
    
    /**
     * @brief Construct a polynomial with coefficients from the range rg.
     * @param rg The range.
     * @param alloc The allocator.
     */
    template<class R> 
    constexpr polynomial(std::from_range_t, R&& rg, const Allocator& alloc = Allocator()) requires(std::ranges::input_range<R> && std::convertible_to<std::ranges::range_reference_t<R>, T>)
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, std::ranges::distance(rg))), _n(std::ranges::distance(rg) - 1){std::uninitialized_copy(rg.begin(), rg.end(), _coef);}
    
    /**
     * @brief Constructs a polynomial with the contents of other.
     * @param other A polynomial of the same type.
     */
    constexpr polynomial(const polynomial& other)
        : _alloc(std::allocator_traits<Allocator>::select_on_container_copy_construction(other.get_allocator())), _coef(std::allocator_traits<Allocator>::allocate(_alloc, other._n + 1)), _n(other._n){std::uninitialized_copy(other.begin(), other.end(), _coef);}
   
    /**
     * @brief Moves a polynomial with the contents of other.
     * @param other A polynomial of the same type.
     */
    constexpr polynomial(polynomial&& other) noexcept
        : _alloc(std::move(other.get_allocator())), _coef(std::move(other._coef)), _n(other._n){other._coef = nullptr;}
    
    /**
     * @brief Constructs a polynomial with the contents of other and the specified allocator.
     * @param other A polynomial of the same type.
     * @param alloc The allocator.
     */
    constexpr polynomial(const polynomial& other, const Allocator& alloc)
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, other._n + 1)), _n(other._n){std::uninitialized_copy(other.begin(), other.end(), _coef);}
    
    /**
     * @brief Move construct a polynomial with the contents of other and the specified allocator.
     * @param other A polynomial of the same type.
     * @param alloc The allocator.
     */
    constexpr polynomial(polynomial&& other, const Allocator& alloc)
        : _alloc(alloc), _n(other._n){
        if(alloc != other.get_allocator()){
            _coef = std::allocator_traits<Allocator>::allocate(_alloc, other._n + 1);
            std::uninitialized_move(other.begin(), other.end(), _coef);
        } else{
            _coef = std::move(other._coef);
            other._coef = nullptr;
        }
    }
    
    /**
     * @brief Construct a polynomial with coeffients copied from an initializer list.
     * @param init The initializer list.
     * @param alloc The allocator.
     */
    constexpr polynomial(std::initializer_list<T> init, const Allocator& alloc = Allocator())
        : _alloc(alloc), _coef(std::allocator_traits<Allocator>::allocate(_alloc, init.size())), _n(init.size() - 1){std::uninitialized_copy(init.begin(), init.end(), _coef);}

    /**
     * @brief Destroy the polynomial.
     */
    constexpr ~polynomial(){
        std::destroy(begin(), end());
        std::allocator_traits<Allocator>::deallocate(_alloc, _coef, _n);
    }

    /**
     * @brief Replaces the coefficients with a copy from other.
     * @param other A polynomial of the same type.
     * @return polynomial& 
     */
    constexpr polynomial& operator=(const polynomial& other){
        if(this != &other){
            pointer _temp;
            auto _a = other._alloc;
            if(std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) _temp = std::allocator_traits<Allocator>::allocate(_a, other._n + 1);
            else _temp = std::allocator_traits<Allocator>::allocate(_alloc, other._n + 1);
            std::uninitialized_copy(other.begin(), other.end(), _temp);
            std::destroy(begin(), end());
            std::allocator_traits<Allocator>::deallocate(_alloc, _coef, _n + 1);
            if(std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) _alloc = other._alloc;
            _n = other._n;
            _coef = _temp;
        }
        return *this;
    }

    /**
     * @brief Replaces the coefficients with the contents from other using move semantics.
     * @param other A polynomial of the same type.
     * @return polynomial& 
     */
    constexpr polynomial& operator=(polynomial&& other) noexcept{
        if(this != &other){
            auto _a = _alloc;
            if(std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value){
                _alloc = other._alloc;
                _coef = other._coef;
            } else if(_alloc != other._alloc){
                pointer _temp = std::allocator_traits<Allocator>::allocate(_alloc, other._n + 1);
                std::uninitialized_move(other.begin(), other.end(), _temp);
                std::destroy(begin(), end());
                std::allocator_traits<Allocator>::deallocate(_alloc, _coef, _n + 1);
                _coef = _temp;
            } else{
                _coef = other._coef;
            }
            _n = other._n;
            other._coef = nullptr;
        }
        return *this;
    }

    /**
     * @brief Replaces the coefficients with the contents of an initializer list.
     * @param ilist The initializer list.
     * @return polynomial& 
     */
    constexpr polynomial& operator=(const std::initializer_list<T>& ilist){
        pointer _temp = std::allocator_traits<Allocator>::allocate(_alloc, ilist.size());
        std::uninitialized_copy(ilist.begin(), ilist.end(), _temp);
        std::destroy(begin(), end());
        std::allocator_traits<Allocator>::deallocate(_alloc, _coef, _n + 1);
        _n = ilist.size() - 1;
        _coef = _temp;
        return *this;
    }

    /**
     * @brief Get the allocator object
     * @return allocator_type 
     */
    constexpr allocator_type get_allocator() const noexcept{return _alloc;}
    
    /**
     * @brief Get the coefficients of the polynomial.
     * @return pointer 
     */
    constexpr pointer coefficients() const noexcept{return _coef;}

    /**
     * @brief Get the degree of the polynomial.
     * @return size_type 
     */
    constexpr size_type degree() const noexcept{return _n;}

    /**
     * @brief Get the coefficient at a specific position.
     * @param pos The position of the coefficient.
     * @return reference 
     */
    constexpr reference at(const size_type& pos){
        if(pos >= _n + 1) throw std::out_of_range("polynomial: access out of bounds");
        return *(_coef + pos);
    }

    /**
     * @brief Get the read-only coefficient at a specific position.
     * @param pos The position of the coefficient.
     * @return const_reference 
     */
    constexpr const_reference at(const size_type& pos) const{
        if(pos >= _n + 1) throw std::out_of_range("polynomial: access out of bounds");
        return *(_coef + pos);
    }

    /**
     * @brief Get the coefficient at a specific position.
     * @param pos The position of the coefficient.
     * @return reference 
     */
    constexpr reference operator[](const size_type& pos){return _coef[pos];}

    /**
     * @brief Get the read-only coefficient at a specific position.
     * @param pos The position of the coefficient.
     * @return const_reference 
     */
    constexpr const_reference operator[](const size_type& pos) const{return _coef[pos];}

    /**
     * @brief Get a reference to the first coefficient.
     * @return reference 
     */
    constexpr reference front(){return *begin();}

    /**
     * @brief Get a read-only reference to the first coefficient.
     * @return const_reference 
     */
    constexpr const_reference front() const{return *cbegin();}

    /**
     * @brief Get a reference to the last coefficient.
     * @return reference 
     */
    constexpr reference back(){return *std::prev(end());}

    /**
     * @brief Get a read-only reference to the last coefficient.
     * @return const_reference 
     */
    constexpr const_reference back() const{return *std::prev(cend());}

    
    constexpr iterator begin(){return iterator(&_coef[0]);}
    constexpr const_iterator begin() const{return const_iterator(&_coef[0]);}
    constexpr const_iterator cbegin() const{return const_iterator(&_coef[0]);}

    constexpr iterator end(){return iterator(&_coef[_n + 1]);}
    constexpr const_iterator end() const{return const_iterator(&_coef[_n + 1]);}
    constexpr const_iterator cend() const{return const_iterator(&_coef[_n + 1]);}

    constexpr reverse_iterator rbegin(){return reverse_iterator(&_coef[_n]);}
    constexpr const_reverse_iterator rbegin() const{return const_reverse_iterator(&_coef[_n]);}
    constexpr const_reverse_iterator crbegin() const{return const_reverse_iterator(&_coef[_n]);}

    constexpr reverse_iterator rend(){return reverse_iterator(&_coef[-1]);}
    constexpr const_reverse_iterator rend() const{return const_reverse_iterator(&_coef[-1]);}
    constexpr const_reverse_iterator crend() const{return const_reverse_iterator(&_coef[-1]);}


};
#endif
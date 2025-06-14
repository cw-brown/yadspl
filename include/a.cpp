template<class __T>
class __reverse_iterator{
public:
    typedef std::random_access_iterator_tag iterator_category;
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
#ifndef RING_HPP
#define RING_HPP

#include <iterator>
#include <memory>
#include <utility>
#include <initializer_list>
#include <ranges>
#include <type_traits>

namespace std{
template<class T, class Allocator = allocator<T>>
class ring{
private:
    template<class D, class Ref, class Ptr, class ContainerPtr>
    class it{
    public:
        typedef contiguous_iterator_tag iterator_category;
        typedef D value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef Ptr pointer;
        typedef Ref reference;
        typedef ContainerPtr container_pointer;
    private:
        container_pointer _buf;
        pointer _ptr;
        size_type _idx;
        difference_type _cnt;
        pointer _sntl;
    public:
        it(): _buf(nullptr), _ptr(nullptr), _idx(0), _cnt(0), _sntl(nullptr){}
        it(container_pointer container, pointer ptr, size_type index, size_type count): _buf(container), _ptr(ptr), _idx(index), _cnt(count), _sntl(&container->data()[-1]){}
        it(const it& other): _buf(other._buf), _ptr(other._ptr), _idx(other._idx), _cnt(other._cnt), _sntl(other._sntl){}
        it(it&& other): _buf(move(other._buf)), _ptr(move(other._ptr)), _idx(other._idx), _cnt(other._cnt), _sntl(move(other._sntl)){
            other._buf = nullptr;
            other._ptr = nullptr;
            other._sntl = nullptr;
        }
        constexpr it& operator=(const it& other){
            if(this != other){
                this->_buf = other._buf;
                this->_ptr = other._ptr;
                this->_idx = other._idx;
                this->_cnt = other._cnt;
                this->_sntl = other._sntl;
            }
            return *this;
        }
        constexpr it& operator=(it&& other){
            if(this != other){
                this->_buf = move(other._buf);
                this->_ptr = move(other._ptr);
                this->_idx = other._idx;
                this->_cnt = other._cnt;
                this->_sntl = move(other._sntl);
                other._buf = nullptr;
                other._ptr = nullptr;
                other._sntl = nullptr;
            }
            return *this;
        }
        
        constexpr reference operator*() const{
            return *this->_ptr;
        }
        constexpr pointer operator->() noexcept{
            return this->_ptr;
        }
        constexpr pointer operator->() const noexcept{
            return this->_ptr;
        }
        constexpr reference operator[](const difference_type& n) noexcept{
            if(this->_idx == 0 && n < 0){
                return this->_buf->data()[this->_buf->max_size() - 1 - n];
            } else{
                return this->_buf->data()[(this->_idx + n) % this->_buf->max_size()];
            }
        }
        constexpr const reference operator[](const difference_type& n) const noexcept{
            if(this->_idx == 0 && n < 0){
                return this->_buf->data()[this->_buf->max_size() - 1 - n];
            } else{
                return this->_buf->data()[(this->_idx + n) % this->_buf->max_size()];
            }
        }
        
        constexpr bool operator==(const it& other) const noexcept{
            if(this->_cnt == 0 && (other._ptr == this->_sntl)){
                return true;
            } else if((this->_idx == other._idx) && (this->_ptr == other._ptr)){
                return true;
            } else if(this->_ptr == this->_sntl && other._ptr == other._sntl){
                return true;
            } else{
                return false;
            }
        }
        constexpr bool operator!=(const it& other) const noexcept{
            return !(*this == other);
        }
        constexpr bool operator< (const it& other) const noexcept{
            if(other._ptr == other._sntl && this->_ptr != other._sntl){
                return true; // a pointer is always less than a sentinel
            } else if(this->_buf->_tail > this->_buf->_head){
                return true;
            } else if(this->_idx < other._idx){
                return true;
            } else{
                return false;
            }
        }
        constexpr bool operator> (const it& other) const noexcept{
            return !(*this < other);
        }
        constexpr bool operator<=(const it& other) const noexcept{
            return *this == other || *this < other ? true : false;
        }
        constexpr bool operator>=(const it& other) const noexcept{
            return *this == other || *this > other ? true : false;
        }

        constexpr it& operator++() noexcept{
            this->_idx = (this->_idx + 1) % this->_buf->max_size();
            this->_ptr = &this->_buf->data()[this->_idx];
            --this->_cnt;
            return *this;
        }
        constexpr it operator++(int) noexcept{
            auto temp = *this;
            this->_idx = (this->_idx + 1) % this->_buf->max_size();
            this->_ptr = &this->_buf->data()[this->_idx];
            --this->_cnt;
            return temp;
        }
        constexpr it& operator--() noexcept{
            this->_idx = this->_idx == 0 ? this->_buf->max_size() - 1 : this->_idx - 1;
            this->_ptr = &this->_buf->data()[this->_idx];
            ++this->_cnt;
            return *this;
        }
        constexpr it operator--(int) noexcept{
            auto temp = *this;
            this->_idx = this->_idx == 0 ? this->_buf->max_size() - 1 : this->_idx - 1;
            this->_ptr = &this->_buf->data()[this->_idx];
            ++this->_cnt;
            return temp;
        }

        constexpr it& operator+=(const difference_type& n) noexcept{
            this->_idx = (this->_idx + n) % this->_buf->max_size();
            this->_ptr = &this->_buf->data()[this->_idx];
            this->_cnt -= n;
            return *this;
        }
        constexpr it& operator-=(const difference_type& n) noexcept{
            this->_idx = this->_idx == 0 ? this->_buf->max_size() - 1 - n : this->_idx - n;
            this->_ptr = &this->_buf->data()[this->_idx];
            this->_cnt += n;
            return *this;
        }
        constexpr friend it operator+(const it& pred, const difference_type& n) noexcept{
            auto temp = pred;
            temp._idx = (pred._idx + n) % pred._buf->max_size();
            temp._ptr = &pred._buf->data()[temp._idx];
            temp._cnt = pred._cnt - n;
            return temp;
        }
        constexpr friend it operator+(const difference_type& n, const it& pred) noexcept{
            auto temp = pred;
            temp._idx = (pred._idx + n) % pred._buf->max_size();
            temp._ptr = &pred._buf->data()[temp._idx];
            temp._cnt = pred._cnt - n;
            return temp;
        }
        constexpr friend it operator-(const it& pred, const difference_type& n) noexcept{
            auto temp = pred;
            temp._idx = pred._idx == 0 ? pred._buf->max_size() - 1 - n : pred._idx - n;
            if(pred._idx == 0 && pred._ptr == pred._sntl){temp._idx = pred._buf->max_size() - n;} // why?
            temp._ptr = &pred._buf->data()[temp._idx];
            temp._cnt = pred._cnt + n;
            return temp;
        }
        constexpr friend difference_type operator-(const it& lhs, const it& rhs) noexcept{
            return abs(lhs._cnt - rhs._cnt);
        }
    };
    template<class D, class Ref, class Ptr, class ContainerPtr>
    class rit{
    public:
        typedef contiguous_iterator_tag iterator_category;
        typedef D value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef Ptr pointer;
        typedef Ref reference;
        typedef ContainerPtr container_pointer;
    private:
        container_pointer _buf;
        pointer _ptr;
        size_type _idx;
        difference_type _cnt;
        pointer _sntl;
    public:
        rit(): _buf(nullptr), _ptr(nullptr), _idx(0), _cnt(0), _sntl(nullptr){}
        rit(container_pointer container, pointer ptr, size_type index, size_type count): _buf(container), _ptr(ptr), _idx(index), _cnt(count), _sntl(&container->data()[-1]){}
        rit(const rit& other): _buf(other._buf), _ptr(other._ptr), _idx(other._idx), _cnt(other._cnt), _sntl(other._sntl){}
        rit(rit&& other): _buf(move(other._buf)), _ptr(move(other._ptr)), _idx(other._idx), _cnt(other._cnt), _sntl(move(other._sntl)){
            other._buf = nullptr;
            other._ptr = nullptr;
            other._sntl = nullptr;
        }
        constexpr rit& operator=(const rit& other){
            if(this != other){
                this->_buf = other._buf;
                this->_ptr = other._ptr;
                this->_idx = other._idx;
                this->_cnt = other._cnt;
                this->_sntl = other._sntl;
            }
            return *this;
        }
        constexpr rit& operator=(rit&& other){
            if(this != other){
                this->_buf = move(other._buf);
                this->_ptr = move(other._ptr);
                this->_idx = other._idx;
                this->_cnt = other._cnt;
                this->_sntl = move(other._sntl);
                other._buf = nullptr;
                other._ptr = nullptr;
                other._sntl = nullptr;
            }
            return *this;
        }
        
        constexpr reference operator*() const{
            return *this->_ptr;
        }
        constexpr pointer operator->() noexcept{
            return this->_ptr;
        }
        constexpr pointer operator->() const noexcept{
            return this->_ptr;
        }
        constexpr reference operator[](const difference_type& n) noexcept{
            if(this->_idx == 0 && n < 0){
                return this->_buf->data()[this->_buf->max_size() - 1 - n];
            } else{
                return this->_buf->data()[(this->_idx + n) % this->_buf->max_size()];
            }
        }
        constexpr const reference operator[](const difference_type& n) const noexcept{
            if(this->_idx == 0 && n < 0){
                return this->_buf->data()[this->_buf->max_size() - 1 - n];
            } else{
                return this->_buf->data()[(this->_idx + n) % this->_buf->max_size()];
            }
        }
        
        constexpr bool operator==(const rit& other) const noexcept{
            if(this->_cnt == 0 && (other._ptr == this->_sntl)){
                return true;
            } else if((this->_idx == other._idx) && (this->_ptr == other._ptr)){
                return true;
            } else if(this->_ptr == this->_sntl && other._ptr == other._sntl){
                return true;
            } else{
                return false;
            }
        }
        constexpr bool operator!=(const rit& other) const noexcept{
            return !(*this == other);
        }
        constexpr bool operator< (const rit& other) const noexcept{
            if(other._ptr == other._sntl && this->_ptr != other._sntl){
                return true; // a pointer is always less than a sentinel
            } else if(this->_buf->_tail > this->_buf->_head){
                return true;
            } else if(this->_idx < other._idx){
                return true;
            } else{
                return false;
            }
        }
        constexpr bool operator> (const rit& other) const noexcept{
            return !(*this < other);
        }
        constexpr bool operator<=(const rit& other) const noexcept{
            return *this == other || *this < other ? true : false;
        }
        constexpr bool operator>=(const rit& other) const noexcept{
            return *this == other || *this > other ? true : false;
        }

        constexpr rit& operator--() noexcept{
            this->_idx = (this->_idx + 1) % this->_buf->max_size();
            this->_ptr = &this->_buf->data()[this->_idx];
            ++this->_cnt;
            return *this;
        }
        constexpr rit operator--(int) noexcept{
            auto temp = *this;
            this->_idx = (this->_idx + 1) % this->_buf->max_size();
            this->_ptr = &this->_buf->data()[this->_idx];
            ++this->_cnt;
            return temp;
        }
        constexpr rit& operator++() noexcept{
            this->_idx = this->_idx == 0 ? this->_buf->max_size() - 1 : this->_idx - 1;
            this->_ptr = &this->_buf->data()[this->_idx];
            --this->_cnt;
            return *this;
        }
        constexpr rit operator++(int) noexcept{
            auto temp = *this;
            this->_idx = this->_idx == 0 ? this->_buf->max_size() - 1 : this->_idx - 1;
            this->_ptr = &this->_buf->data()[this->_idx];
            --this->_cnt;
            return temp;
        }

        constexpr rit& operator-=(const difference_type& n) noexcept{
            this->_idx = (this->_idx + n) % this->_buf->max_size();
            this->_ptr = &this->_buf->data()[this->_idx];
            this->_cnt += n;
            return *this;
        }
        constexpr rit& operator+=(const difference_type& n) noexcept{
            this->_idx = this->_idx == 0 ? this->_buf->max_size() - 1 - n : this->_idx - n;
            this->_ptr = &this->_buf->data()[this->_idx];
            this->_cnt -= n;
            return *this;
        }
        constexpr friend rit operator+(const rit& pred, const difference_type& n) noexcept{
            auto temp = pred;
            temp._idx = pred._idx == 0 ? pred._buf->max_size() - 1 - n : pred._idx - n;
            if(pred._idx == 0 && pred._ptr == pred._sntl){temp._idx = pred._buf->max_size() - n;} // why?
            temp._ptr = &pred._buf->data()[temp._idx];
            temp._cnt = pred._cnt - n;
            return temp;
        }
        constexpr friend rit operator+(const difference_type& n, const rit& pred) noexcept{
            auto temp = pred;
            temp._idx = pred._idx == 0 ? pred._buf->max_size() - 1 - n : pred._idx - n;
            if(pred._idx == 0 && pred._ptr == pred._sntl){temp._idx = pred._buf->max_size() - n;} // why?
            temp._ptr = &pred._buf->data()[temp._idx];
            temp._cnt = pred._cnt - n;
            return temp;
        }
        constexpr friend rit operator-(const rit& pred, const difference_type& n) noexcept{
            auto temp = pred;
            temp._idx = (pred._idx + n) % pred._buf->max_size();
            temp._ptr = &pred._buf->data()[temp._idx];
            temp._cnt = pred._cnt + n;
            return temp;
        }
        constexpr friend difference_type operator-(const rit& lhs, const rit& rhs) noexcept{
            return abs(lhs._cnt - rhs._cnt);
        }
    };
public:
    typedef T value_type;
    typedef Allocator allocator_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef it<T, T&, T*, ring<T>*> iterator;
    typedef it<T, const T&, const T*, const ring<T>*> const_iterator;
    typedef rit<T, T&, T*, ring<T>*> reverse_iterator;
    typedef rit<T, const T&, const T*, const ring<T>*> const_reverse_iterator;
private:
    size_type _max_size;
    size_type _head;
    size_type _tail;
    size_type _size;
    allocator_type _alloc;
    pointer _buffer;
    constexpr void _M_range_check(size_type __n) const{
        if(__n >= (this->_head + this->_max_size - this->_tail)){
            throw out_of_range("ring::_M_range_check: __n (which is " + to_string(__n) + ") >= this->max_size() (which is " + to_string(this->_max_size) + "), or references uninitialized data beyond this->_head (which is " + to_string(this->_head) + ")");
        }
    }
    constexpr void _M_resize_check(size_type __n) const{
        if(__n > this->_max_size){
            throw out_of_range("ring::_M_resize_check: __n (which is " + to_string(__n) + ") > this->max_size() (which is " + to_string(this->_max_size) + ")");
        }
    }
    constexpr bool _full() const{
        return this->_size == this->_max_size;
    }
    constexpr void _incr(){
        if(this->_full()){
            this->_tail = (this->_tail + 1) % this->_max_size;
        } else{
            ++this->_size;
        }
        this->_head = (this->_head + 1) % this->_max_size;
    };
public:
    /**
     * @brief Default constructor that makes no elements nor allocates memory.
     */
    constexpr ring()
    : ring(Allocator()){}

    /**
     * @brief Constructs a ring with no elements.
     */
    constexpr explicit ring(const Allocator& alloc) noexcept
    : _max_size(0), _head(0), _tail(0), _size(0), _alloc(alloc), _buffer(allocator_traits<Allocator>::allocate(this->_alloc, 0)){}
    
    /**
     * @brief Construct a new ring object with default constructed elements.
     * @param count The number of elements of the ring.
     * @param alloc An allocator.
     */
    constexpr explicit ring(size_type count, const Allocator& alloc = Allocator())
    : _max_size(count), _head(0), _tail(0), _size(0), _alloc(alloc), _buffer(allocator_traits<Allocator>::allocate(this->_alloc, count)){}

    /**
     * @brief Constructs a ring with elements of a specific value.
     * @param count The number of elements of the ring.
     * @param value The value to assign to each element.
     * @param alloc An allocator.
     */
    constexpr explicit ring(size_type count, const_reference value, const Allocator& alloc = Allocator())
    : _max_size(count), _head(0), _tail(0), _size(0), _alloc(alloc), _buffer(allocator_traits<Allocator>::allocate(this->_alloc, count)){
        for(size_type i = 0; i < count; ++i){
            construct_at(this->_buffer + i, value);
            this->_incr();
        }
    }

    /**
     * @brief Builds a ring from a range.
     * @param first An input iterator.
     * @param last An input iterator.
     * @param alloc An allocator.
     */
    template<class InputIt>
    constexpr ring(InputIt first, InputIt last, const Allocator& alloc = Allocator()) requires (input_iterator<InputIt>)
    : _max_size(distance(first, last)), _head(0), _tail(0), _size(0), _alloc(alloc), _buffer(allocator_traits<Allocator>::allocate(this->_alloc, distance(first, last))){
        for(; first != last; ++first){
            construct_at(this->_buffer + this->_head, *first);
            this->_incr();
        }
    }

    /**
     * @brief Builds a ring from a container compatible range object.
     * @param rg std::ranges::input_range container compatible range.
     * @param alloc An allocator.
     */
    template<class R>
    constexpr ring(from_range_t, R&& rg, const Allocator& alloc = Allocator()) requires (ranges::input_range<R> && convertible_to<ranges::range_reference_t<R>, T>)
    : _max_size(ranges::distance(rg)), _head(0), _tail(0), _size(0), _alloc(alloc), _buffer(allocator_traits<Allocator>::allocate(this->_alloc, ranges::distance(rg))){
        for(auto&& v: rg){
            construct_at(this->_buffer + this->_head, v);
            this->_incr();
        }
    }
    
    /**
     * @brief Copy construct a ring object.
     * @param other A ring with identical element and allocator type
     */
    constexpr ring(const ring& other)
    : _max_size(other._max_size), _head(other._head), _tail(other._tail), _size(other._size){
        this->_alloc = allocator_traits<Allocator>::select_on_container_copy_construction(other.get_allocator());
        this->_buffer = allocator_traits<Allocator>::allocate(this->_alloc, this->_max_size);
        for(size_type i = 0; i < other._max_size; ++i){
            construct_at(this->_buffer + i, other._buffer[i]);
        }
    }
    
    /**
     * @brief Move construct a ring object.
     * @param other A ring of identical element and allocator types.
     */
    constexpr ring(ring&& other) noexcept
    : _max_size(other._max_size), _head(other._head), _tail(other._tail), _size(other._size){
        this->_alloc = move(other.get_allocator());
        this->_buffer = other._buffer;
        other._buffer = nullptr;
    }

    /**
     * @brief Copy constructs a ring object with a different allocator.
     * @param other Another ring with identical element and allocator type.
     * @param alloc An allocator.
     */
    constexpr ring(const ring& other, const type_identity_t<Allocator>& alloc)
    : _max_size(other._max_size), _head(other._head), _tail(other._tail), _size(other._size), _alloc(alloc){
        this->_buffer = allocator_traits<Allocator>::allocate(this->_alloc, this->_max_size);
        for(size_type i = 0; i < other._max_size; ++i){
            construct_at(this->_buffer + i, other._buffer[i]);
        }
    }

    /**
     * @brief Move construct a ring object with a specified allocator.
     * @param other A ring of identical element and allocator types.
     * @param alloc An allocator.
     */
    constexpr ring(ring&& other, const type_identity_t<Allocator>& alloc)
    : _max_size(other._max_size), _head(other._head), _tail(other._tail), _size(other._size), _alloc(alloc){
        if(alloc == other.get_allocator()){
            this->_buffer = other._buffer;
            other._buffer = nullptr;
        } else{
            this->_buffer = allocator_traits<Allocator>::allocate(this->_alloc, this->_max_size);
            for(size_type i = 0; i < other._max_size; ++i){
                construct_at(this->_buffer + i, move_if_noexcept(other._buffer[i]));
            }
        }
    }

    /**
     * @brief Build a ring from an initializer list.
     * @param init The initializer list.
     * @param alloc An allocator.
     */
    constexpr ring(initializer_list<T> init, const Allocator& alloc = Allocator())
    : _max_size(init.size()), _head(0), _tail(0), _size(0), _alloc(alloc), _buffer(allocator_traits<Allocator>::allocate(this->_alloc, init.size())){
        for(auto&& v : init){
            this->_buffer[_head] = v;
            this->_incr();
        }
    }

    /**
     * @brief Destroy the ring object and deallocate all memory.
     */
    constexpr ~ring(){
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
    }

    /**
     * @brief Ring copy assignment operation.
     * @param other A ring of identical element and allocator types.
     * @return ring&
     */
    constexpr ring& operator=(const ring& other){
        if(this != &other){
            pointer _temp;
            auto _a = other._alloc;
            if(allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) _temp = allocator_traits<Allocator>::allocate(_a, other._max_size);
            else _temp = allocator_traits<Allocator>::allocate(this->_alloc, other._max_size);
            uninitialized_copy(other._buffer, other._buffer + other._max_size, _temp);
            destroy(this->begin(), this->end());
            allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
            if(allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) this->_alloc = other._alloc;
            this->_max_size = other._max_size;
            this->_head = other._head;
            this->_tail = other._tail;
            this->_size = other._size;
            this->_buffer = _temp;
        }
        return *this;
    }
    
    /**
     * @brief Ring move assignment operation.
     * @param other A ring of identical element and allocator types.
     * @return ring& 
     */
    constexpr ring& operator=(ring&& other) noexcept{
        if(this != &other){
            destroy(this->begin(), this->end());
            allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
            if(allocator_traits<Allocator>::propagate_on_container_move_assignment::value){
                this->_alloc = other._alloc;
                this->_buffer = other._buffer;
            } else if(this->_alloc != other._alloc){
                pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, other._size);
                uninitialized_move(other._buffer, other._buffer + other._max_size, _temp);
                this->_buffer = _temp;
            } else{
                this->_buffer = other._buffer;
            }   
            this->_max_size = other._max_size;
            this->_head = other._head;
            this->_tail = other._tail;
            this->_size = other._size;
            other._buffer = nullptr;
        }
        return *this;
    }
    
    /**
     * @brief Ring aggregate assignment operation.
     * @param ilist The aggregate list.
     * @return ring& 
     */
    constexpr ring& operator=(initializer_list<value_type> ilist){
        pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, ilist.size());
        uninitialized_copy(ilist.begin(), ilist.end(), _temp);
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
        this->_max_size = ilist.size();
        this->_head = 0;
        this->_tail = 0;
        this->_size = ilist.size();
        this->_buffer = _temp;
        return *this;
    }
    
    /**
     * @brief Assigns a value to a ring.
     * @param count Number of elements to assign.
     * @param value Value to be assigned.
     */
    constexpr void assign(size_type count, const_reference value){
        pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, count);
        uninitialized_fill_n(_temp, count, value);
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
        this->_max_size = count;
        this->_head = 0;
        this->_tail = 0;
        this->_size = count;
        this->_buffer = _temp;
    }
    
    /**
     * @brief Assigns a range to a ring.
     * @param first An input iterator.
     * @param last An input iterator.
     */
    template<class InputIt>
    constexpr void assign(InputIt first, InputIt last) requires (input_iterator<InputIt>){
        pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, distance(first, last));
        uninitialized_copy(first, last, _temp);
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
        this->_max_size = distance(first, last);
        this->_head = 0;
        this->_tail = 0;
        this->_size = distance(first, last);
        this->_buffer = _temp;
    }
    
    /**
     * @brief Assign an initializer list to a ring.
     * @param ilist The initializer list.
     */
    constexpr void assign(initializer_list<value_type> ilist){
        pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, ilist.size());
        uninitialized_copy(ilist.begin(), ilist.end(), _temp);
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
        this->_max_size = ilist.size();
        this->_head = 0;
        this->_tail = 0;
        this->_size = ilist.size();
        this->_buffer = _temp;
    }
    
    /**
     * @brief Assign a range to the ring.
     * @param rg A std::ranges::input_range container compatible range.
     */
    template<class R>
    constexpr void assign_range(R&& rg) requires(ranges::input_range<R> && convertible_to<ranges::range_reference_t<R>, T>){
        pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, ranges::distance(rg));
        uninitialized_copy(rg.begin(), rg.end(), _temp);
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
        this->_max_size = ranges::distance(rg);
        this->_head = 0;
        this->_tail = 0;
        this->_size = ranges::distance(rg);
        this->_buffer = _temp;
    }

    /**
     * @brief Get the allocator object.
     * 
     * @return allocator_type
     */
    constexpr allocator_type get_allocator() const noexcept{
        return this->_alloc;
    }

    /**
     * @brief Provides access to the data contained in the ring, starting from the tail.
     * @param pos The index of the element to obtain, counted from the tail of the ring.
     * @return reference 
     * @throw std::out_of_range If @a pos is an invalid index.
     */
    constexpr reference at(size_type pos){
        this->_M_range_check(pos);
        return this->_buffer[(this->_tail + pos) % this->_max_size];
    }
    
    /**
     * @brief Provides read-only access to the data contained in the ring, starting from the tail.
     * @param pos The index of the element to obtain, counted from the tail of the ring.
     * @return const_reference 
     * @throw std::out_of_range If @a pos is an invalid index.
     */
    constexpr const_reference at(size_type pos) const{
        this->_M_range_check(pos);
        return this->_buffer[(this->_tail + pos) % this->_max_size];
    }
    
    /**
     * @brief Provides access to the data contained in the ring, starting from the tail.
     * @param pos The index of the element to obtain, counted from the tail of the ring.
     * @return reference 
     * @throw std::out_of_range If @a pos is an invalid index.
     */
    constexpr reference operator[](size_type pos){
        this->_M_range_check(pos);
        return this->_buffer[(this->_tail + pos) % this->_max_size];
    }
    
    /**
     * @brief Provides read-only access to the data contained in the ring, starting from the tail.
     * @param pos The index of the element to obtain, counted from the tail of the ring.
     * @return const_reference 
     * @throw std::out_of_range If @a pos is an invalid index.
     */
    constexpr const_reference operator[](size_type pos) const{
        this->_M_range_check(pos);
        return this->_buffer[(this->_tail + pos) % this->_max_size];
    }
    
    /**
     * @brief Returns a reference to the data at the tail of the ring. Equivalent to *begin().
     * @return reference 
     */
    constexpr reference front(){
        return this->_buffer[this->_tail];
    }
    
    /**
     * @brief Returns a read-only reference to the data at the tail of the ring. Equivalent to *begin().
     * @return reference 
     */
    constexpr const_reference front() const{
        return this->_buffer[this->_tail];
    }
    
    /**
     * @brief Returns a reference to the data at the head of the ring. Equivalent to *(end() - 1).
     * @return reference 
     */
    constexpr reference back(){
        return this->_buffer[this->_head == 0 ? this->_max_size - 1 : this->_head - 1];
    }
    
    /**
     * @brief Returns a read-only reference to the data at the head of the ring. Equivalent to *(end() - 1).
     * @return reference 
     */
    constexpr const_reference back() const{
        return this->_buffer[this->_head == 0 ? this->_max_size - 1 : this->_head - 1];
    }
    
    /**
     * @brief Provides access to the internal contiguous array.
     * @return pointer 
     */
    constexpr pointer data() noexcept{
        return this->_buffer;
    }
    
    /**
     * @brief Provides read-only access to the internal contiguous array.
     * @return const_pointer 
     */
    constexpr const_pointer data() const noexcept{
        return this->_buffer;
    }

    /**
     * @brief Returns an iterator to the normal beginning (tail) of the ring.
     * @return iterator 
     */
    constexpr iterator begin() noexcept{
        if(this->_full()){
            return iterator(this, &this->_buffer[this->_tail], this->_tail, this->_max_size);
        } else if(this->empty()){
            return iterator(this, &this->_buffer[0], this->_tail, 0);
        } else{
            return iterator(this, &this->_buffer[this->_tail], this->_tail, this->_size);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the normal beginning (tail) of the ring.
     * @return const_iterator 
     */
    constexpr const_iterator begin() const noexcept{
        if(this->_full()){
            return const_iterator(this, &this->_buffer[this->_tail], this->_tail, this->_max_size);
        } else if(this->empty()){
            return const_iterator(this, &this->_buffer[0], this->_tail, 0);
        } else{
            return const_iterator(this, &this->_buffer[this->_tail], this->_tail, this->_size);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the normal beginning (tail) of the ring.
     * @return const_iterator 
     */
    constexpr const_iterator cbegin() const noexcept{
        if(this->_full()){
            return const_iterator(this, &this->_buffer[this->_tail], this->_tail, this->_max_size);
        } else if(this->empty()){
            return const_iterator(this, &this->_buffer[0], this->_tail, 0);
        } else{
            return const_iterator(this, &this->_buffer[this->_tail], this->_tail, this->_size);
        }
    }

    /**
     * @brief Returns an iterator to the normal end (head) of the ring. When the ring is full or empty, will be a sentinel.
     * @return iterator 
     */
    constexpr iterator end() noexcept{
        if(this->_full() || this->empty()){
            return iterator(this, &this->_buffer[-1], this->_head, 0);
        } else{
            return iterator(this, &this->_buffer[this->_head], this->_head, 0);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the normal end (head) of the ring. When the ring is full or empty, will be a sentinel.
     * @return const_iterator 
     */
    constexpr const_iterator end() const noexcept{
        if(this->_full() || this->empty()){
            return const_iterator(this, &this->_buffer[-1], this->_head, 0);
        } else{
            return const_iterator(this, &this->_buffer[this->_head], this->_head, 0);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the normal end (head) of the ring. When the ring is full or empty, will be a sentinel.
     * @return const_iterator 
     */
    constexpr const_iterator cend() const noexcept{
        if(this->_full() || this->empty()){
            return const_iterator(this, &this->_buffer[-1], this->_head, 0);
        } else{
            return const_iterator(this, &this->_buffer[this->_head], this->_head, 0);
        }
    }

    /**
     * @brief Returns an iterator to the reversed beginning (head - 1) of the ring.
     * @return reverse_iterator 
     */
    constexpr reverse_iterator rbegin() noexcept{
        if(this->_full()){
            size_type loc = this->_head == 0 ? this->_max_size - 1 : this->_head - 1;
            return reverse_iterator(this, &this->_buffer[loc], loc, this->_max_size);
        } else if(this->empty()){
            return reverse_iterator(this, &this->_buffer[0], this->_head, 0);
        } else{
            size_type loc = this->_head == 0 ? this->_max_size - 1 : this->_head - 1;
            return reverse_iterator(this, &this->_buffer[loc], loc, this->_size);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the reversed beginning (head - 1) of the ring.
     * @return const_reverse_iterator 
     */
    constexpr const_reverse_iterator rbegin() const noexcept{
        if(this->_full()){
            size_type loc = this->_head == 0 ? this->_max_size - 1 : this->_head - 1;
            return const_reverse_iterator(this, &this->_buffer[loc], loc, this->_max_size);
        } else if(this->empty()){
            return const_reverse_iterator(this, &this->_buffer[0], this->_head, 0);
        } else{
            size_type loc = this->_head == 0 ? this->_max_size - 1 : this->_head - 1;
            return const_reverse_iterator(this, &this->_buffer[loc], loc, this->_size);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the reversed beginning (head - 1) of the ring.
     * @return const_reverse_iterator 
     */
    constexpr const_reverse_iterator crbegin() const noexcept{
        if(this->_full()){
            size_type loc = this->_head == 0 ? this->_max_size - 1 : this->_head - 1;
            return const_reverse_iterator(this, &this->_buffer[loc], loc, this->_max_size);
        } else if(this->empty()){
            return const_reverse_iterator(this, &this->_buffer[0], this->_head, 0);
        } else{
            size_type loc = this->_head == 0 ? this->_max_size - 1 : this->_head - 1;
            return const_reverse_iterator(this, &this->_buffer[loc], loc, this->_size);
        }
    }

    /**
     * @brief Returns an iterator to the reversed end (tail - 1) of the ring. When the ring is full or empty, will be a sentinel.
     * 
     * @return reverse_iterator 
     */
    constexpr reverse_iterator rend() noexcept{
        size_type loc = this->_tail == 0 ? this->_max_size - 1 : this->_tail - 1;
        if(this->_full() || this->empty()){
            return reverse_iterator(this, &this->_buffer[-1], loc, 0);
        } else{
            return reverse_iterator(this, &this->_buffer[loc], loc, 0);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the reversed end (tail - 1) of the ring. When the ring is full or empty, will be a sentinel.
     * 
     * @return const_reverse_iterator 
     */
    constexpr const_reverse_iterator rend() const noexcept{
        size_type loc = this->_tail == 0 ? this->_max_size - 1 : this->_tail - 1;
        if(this->_full() || this->empty()){
            return const_reverse_iterator(this, &this->_buffer[-1], loc, 0);
        } else{
            return const_reverse_iterator(this, &this->_buffer[loc], loc, 0);
        }
    }
    
    /**
     * @brief Returns a read-only iterator to the reversed end (tail - 1) of the ring. When the ring is full or empty, will be a sentinel.
     * 
     * @return const_reverse_iterator 
     */
    constexpr const_reverse_iterator crend() const noexcept{
        size_type loc = this->_tail == 0 ? this->_max_size - 1 : this->_tail - 1;
        if(this->_full() || this->empty()){
            return const_reverse_iterator(this, &this->_buffer[-1], loc, 0);
        } else{
            return const_reverse_iterator(this, &this->_buffer[loc], loc, 0);
        }
    }

    /**
     * @brief Returns true if the ring is empty.
     * @return true 
     * @return false 
     */
    constexpr bool empty() const{
        return this->_size == 0;
    }
    
    /**
     * @brief Returns the current number of elements that are stored.
     * @return size_type 
     */
    constexpr size_type size() const noexcept{
        return this->_size;
    }
    
    /**
     * @brief Returns the maximum number of elements the ring can hold.
     * @return size_type 
     */
    constexpr size_type max_size() const noexcept{
        return this->_max_size;
    }
    
    /**
     * @brief Attempts to shrink the ring such that it now only has memory for the total number of stored elements, i.e. to reduce max_size() to size().
     */
    constexpr void shrink_to_fit(){
        if(this->_full()) return;
        pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, this->_size);
        uninitialized_move(this->begin(), this->end(), _temp);
        destroy(this->begin(), this->end());
        allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
        this->_buffer = _temp;
        this->_max_size = this->_size;
        this->_head = 0;
        this->_tail = 0;
    }

    /**
     * @brief Erases all elements. If the elements are pointers, will not erase the pointed-to memory.
     */
    constexpr void clear() noexcept{
        destroy(this->begin(), this->end());
        this->_size = 0;
        this->_head = 0;
        this->_tail = 0;
    }    
    
    /**
     * @brief Copies a value into the head of the ring.
     * @param value Data to be copied.
     */
    constexpr void push_back(const T& value){
        allocator_traits<Allocator>::destroy(this->_alloc, this->_buffer + this->_head);
        allocator_traits<Allocator>::construct(this->_alloc, this->_buffer + this->_head, value);
        this->_incr();
    }
    
    /**
     * @brief Moves a value into the head of the ring.
     * @param value Data to be moved.
     */
    constexpr void push_back(T&& value) noexcept{
        allocator_traits<Allocator>::destroy(this->_alloc, this->_buffer + this->_head);
        allocator_traits<Allocator>::construct(this->_alloc, this->_buffer + this->_head, value);
        this->_incr();
    }
    
    /**
     * @brief Constructs elements to the head of the ring. The ring will automatically roll over old data, and will not expand to accommodate the new data.
     * @param args Arguments to forward to the constructor of the element.
     * @return reference
     */
    template<class...Args> 
    constexpr reference emplace_back(Args&&... args){
        allocator_traits<Allocator>::destroy(this->_alloc, this->_buffer + this->_head);
        allocator_traits<Allocator>::construct(this->_alloc, this->_buffer + this->_head, args...);
        reference _data = *(this->_buffer + this->_head);
        this->_incr();
        return _data;
    }
    
    /**
     * @brief Appends a container compatible range to the head of the ring. If the current number of elements plus the number of elements in the range
     *  exceed max_size(), this operation will automatically expand the ring to accommodate the data.
     * @param rg std::ranges::input_range container compatible range to append.
     */
    template<class R>
    constexpr void append_range(R&& rg) requires(ranges::input_range<R> && convertible_to<ranges::range_reference_t<R>, T>){
        if(this->_size + ranges::distance(rg) <= this->_max_size){
            uninitialized_copy(rg.begin(), rg.end(), this->_buffer + this->_head);
            this->_head += ranges::distance(rg);
            this->_size += ranges::distance(rg);
        } else{
            size_type _new_max = this->_size + ranges::distance(rg);
            pointer _temp = allocator_traits<Allocator>::allocate(this->_alloc, _new_max);
            uninitialized_copy(this->begin(), this->end(), _temp);
            uninitialized_copy(rg.begin(), rg.end(), _temp + this->_size);
            destroy(this->begin(), this->end());
            allocator_traits<Allocator>::deallocate(this->_alloc, this->_buffer, this->_max_size);
            this->_buffer = _temp;
            this->_max_size = _new_max;
            this->_size = _new_max;
            this->_head = 0;
            this->_tail = 0;
        }
    }

    /**
     * @brief Destroys the tail element of the ring. If empty() is @a true, the behavior is undefined, though the ring will perform no operation on the underlying data.
     */
    constexpr void pop_front(){
        if(this->empty()) return;
        allocator_traits<Allocator>::destroy(this->_alloc, this->_buffer + this->_tail);
        --this->_size;
        this->_tail = (this->_tail + 1) % this->_max_size;
    }
    
    /**
     * @brief Resizes the ring to contain count number of elements, starting from the tail.
     * @param count The total number of elements to resize to.
     */
    constexpr void resize(size_type count){
        if(this->_size > count){
            for(size_type i = 0; i < this->_size - count; ++i){
                allocator_traits<Allocator>::destroy(this->_alloc, &*(this->rbegin() + i));
            }
            this->_head = (this->_tail + count) % this->_max_size;
            this->_size = count;
        } else if(this->_size < count){
            this->_M_resize_check(count);
            for(size_type i = 0; i < count - this->_size; ++i){
                allocator_traits<Allocator>::construct(this->_alloc, &*(this->end() + i));
            }
            this->_head = (this->_head + (count - this->_size)) % this->_max_size;
            this->_size = count;
        } else return;
    }
    
    /**
     * @brief Resizes the ring to contain count number of elements. If max_size() is less than count, additional copies of value are appended.
     * @param count The total number of elements to resize to.
     * @param value The value to initialize new elements with.
     */
    constexpr void resize(size_type count, const value_type& value){
        if(this->_size > count){
            for(size_type i = 0; i < this->_size - count; ++i){
                allocator_traits<Allocator>::destroy(this->_alloc, &*(this->rbegin() + i));
            }
            this->_head = (this->_tail + count) % this->_max_size;
            this->_size = count;
        } else if(this->_size < count){
            this->_M_resize_check(count);
            for(size_type i = 0; i < count - this->_size; ++i){
                allocator_traits<Allocator>::construct(this->_alloc, &*(this->end() + i), value);
            }
            this->_head = (this->_head + (count - this->_size)) % this->_max_size;
            this->_size = count;
        } else return;
    }
    
    /**
     * @brief Swaps the contents and capacity of two rings. Does not invoke any move, copy, or swap operations of the individual elements.
     * @param other A ring of identical element and allocator types.
     */
    constexpr void swap(ring& other) noexcept(allocator_traits<Allocator>::is_always_equal::value){
        using std::swap;
        swap(this->_max_size, other._max_size);
        swap(this->_head, other._head);
        swap(this->_tail, other._tail);
        swap(this->_size, other._size);
        swap(this->_buffer, other._buffer);
        if(allocator_traits<Allocator>::propagate_on_container_swap::value) swap(this->_alloc, other._alloc);
    }

    /**
     * @brief Swaps the contents and capacity of two rings. Does not invoke any move, copy, or swap operations of the individual elements.
     * 
     * @param lhs Container to swap.
     * @param rhs Container to swap.
     */
    constexpr friend void swap(ring& lhs, ring& rhs) noexcept(noexcept(lhs.swap(rhs))){
        lhs.swap(rhs);
    }

    /**
     * @brief Lexicographically compares two rings.
     * @param other A ring of identical element and allocator types.
     * @return true 
     * @return false 
     */
    constexpr bool operator==(const ring& other) const noexcept requires(equality_comparable<T>){
        return equal(this->begin(), this->end(), other.begin(), other.end());
    }
    /**
     * @brief Lexicographically three way compares two rings.
     * @param other A ring of identical element and allocator types.
     */
    constexpr strong_ordering operator<=>(const ring& other) const noexcept requires(equality_comparable<T>){
        return lexicographical_compare_three_way(this->begin(), this->end(), other.begin(), other.end());
    }
};
}
#endif
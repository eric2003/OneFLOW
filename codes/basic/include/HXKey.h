/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#pragma once
#include "HXVector.h"
#include <algorithm>

BeginNameSpace( ONEFLOW )

template<typename T = int>
class HXKey
{
public:
    using value_type     = T;
    //using container_type = std::vector<T>;
    using container_type = HXVector<T>;

    container_type data;   // 始终保持有序

    HXKey() = default;

    explicit HXKey(const container_type& nodes)
        : data(nodes)
    {
        Normalize();
    }

    explicit HXKey(container_type&& nodes)
        : data(std::move(nodes))
    {
        Normalize();
    }

    // 方便两点边
    //HXKey(T a, T b)
    //{
    //    data = {a, b};
    //    Normalize();
    //}

    // 两点快捷构造（对边很有用）
    HXKey(const T& a, const T& b)
    {
        data.resize(2);
        data[0] = a;
        data[1] = b;
        Normalize();
    }

    void Normalize()
    {
        std::sort(data.begin(), data.end());
        // 如果业务上不允许重复节点，可打开下一行
        // data.erase(std::unique(data.begin(), data.end()), data.end());
    }

    bool operator<(const HXKey& rhs) const
    {
        return data < rhs.data;
    }

    bool operator==(const HXKey& rhs) const
    {
        return data == rhs.data;
    }
};

//template<typename T = int>
//class HXKey
//{
//public:
//    using value_type     = T;
//    using container_type = std::vector<T>;
//    container_type data;   // 始终保持有序
//public:
//    std::size_t size, id;
//    HXVector< T > data;
//public:
//    HXKey();
//    HXKey( std::size_t size, std::size_t id = 0 );
//    HXKey( const HXKey & rhs );
//    HXKey & operator = ( const HXKey & rhs );
//    ~HXKey( void  );
//public:
//    bool operator < ( const HXKey & rhs ) const;
//};
//
//template < typename T >
//HXKey<T>::HXKey()
//{
//    this->size = 0;
//    this->id   = 0;
//}
//
//template < typename T >
//HXKey<T>::HXKey( std::size_t size, std::size_t id )
//{
//    this->size = size;
//    this->id   = id;
//    this->data.resize( size );
//}
//
//template < typename T >
//HXKey<T>::HXKey( const HXKey<T> & rhs )
//{
//    this->size = rhs.size;
//    this->id   = rhs.id;
//    this->data.resize( size );
//
//    for ( std::size_t i = 0; i < size; ++ i )
//    {
//        this->data[ i ] = rhs.data[ i ];
//    }
//}
//
//template < typename T >
//HXKey<T> & HXKey<T>::operator = ( const HXKey<T> & rhs )
//{
//    if ( this == & rhs ) return * this;
//
//    if ( this->size != rhs.size )
//    {
//        this->size = rhs.size;
//        this->data.resize( this->size );
//    }
//
//    this->id   = rhs.id;
//    for ( std::size_t i = 0; i < size; ++ i )
//    {
//        this->data[ i ] = rhs.data[ i ];
//    }
//
//    return * this;
//}
//
//template < typename T >
//HXKey<T>::~HXKey()
//{
//}
//
//template < typename T >
//bool HXKey<T>::operator < ( const HXKey<T> & rhs ) const
//{
//    if ( this->size != rhs.size ) return this->size < rhs.size;
//    for ( std::size_t i = 0; i < this->size; ++ i )
//    {
//        if ( this->data[ i ] != rhs.data[ i ] )
//        {
//            return this->data[ i ] < rhs.data[ i ];
//        }
//    }
//    return false;
//}

EndNameSpace

#ifndef PC_H
#define PC_H 1
#include <functional>
#include <ext/hash_map>

namespace __gnu_cxx
{
        template<> struct hash< std::string >
        {
                size_t operator()( const std::string& x ) const
                {
                        return hash< const char* >()( x.c_str() );
                }
        };
}

#endif

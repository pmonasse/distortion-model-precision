set(SRC ccmath_svd.cpp
        homography.cpp homography.h
        matrix.cpp matrix.h
        numerics.cpp numerics.h
        vector.cpp)

add_library(Numerics ${SRC})

if(UNIX)
    set_target_properties(Numerics PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
endif()

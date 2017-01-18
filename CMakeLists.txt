cmake_minimum_required(VERSION 2.8)

add_subdirectory(src/libNumerics)

set(SRC
    src/cmdLine.h
    src/model.cpp src/model.h
    src/modelLensfun.cpp src/modelLensfun.h
    src/modelDivision.cpp src/modelDivision.h
    src/modelFOV.cpp src/modelFOV.h
    src/modelPoly.cpp src/modelPoly.h
    src/modelRadial.cpp src/modelRadial.h
    src/modelRational.cpp src/modelRational.h
    src/main.cpp)

add_executable(precision_analysis ${SRC})
target_link_libraries(precision_analysis Numerics)

if(UNIX)
    set_target_properties(precision_analysis PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
endif()

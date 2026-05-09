cmake_minimum_required(VERSION 3.10)

# Define the project
project(CB2 VERSION 1.0)

# Define C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Source files
set(SOURCES
    boost_1.cpp
    global_functions.cpp
    main.cpp
    model1.cpp
    model2.cpp
    model2b.cpp
)

# Header files
set(HEADERS
    boost_1.h
    global_functions.h
    global_variables.h
    model1.h
    model2.h
    model2b.h
)

# Add the executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Set configurations for Win32 and x64
set(CPLEX_INCLUDE_DIR "C:/Program Files/IBM/ILOG/CPLEX_Studio2211") 
set(BOOST_INCLUDE_DIR "C:/source/boost/boost_1_82_0")
set(CPLEX_LIB_DIR "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/lib/x64_windows_msvc14/stat_mda")
set(CONCERT_LIB_DIR "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/concert/lib/x64_windows_msvc14/stat_mda")

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${BOOST_INCLUDE_DIR}
    ${CPLEX_INCLUDE_DIR}/cplex/include
    ${CPLEX_INCLUDE_DIR}/concert/include
)

# Link CPLEX libraries
target_link_directories(${PROJECT_NAME} PRIVATE ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE cplex2211 ilocplex concert)

# Debug configuration
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MD /D _CONSOLE /D WIN32 /D IL_STD /D _ALLOW_ITERATOR_DEBUG_LEVEL_MISMATCH /D _CRT_SECURE_NO_WARNINGS /W3 /Zc:__cplusplus")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /DEBUG")

# Release configuration
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /O2 /D NDEBUG /D _CONSOLE /W3")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /OPT:REF /OPT:ICF")

# Specify target properties for Debug and Release configurations
if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:DEBUG>:/Zi /SDL>")
    target_compile_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:RELEASE>:/GL /Oi>")
    target_link_options(${PROJECT_NAME} PRIVATE "$<$<CONFIG:RELEASE>:/LTCG>")
endif()

# Link incrementally in Debug, disable in Release
set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS_DEBUG "/INCREMENTAL:YES"
    LINK_FLAGS_RELEASE "/INCREMENTAL:NO"
)
if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
endif()
# Set up for a Console application
set_target_properties(${PROJECT_NAME} PROPERTIES
    WIN32_EXECUTABLE OFF
)

# Additional dependencies for Debug|x64
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        cplex2211.lib
        ilocplex.lib
        concert.lib
    )
endif()

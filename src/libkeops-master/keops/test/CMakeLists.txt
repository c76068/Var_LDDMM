find_package(Threads)
find_package(GTest)

if(Threads_FOUND) # skip test if pthread is not found
    if(NOT GTest_FOUND)

        # Tell the user we are going to d/l GTest
        message(STATUS "Use ExternalProject to get the missing GTest lib. It may increase compilation time.")

        # Enable ExternalProject CMake module
        INCLUDE(ExternalProject)
    
        # Set default ExternalProject root directory
        SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/third_party)
    
        # Add gtest
        # http://stackoverflow.com/questions/9689183/cmake-googletest
        ExternalProject_Add(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            INSTALL_COMMAND ""
            # Wrap download, configure and build steps in a script to log output
            LOG_DOWNLOAD ON
            LOG_CONFIGURE ON
            LOG_BUILD ON
        )
    
        ExternalProject_Get_Property(googletest source_dir)
        set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)
        Include_Directories(${GTEST_INCLUDE_DIRS})
    
        ExternalProject_Get_Property(googletest binary_dir)
        set(GTEST_LIBRARY_PATH ${binary_dir}/googlemock/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
        set(GTEST_LIBRARY gtest)
        add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
        set_target_properties(${GTEST_LIBRARY} PROPERTIES
            IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
        )
    
        add_dependencies(${GTEST_LIBRARY} googletest)
    
        set(GTEST_MAIN_LIBRARY_PATH ${binary_dir}/googlemock/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
        set(GTEST_MAIN_LIBRARY gtest_main)
        add_library(${GTEST_MAIN_LIBRARY} UNKNOWN IMPORTED)
        set_target_properties(${GTEST_MAIN_LIBRARY} PROPERTIES
            IMPORTED_LOCATION ${GTEST_MAIN_LIBRARY_PATH}
        )
        add_dependencies(${GTEST_MAIN_LIBRARY} googletest)
    
    endif()
    
    
    ### - unit test
    if (USE_CUDA)
        CUDA_add_executable(
            gtest_grad1conv
            ${CMAKE_CURRENT_SOURCE_DIR}/gtest_grad1conv.cu
        )
    
        target_link_libraries(
            gtest_grad1conv
            radial_kernels_grad1conv
            gtest
            pthread
        )
    
        #--------------------------------------------------------
        Add_library(
            testic OBJECT 
            ${CMAKE_CURRENT_SOURCE_DIR}/gtest_conv.cpp
        )
    
        target_compile_options(
            testic BEFORE
            PRIVATE -include formula.h
        )
    
        Add_dependencies(
            testic
            gtest
        )

        Add_executable(
            gtest_conv
            $<TARGET_OBJECTS:testic>
        )
    
        target_link_libraries(
            gtest_conv
            keops
            gtest
            pthread
        )
    
        set_target_properties( # pass needed options to add_executable()
            gtest_conv PROPERTIES
            LINKER_LANGUAGE CXX
        )
    
    endif()

endif() # Threads_FOUND


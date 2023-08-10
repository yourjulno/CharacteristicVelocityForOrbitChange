add_subdirectory (googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})


file(GLOB_RECURSE TEST_SOURCES test_*.cpp)  # В переменной TEST_SOURCES лежат все тесты, имена файлов которые начинаются с test_


foreach(TEST_FILE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    set(EXECUTABLE_NAME run_${TEST_FILE_NAME})
    #    set(GCC_COMPILE_RELEASE_OPTIONS "-O3")
    add_executable(${EXECUTABLE_NAME} ${TEST_FILE}

    )
    target_link_options(${EXECUTABLE_NAME} PRIVATE -static-libstdc++)
    target_link_libraries(
            ${EXECUTABLE_NAME}
            CharacteristicVelocityForOrbitChange	# Подключаем проект
            gtest		# Подключаем библиотеку gtest
            gtest_main	# Подключаем библиотеку gtest
    )
    #    target_compile_options(${EXECUTABLE_NAME} -O3)
    message("TEST_FILE: ${TEST_FILE}")
    add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})
endforeach()
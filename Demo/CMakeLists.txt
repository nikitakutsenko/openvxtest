find_package(OpenCV REQUIRED PATHS "/home/nikita/OpenVX_MIET/vcpkg/installed/x64-linux/share/opencv/")


file(GLOB SRCS Kernels/*.cpp)
add_executable(demo demo.cpp ${SRCS})
target_link_libraries(demo openvx_kernels ${OpenCV_LIBS})

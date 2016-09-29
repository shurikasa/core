find_package(PETSc REQUIRED)

set(SOURCES
app.cc
assemble.cc
integrate.cc
linsys.cc
post.cc
pre.cc
utils.cc
)


set(HEADERS
app.h
integrate.h
linsys.h
utils.h
)

include_directories(${CMAKE_SOURCE_DIR}/apf ${CMAKE_SOURCE_DIR}/can ${CMAKE_SOURCE_DIR}/gmi ${CMAKE_SOURCE_DIR}/pcu ${PETSC_INCLUDES})

add_library(pe ${SOURCES})

target_include_directories(pe PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

add_executable(pe_exec main.cc)
target_link_libraries(pe_exec core pe ${PETSC_LIBRARIES})

bob_export_target(pe_exec)
bob_end_subdir()


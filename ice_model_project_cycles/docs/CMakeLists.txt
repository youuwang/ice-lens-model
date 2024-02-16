
find_package(Doxygen
  REQUIRED dot)


if (DOXYGEN_FOUND)
  # to set other options, read: https://cmake.org/cmake/help/v3.9/module/FindDoxygen.html

  set(DOXYGEN_INPUT_DIR ${CMAKE_SOURCE_DIR}/src)
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
  set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  set(DOXYGEN_XML_DIR ${DOXYGEN_OUTPUT_DIR}/xml)
  set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  #Replace variables inside @@ with the current values
  configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
  add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
    DEPENDS ${DOXYGEN_INPUT_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
    MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
    COMMENT "Generating documention with Doxygen")
  
  add_custom_target(doc ALL DEPENDS ${DOXYGEN_INDEX_FILE})
endif(DOXYGEN_FOUND)


find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

add_custom_target(Sphinx ALL
  COMMAND ${SPHINX_EXECUTABLE} -b html
  -Dbreathe_projects.ice_model=${DOXYGEN_XML_DIR}
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS doc
  COMMENT "Generating documentation with Sphinx")

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(VisualImageRegistration)

set(Slicer_DIR T:/Slicer-build/Slicer-build)
set(BRAINSCommonLib_SRC  ${Slicer_DIR}/../BRAINSTools/BRAINSCommonLib)
set(BRAINSCommonLib_DIR  ${Slicer_DIR}/E/BRAINSTools/BRAINSCommonLib)
#set(BRAINSFitLib_SRC  ${Slicer_DIR}/../BRAINSTools/BRAINSFit)
#set(BRAINSFitLib_DIR  ${Slicer_DIR}/E/BRAINSTools/BRAINSFit)
set(ITK_DIR ${Slicer_DIR}/../ITK-build)



FIND_PACKAGE(ITK)

IF(ITK_FOUND) 
    INCLUDE(${ITK_USE_FILE})

ELSE(ITK_FOUND) 
    MESSAGE(FATAL_ERROR 
            "ITK not found. Please set ITK_DIR.") 
ENDIF(ITK_FOUND)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)



include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${BRAINSCommonLib_SRC}
  ${BRAINSCommonLib_DIR}
  #${BRAINSFitLib_SRC}
  #${BRAINSFitLib_DIR}
  )
  message(${CMAKE_CURRENT_BINARY_DIR})




set(${PROJECT_NAME}_SRCS
BRAINSFit.cxx
#BRAINSFitHelper.h
#BRAINSFitHelper.cxx
#genericRegistrationHelper.h
#genericRegistrationHelper.hxx
#BRAINSCommonLibWin32Header.h
#BRAINSThreadControl.h
#BRAINSThreadControl.cxx
#BRAINSFitHelperTemplate.h
#BRAINSFitHelperTemplate.hxx
#itkIO.h
#Imgmath.h
)

set(lib_link
    ${ITK_LIBRARIES}
    #${Slicer_DIR}/lib/Visual3D-4.11/Release/BRAINSCommonLib.lib
	${Slicer_DIR}/lib/Visual3D-4.11/DeBug/BRAINSCommonLib.lib
)

ADD_EXECUTABLE(${PROJECT_NAME} ${${PROJECT_NAME}_SRCS} )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${lib_link})
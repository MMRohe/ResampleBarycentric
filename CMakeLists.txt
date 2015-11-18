PROJECT(ResampleBarycentric)

# Define the minimum CMake version needed
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(NOT VTK_DIR)
  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ENDIF(NOT VTK_DIR)
INCLUDE(${VTK_USE_FILE})


FIND_PACKAGE(TTK)
IF(TTK_FOUND)
INCLUDE(${TTK_USE_FILE})
ENDIF(TTK_FOUND)


INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/RegistrationProgrammingInterface
    ${PROJECT_SOURCE_DIR}/LogDemons
    ${PROJECT_SOURCE_DIR}/tclap

    ${PROJECT_SOURCE_DIR}/LCClogDemons
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

SET(ResampleBarycentric_Headers
  ${PROJECT_SOURCE_DIR}/LCClogDemons/itkMultiResolutionLCCDeformableRegistration.txx
  ${PROJECT_SOURCE_DIR}/LCClogDemons/itkMultiResolutionLCCDeformableRegistration.h
  ${PROJECT_SOURCE_DIR}/LogDemons/itkMultiResolutionLogDomainDeformableRegistration.hxx
  ${PROJECT_SOURCE_DIR}/LogDemons/itkMultiResolutionLogDomainDeformableRegistration.h
  ${PROJECT_SOURCE_DIR}/RegistrationProgrammingInterface/rpiRegistrationMethod.hxx
  ${PROJECT_SOURCE_DIR}/RegistrationProgrammingInterface/rpiRegistrationMethod.cxx
${PROJECT_SOURCE_DIR}/LCClogDemons/rpiLCClogDemons.hxx
${PROJECT_SOURCE_DIR}/LCClogDemons/rpiLCClogDemons.cxx
${PROJECT_SOURCE_DIR}/LCClogDemons/itkLCCDeformableRegistrationFilter.txx
${PROJECT_SOURCE_DIR}/LCClogDemons/itkLCCDeformableRegistrationFilter.h
${PROJECT_SOURCE_DIR}/LCClogDemons/itkSymmetricLCClogDemonsRegistrationFilter.h
${PROJECT_SOURCE_DIR}/LCClogDemons/itkSymmetricLCClogDemonsRegistrationFilter.txx
${PROJECT_SOURCE_DIR}/LCClogDemons/itkLocalCriteriaOptimizer.h
${PROJECT_SOURCE_DIR}/LCClogDemons/itkLocalCriteriaOptimizer.txx

  )

SET(ResampleBarycentric_Srcs
  ${PROJECT_SOURCE_DIR}/Main.cxx
)

ADD_EXECUTABLE( ResampleBarycentric ${ResampleBarycentric_Headers} ${ResampleBarycentric_Srcs})



Set( VTK_DEPENDANCIES "" )
If( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        vtkCommon
        vtkGraphics
        vtkIO
        vtkRendering
        vtkFiltering
    )
Else( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        ${VTK_LIBRARIES}
    )
Endif( "${VTK_MAJOR_VERSION}" LESS 6 )


Target_Link_Libraries(ResampleBarycentric
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
  ${TTK_LIBRARIES}

)

# Dependencies (only for VTK < 6) :
If( "${VTK_MAJOR_VERSION}" LESS 6 )

Add_Dependencies(ResampleBarycentric
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Endif( "${VTK_MAJOR_VERSION}" LESS 6 )

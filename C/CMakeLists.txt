cmake_minimum_required(VERSION 3.12)
project(CATCH22 C)

set(CMAKE_C_STANDARD 11)

# Add your source files
set(SOURCE_FILES
    main.c
    CO_AutoCorr.c
    DN_HistogramMode_10.c
    DN_HistogramMode_5.c
    DN_OutlierInclude.c
    FC_LocalSimple.c
    IN_AutoMutualInfoStats.c
    MD_hrv.c
    PD_PeriodicityWang.c
    SB_BinaryStats.c
    SB_CoarseGrain.c
    SB_MotifThree.c
    SB_TransitionMatrix.c
    SC_FluctAnal.c
    SP_Summaries.c
    DN_Mean.c
    DN_Spread_Std.c
    butterworth.c
    fft.c
    helper_functions.c
    histcounts.c
    splinefit.c
    stats.c
)

# Add the executable
add_executable(run_features ${SOURCE_FILES})

# Add OpenMP support
if(OpenMP_CXX_FOUND)
    target_link_libraries(run_features PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add FFTW library
find_package(FFTW3 REQUIRED)
target_include_directories(run_features PRIVATE ${FFTW3_INCLUDE_DIRS})
target_link_libraries(run_features PRIVATE ${FFTW3_LIBRARIES})

# Specify include directories
target_include_directories(run_features PRIVATE /usr/local/include)

# Specify library directories
target_link_directories(run_features PRIVATE /usr/local/lib)

# Specify optimization level
target_compile_options(run_features PRIVATE -O3)

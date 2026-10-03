# Enable the system UCRT integration by default on supported toolchains.
set(_astraxis_vc_ltl_default OFF)
if(WIN32 AND MSVC)
    set(_astraxis_vc_ltl_default ON)
endif()
option(ASTRAXIS_USE_VC_LTL "Use VC-LTL5 with the Windows 10 system UCRT in non-Debug builds" ${_astraxis_vc_ltl_default})

if(NOT ASTRAXIS_USE_VC_LTL)
    return()
endif()

if(NOT WIN32 OR NOT MSVC)
    message(FATAL_ERROR "ASTRAXIS_USE_VC_LTL requires Windows and MSVC")
endif()

include(FetchContent)
FetchContent_Declare(vc_ltl
    URL https://github.com/Chuyu-Team/VC-LTL5/releases/download/v5.3.1/VC-LTL-Binary.7z
    URL_HASH SHA256=7a18799ed3aa84a225610a5447a56bc534c5c98ccb8dec05caba0e3f633431ad
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(vc_ltl)

set(VC_LTL_Root "${vc_ltl_SOURCE_DIR}")
set(WindowsTargetPlatformMinVersion "10.0.10240.0")
set(SupportLTL "true")
# Ask upstream to avoid applying unconditional directory-wide settings.
set(VC_LTL_EnableCMakeInterface ON)
include("${VC_LTL_Root}/config/config.cmake")

if(NOT InternalSupportLTL STREQUAL "ucrt")
    message(FATAL_ERROR "VC-LTL failed to select the Windows 10 UCRT runtime")
endif()

# Visual Studio generates Debug and Release together. Gate both headers and
# libraries per configuration, including dependencies created in subdirectories.
foreach(_vc_ltl_include IN LISTS VC_LTL_Include)
    include_directories(BEFORE SYSTEM "$<$<NOT:$<CONFIG:Debug>>:${_vc_ltl_include}>")
endforeach()
link_directories("$<$<NOT:$<CONFIG:Debug>>:${VC_LTL_Library}>")

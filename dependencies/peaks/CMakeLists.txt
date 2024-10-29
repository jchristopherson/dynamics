include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    peaks
    GIT_REPOSITORY "https://github.com/jchristopherson/peaks"
    GIT_TAG main
)

FetchContent_MakeAvailable(peaks)

set(peaks_INCLUDE_DIR ${peaks_BINARY_DIR}/include)
set(peaks_INCLUDE_DIR ${peaks_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(peaks_LIBRARY peaks)
set(peaks_LIBRARY ${peaks_LIBRARY} PARENT_SCOPE)

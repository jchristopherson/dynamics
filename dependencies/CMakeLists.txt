# Get FERROR
add_subdirectory(ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)

# Get LINALG
add_subdirectory(linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)

# Get SPECTRUM
add_subdirectory(spectrum)
set(spectrum_LIBRARY ${spectrum_LIBRARY} PARENT_SCOPE)
set(spectrum_INCLUDE_DIR ${spectrum_INCLUDE_DIR} PARENT_SCOPE)

# Get FFTPACK
add_subdirectory(fftpack)
set(fftpack_LIBRARY ${fftpack_LIBRARY} PARENT_SCOPE)
set(fftpack_INCLUDE_DIR ${fftpack_INCLUDE_DIR} PARENT_SCOPE)

# Get DIFFEQ
add_subdirectory(diffeq)
set(diffeq_LIBRARY ${diffeq_LIBRARY} PARENT_SCOPE)
set(diffeq_INCLUDE_DIR ${diffeq_INCLUDE_DIR} PARENT_SCOPE)

# Get NONLIN
add_subdirectory(nonlin)
set(nonlin_LIBRARY ${nonlin_LIBRARY} PARENT_SCOPE)
set(nonlin_INCLUDE_DIR ${nonlin_INCLUDE_DIR} PARENT_SCOPE)

# Get FSTATS
add_subdirectory(fstats)
set(fstats_LIBRARY ${fstats_LIBRARY} PARENT_SCOPE)
set(fstats_INCLUDE_DIR ${fstats_INCLUDE_DIR} PARENT_SCOPE)

# Get PEAKS
add_subdirectory(peaks)
set(peaks_LIBRARY ${peaks_LIBRARY} PARENT_SCOPE)
set(peaks_INCLUDE_DIR ${peaks_INCLUDE_DIR} PARENT_SCOPE)

program prog081

#ifdef _OPENMP
    integer, external :: omp_get_num_threads
    external :: omp_set_num_threads
#endif
#ifdef _OPENMP
    print '("The number of threads is ",i0)', omp_get_num_threads()
    print '("Setting the number of threads to 4")'
    call omp_set_num_threads(4)
    print '("The number of threads is now ",i0)', omp_get_num_threads()
    ! Prints 1 because were are not in a parallel region.
#else
    print '("Non-parallel version")'
#endif

end program prog081
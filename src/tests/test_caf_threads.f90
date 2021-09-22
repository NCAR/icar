program test_caf_threads

#if defined(_OPENMP)
    use omp_lib
#endif

    implicit none

    integer :: tid


    !$omp parallel sections private(tid)

    !$omp section

#if defined(_OPENMP)
        tid = omp_get_thread_num()
#else
        tid = 1
#endif

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(), "    Thread:",tid, " is hanging out!"
        !$omp end critical(print_lock)
        call sleep(1)
    !$omp section

#if defined(_OPENMP)
        tid = omp_get_thread_num()
#else
        tid = 1
#endif

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(),"    Thread:",tid, " is shooting the breeze!"
        !$omp end critical(print_lock)
        call sleep(1)
    !$omp section

#if defined(_OPENMP)
        tid = omp_get_thread_num()
#else
        tid = 1
#endif

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(), "    Thread:",tid, " is waiting..."
        !$omp end critical(print_lock)
        call sleep(1)
        sync all

    !$omp section

#if defined(_OPENMP)
        tid = omp_get_thread_num()
#else
        tid = 1
#endif

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(), "    Thread:",tid, " is not waiting!"
        !$omp end critical(print_lock)
    !$omp end parallel sections

    print '(A,I2,A)', "Image",this_image(), " Completed"

end program test_caf_threads

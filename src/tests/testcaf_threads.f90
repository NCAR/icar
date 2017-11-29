program test_caf_threads

    use omp_lib
    implicit none

    !$omp parallel sections

    !$omp section

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(), "    Thread:",omp_get_thread_num(), " is hanging out!"
        !$omp end critical(print_lock)

    !$omp section

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(),"    Thread:",omp_get_thread_num(), " is shooting the breeze!"
        !$omp end critical(print_lock)

    !$omp section

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(), "    Thread:",omp_get_thread_num(), " is waiting..."
        !$omp end critical(print_lock)

        sync all

    !$omp section

        !$omp critical(print_lock)
            print '(A,I2,A,I2,A)', "Image:",this_image(), "    Thread:",omp_get_thread_num(), " is not waiting!"
        !$omp end critical(print_lock)

    !$omp end parallel sections

    print '(A,I2,A)', "Image",this_image(), " Completed"

end program test_caf_threads

        integer, intent(in) :: source_image, first_image, last_image
        integer :: dest_image

        !call co_broadcast(coarray, source_image)
        !return

        if (first_image==last_image) return

        if (source_image/=first_image) then
            dest_image=first_image

            if (this_image()==source_image) then
                ! This is a "put" and needs to come before the sync
                ! coarray(:,:,:,:)[dest_image] = coarray
                sync images(dest_image)
            elseif (this_image()== dest_image) then
                ! This is a "get" and needs to come after the sync
                sync images(source_image)

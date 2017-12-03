        endif

        if (this_image()<source_image) then
            call broadcast(coarray, dest_image, dest_image, source_image-1)
        else
            call broadcast(coarray, source_image, source_image, last_image)
        endif
        else
        dest_image = ((last_image-first_image)+1)/2 + first_image

        if (this_image()==source_image) then
            ! This is a "put" and needs to come before the sync
            ! coarray(:,:,:,:)[dest_image] = coarray
            sync images(dest_image)
        elseif (this_image()== dest_image) then
            ! This is a "get" and needs to come after the sync
            sync images(source_image)

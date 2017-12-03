        endif

        if (this_image()<dest_image) then
            call broadcast(coarray, source_image, first_image, dest_image-1)
        else
            call broadcast(coarray, dest_image, dest_image, last_image)
        endif
        endif

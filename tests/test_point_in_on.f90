program test_point_in_on

    use geo, only : point_is_on_line, point_in_poly
    implicit none

    real :: x,y
    real :: poly(2,4)
    real :: x0,y0,x1,y1

    logical :: passing
    passing=.True.

    x0=0; y0=0; x1=1; y1=1
    x=0.5; y=0.5
    if (point_is_on_line(x,y,x0,y0,x1,y1)) then
        print*, "Correctly identified ",x, y, "on line."
    else
        passing=.false.
        print*, "Point_is_on_line Failed", x, y
    endif

    x=0.6; y=0.5
    if (point_is_on_line(x,y,x0,y0,x1,y1)) then
        print*, "Point_is_on_line Failed", x, y
        passing=.false.
    else
        print*, "Correctly identified ",x, y, "off line."
    endif


    poly(1,:) = [0,0,1,1]
    poly(2,:) = [0,1,1,0]

    x=0.5; y=0.5
    passing = passing.and.test_poly(x,y,poly,.True.)

    x=1.01; y=0.5
    passing = passing.and.test_poly(x,y,poly,.False.)


    poly(1,:) = [0.5,0.0,0.5,1.0]
    poly(2,:) = [0.0,0.5,1.0,0.5]

    x=0.5; y=0.5
    passing = passing.and.test_poly(x,y,poly,.True.)
    x=0.51; y=0.51
    passing = passing.and.test_poly(x,y,poly,.True.)
    x=0.25; y=0.25
    passing = passing.and.test_poly(x,y,poly,.True.)
    x=0.2; y=0.25
    passing = passing.and.test_poly(x,y,poly,.False.)
    x=0.2; y=0.75
    passing = passing.and.test_poly(x,y,poly,.False.)
    x=0.25; y=0.75
    passing = passing.and.test_poly(x,y,poly,.True.)

    if (passing) then
        print*, "Passed"
    else
        print*, "Failed"
    endif

contains

    function test_poly(x,y,poly, expected) result(passed)
        implicit none
        real,       intent(in) :: x,y
        real,       intent(in) :: poly(:,:)
        logical,    intent(in) :: expected
        logical :: passed

        if (point_in_poly(x,y,poly) == expected) then
            passed=.True.
        else
            passed=.False.
            print*, "point_in_poly Failed", x, y
            print*, poly(1,:)
            print*, poly(2,:)
        endif

    end function test_poly

end program test_point_in_on

program test_point_in_on

    use geo, only : point_is_on_line, point_in_poly
    implicit none

    real :: x,y
    real :: poly(2,4), large_poly(2,8)
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

    print*, "Testing diamond polygon edge cases"
    print*, " x = polygon vertex "
    print*, " o = test point     "
    print*, ""

    x=0.5; y=0.5
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        /   \    "
    print*, "       x  o  x   "
    print*, "        \   /    "
    print*, "         \ /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.True.)

    x=-0.1; y=0.5
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        /   \    "
    print*, "   o   x     x   "
    print*, "        \   /    "
    print*, "         \ /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.False.)

    x=0.1; y=0.0
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        /   \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "         \ /     "
    print*, "   o      x      "
    passing = passing.and.test_poly(x,y,poly,.False.)

    x=0.51; y=0.51
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        / o \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "         \ /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.True.)

    x=0.25; y=0.25
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        /   \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "         o /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.True.)

    x=0.2; y=0.25
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        /   \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "        o\ /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.False.)

    x=0.9; y=0.25
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         / \     "
    print*, "        /   \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "         \ / o   "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.False.)

    x=0.2; y=0.75
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "        o/ \     "
    print*, "        /   \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "         \ /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.False.)

    x=0.25; y=0.75
    print*, ""
    print*, "Testing: ", x, y
    print*, "          x      "
    print*, "         o \     "
    print*, "        /   \    "
    print*, "       x     x   "
    print*, "        \   /    "
    print*, "         \ /     "
    print*, "          x      "
    passing = passing.and.test_poly(x,y,poly,.True.)

    print*, "Testing larger polygon edge cases"
    large_poly(1,:) = [0.0, 0.0, 1.0, 1.0, 0.75, 0.75, 0.25, 0.25]
    large_poly(2,:) = [0.0, 1.0, 1.0, 0.0, 0.00, 0.50, 0.50, 0.00]

    x=0.25; y=0.5
    print*, ""
    print*, "Testing: ", x, y
    print*, "  x-------------x"
    print*, "  |             |"
    print*, "  |  o x---x    |"
    print*, "  |    |   |    |"
    print*, "  x----x   x----x"
    passing = passing.and.test_poly(x,y,large_poly,.True.)


    print*, ""
    large_poly(1,:) = [0.0, 0.0, 1.0, 1.0,  0.75, 0.75, 0.25, 0.25]
    large_poly(2,:) = [0.0, 1.0, 1.0, 0.75, 0.75, 0.50, 0.50, 0.00]

    x=0.25; y=0.5
    print*, ""
    print*, "Testing: ", x, y
    print*, "  x--------------x"
    print*, "  |              |"
    print*, "  |         x----x"
    print*, "  |  o x----x     "
    print*, "  |    |          "
    print*, "  |    |          "
    print*, "  x----x          "

    passing = passing.and.test_poly(x,y,large_poly,.True.)

    print*, ""
    print*, "Final Result:"
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
            print*, "point_in_poly Failed for:"
            print*, " x=",x,"    y=",y
            print*, "Polygon x vertices:"
            print*, poly(1,:)
            print*, "Polygon y vertices:"
            print*, poly(2,:)
        endif

    end function test_poly

end program test_point_in_on

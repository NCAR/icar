program test_point_in_on

    use geo, only : point_is_on_line, point_in_poly
    implicit none

    real :: x,y
    real :: poly(2,4), large_poly(2,8), triangle(2,3)
    real :: x0,y0,x1,y1
    integer :: i,j

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


    x = 35.4495850
    y = 34.0271034
    print*, "Testing problem grid cell from a simulation."
    print*, "[",x,",", y,"]"
    print*, ""

    ! These are the coordinates of the bounding box of grid cells for which find_surrounding was failing before
    ! poly(1,:) = [35.40149, 35.46729, 35.40778, 35.47363]
    ! poly(2,:) = [33.99421, 33.98899, 34.04877, 34.04354]

    triangle(1,:) = [35.4077759, 35.4736328, 35.4375458]
    triangle(2,:) = [34.0487747, 34.0435371, 34.0188789]

    passing = passing.and.test_poly(x,y,triangle,.True.)

    ! Note that this will actually fail because the precision is "too" high,
    ! but if the points are in a different order it will fail for the other viable set as well! 1e-7 is too high a precision for lat/lon as single precision reals
    !
    ! print*, point_in_poly(x,y,triangle, precision=1e-7)
    ! result = .False. ^^^^

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


    large_poly(1,1:3) = [281.6548,281.7756,281.9876]
    large_poly(2,1:3) = [34.91632, 35.36251,35.08936]

    x=281.98; y=35.0
    print*, ""
    print*, "Testing: ", x, y
    print*, "  x-----x"
    print*, "  |    / "
    print*, "  |   /o "
    print*, "  |  /   "
    print*, "  | /    "
    print*, "  |/     "
    print*, "  x      "
    passing = passing.and.test_poly(x,y,large_poly(:,1:3),.False.)


    large_poly(1,1:3) = [281.6548,282.1,281.9876]
    large_poly(2,1:3) = [34.91632, 34.8,35.08936]

    x=281.98; y=35.0
    print*, ""
    print*, "Testing: ", x, y
    print*, "        x "
    print*, "       /| "
    print*, "      /o| "
    print*, "     /  | "
    print*, "    /   | "
    print*, "   /    | "
    print*, "  x-----x"
    passing = passing.and.test_poly(x,y,large_poly(:,1:3),.True.)

    large_poly(1,1:3) =[290.0750,289.8967,289.7062]
    large_poly(2,1:3) =[38.79078,38.35439,38.64164]
    x=289.9; y=38.72
    print*, ""
    print*, "Testing: ", x, y
    print*, "        x "
    print*, "       /| "
    print*, "      o | "
    print*, "     /  | "
    print*, "    /   | "
    print*, "   /    | "
    print*, "  x-----x"
    passing = passing.and.test_poly(x,y,large_poly(:,1:3),.True.)



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

        if (point_in_poly(x,y,poly) .eqv. expected) then
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

rei
fetch 2bl2, async=0
load c:/Users/Nick/Downloads/2bl2_2fofc.dsn6

hi ////LHG`
isosurface vm, 2bl2_2fofc, 0.5, ////LHG`, carve=2

remove ////UMQ+NA+HOH`
util.cbc //A+B+C+D+E+F+G+H+I+J

set_view (\
    -0.225055233,   -0.239427328,   -0.944470525,\
    -0.060426727,    0.970901608,   -0.231726825,\
     0.972467601,    0.004920267,   -0.232972577,\
    -0.000208192,   -0.000080936, -146.140426636,\
    32.263572693,   36.124916077,   34.739185333,\
    36.761730194,  255.552764893,  -20.000000000 )

set ray_opaque_background, 
ray 2000, 1600
png hw4.py
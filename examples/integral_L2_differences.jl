# integral L2 differences can be used for edge detection, and in general for detecting patches with a high variability. 

using TestImages, ImageView, MarcIntegralArrays

img = TestImages.testimage( "fabio_color_256" );

scale = (2,2)
L2s   = MarcIntegralArrays.localL2avg( img, scale );
ImageView.imshow( L2s );


# JPEG-Compression-in-C
Implementation of Standard JPEG Image Compression algorithm in C.

-> Input raw image named 'barbara_raw_input' is given as an input to the program.
<br>
-> Running the program will ask for input raw image name. Provide name 'barbara_raw_input'.
<br>
-> When asked for size of image, provide 512 two times seperating by space as our input raw image has dimension 512*512.Input takes number of columns and rows respectively.
<br>
-> Output reconstructed image is constructed and named 'barbara_raw_input_jpeg_reconstructed'.
<br>
-> The txt file named 'compressed_image' is the actual compressed image that we store after applying DCT(Discrete Cosine Transform).
<br>
-> Use software like 'imageJ' to open raw image and see results.
<br>
-> 'Output_Screenshot.png' shows input image(on left) and output image(on right) obtained after reconstructing from 'compressed_image' file.
<br>
-> You will be able to see very minor loss of information in terms of blurring the reconstructed image as JPEG compression is not fully lossless compression. Still, almost 50% size is reduced. You can verify this by checking size of files 'barbara_raw_input' and 'compressed_image' which are 262KB and 144KB respectively.



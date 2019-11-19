# Whi[s]kerMan Documentation

### Open WhiskerMan 

1. Extract Whiskerman  to a folder (e.g. "C:\user\whiskerman ").

2. Within Matlab, add this folder to the current path (e.g addpath(‘C:\ user\whiskerman’)).

3. Make a  folder that contains the horizontal view videos (e.g. cd ‘C:\ExampleH’). If doing 3D tracking also, within this folder, create a sub-folder called “vertical_view” (the name must be exactly this) containing the vertical view videos.

4. Within Matlab, “cd” to the folder containing the horizontal view videos and enter ‘WhiskerMan’ on the command line.  A GUI resembling the following should appear.



![open](https://github.com/PetersenLab/WhiskerMan/blob/master/doc/screenshots/Screenshot_open.png)

# 2D Whisker Tracking

## Tracking one whisker:

1. In the Whiskerman gui … Select “2D tracking”

2. Select the video to be tracked using the “Choose H video” menu.

3. Select “Continuous tracking” (otherwise only the current frame is tracked)

4. Press “Track whiskers” to start the tracking of the first frame. 

5. Whiskerman will prompt you to initialise the control points of the whiskers highlighting the position of the cursor (see the paper for explanation). 
To define the control points, select, using mouse clicks, the first, second and third control point of the Bezier curve in the horizontal view (starting with that closest to the snout).
If parameters are set appropriately for your video, Whiskerman will track each frame until either the end of the video is reached or the algorithm’s quality control criteria fail in a given frame.

6. To halt tracking at any point, deselect “Continuous tracking”.

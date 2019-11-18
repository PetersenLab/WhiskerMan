
## Tracking one whisker:

1. In the Whiskerman gui … Select “2D tracking”

2. Select the video to be tracked using the “Choose H video” menu.

3. Select “Continuous tracking” (otherwise only the current frame is tracked)

4. Press “Track whiskers” to start the tracking of the first frame. 

5. Whiskerman will prompt you to initialise the control points of the whiskers highlighting the position of the cursor (see the paper for explanation). 
To define the control points, select, using mouse clicks, the first, second and third control point of the Bezier curve in the horizontal view (starting with that closest to the snout).
If parameters are set appropriately for your video, Whiskerman will track each frame until either the end of the video is reached or the algorithm’s quality control criteria fail in a given frame.

6. To halt tracking at any point, deselect “Continuous tracking”.

**Common problem:** 

If Whiskerman tracks the first frame but then stops, a common reason for this is that the parameter “current whisker energy threshold” is too low.  See the paper for explanation of how Whiskerman evaluates the quality of its tracking solutions (using a cost function) and see the hints below for something to try that often works.

## Tracking multiple whiskers:

1. In the whiskerman gui … Select “2D tracking”

2. Select the video to be tracked using the “Choose H video” menu.

3. Define the whiskers to be tracked in the “identification of the whiskers” section. To do this, write the name of the first whisker and define the energy threshold for the whisker. 

4. To add another whisker, press “add whisker” and complete the information of the following whisker.

5. Press “Track whiskers” to start the tracking of the first frame. Whiskerman will prompt the user to initialise the control points of the whiskers highlighting the position of the cursor. To define the control points, select, using mouse clicks, the first, second and third control point of the Bezier curve in the horizontal view (starting with that closest to the snout). If more than one whisker were defined, follow the same procedure for the next whiskers.
Whiskerman will show the tracking of the first frame and will move to the next frame. To track the current frame, press “Track whiskers”. Alternatively, select “Continuous tracking” and then press “Track whiskers”, in which case the tracking proceeds automatically from frame to frame until: a) “Continuous tracking” is deselected, b) the end of the video is reached or c) tracking fails for all selected whiskers.

**Hints:**

- If the tracker halts after tracking just one frame, it may be that the “current whisker energy threshold” parameter is too low.  
Whiskerman evaluates the quality of its solution for a given frame by computing an “energy” value (see the paper for detailed explanation).  

- If this energy is less than threshold, tracking stops.  The radio button “whisker selected” will be deselected.  
To restart tracking, increase the value of the threshold, (re)select “whisker selected”, and then click “Track whiskers”. 

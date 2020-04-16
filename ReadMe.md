# Adapt Ride

This program takes the ride recordered data by SimRa App and generate a CSV file with the data that will be used by the trained Neural Network to find the incidents that may have been occurred during the ride

**Basic Configuration**

-*ridePath*: the path of the file that needs to be transformed

-*outputPath*: the path of the generated CSV file

**Optional Configuration**

-*avoidSeconds*: this variable sets the number of seconds that we want to avoid from the begining and end of the ride (if set to '0' the program takes all the ride). With this feature, initial and end readings will be omitted (they usually contain nosy data due to the manipulation of the phone).

**Advanced Configuration**

This 2 parameters must be set with the same values as the ones used to obtain the dataset (using AdaptDatabase) to train the neural network

-*WINDOWFRAME*: This defines the window length (in milliseconds) used to window the ride and obtain the statistics for each of the window

-*WINDOWSPLIT*: This defines in how many parts the window is splitted and therefore, the value of the overlapped window which is *WINDOWFRAME/WINDOWSPLIT*. For instance, WINDOWFRAME = 6000 and WINDOWSPLIT = 3 -> 1st window: 0-6 sec | 2nd window: 2-8 sec | ...

**Considerations**: The ride file must have at least one line in the incidents description header in order to indicate fields 'Bike' and 'pLoc'

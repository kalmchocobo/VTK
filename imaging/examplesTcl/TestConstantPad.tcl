catch {load vtktcl}
# Make an image larger by repeating the data.  Tile.

source vtkImageInclude.tcl

vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
#reader ReleaseDataFlagOff
#reader DebugOn

vtkImageConstantPad pad
pad SetInput [reader GetOutput]
pad SetOutputWholeExtent -300 355 -300 370 0 92
pad SetConstant 800
pad SetNumberOfThreads 1
pad ReleaseDataFlagOff

vtkImageViewer viewer
viewer SetInput [pad GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 1200
viewer SetColorLevel 600
#viewer DebugOn


# make interface
source WindowLevelInterface.tcl





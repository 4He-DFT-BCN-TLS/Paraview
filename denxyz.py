#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
denxyzvtk = LegacyVTKReader(FileNames=['/users/p1039/barranco/tmpdir/Rb5pP_3o2-relaxation/eta-0.1/Paraview/denxyz.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1724, 1080]

# show data in view
denxyzvtkDisplay = Show(denxyzvtk, renderView1)
# trace defaults for the display properties.
denxyzvtkDisplay.Representation = 'Outline'
denxyzvtkDisplay.ColorArrayName = ['POINTS', '']
denxyzvtkDisplay.ScalarOpacityUnitDistance = 0.7359123379451757

# reset view to fit data
renderView1.ResetCamera()

# hide data in view
Hide(denxyzvtk, renderView1)

# Properties modified on renderView1
renderView1.Background = [0.0, 0.0, 0.0]

# create a new 'Contour'
contour1 = Contour(Input=denxyzvtk)
contour1.ContourBy = ['POINTS', 'dens']
contour1.Isosurfaces = [0.1268286705017247]
contour1.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour1
contour1.Isosurfaces = [0.0109]

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.ColorArrayName = [None, '']

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# change solid color
contour1Display.DiffuseColor = [0.3333333333333333, 0.6666666666666666, 1.0]

# Properties modified on contour1Display
contour1Display.Specular = 0.01

# Properties modified on contour1Display
contour1Display.Specular = 0.04

# Properties modified on contour1Display
contour1Display.Specular = 0.05

# Properties modified on contour1Display
contour1Display.Specular = 0.07

# Properties modified on contour1Display
contour1Display.Specular = 0.08

# Properties modified on contour1Display
contour1Display.Specular = 0.1

# Properties modified on contour1Display
contour1Display.Specular = 0.11

# Properties modified on contour1Display
contour1Display.Specular = 0.14

# Properties modified on contour1Display
contour1Display.Specular = 0.15

# Properties modified on contour1Display
contour1Display.Specular = 0.17

# Properties modified on contour1Display
contour1Display.Specular = 0.2

# Properties modified on contour1Display
contour1Display.Specular = 0.22

# Properties modified on contour1Display
contour1Display.Specular = 0.24

# Properties modified on contour1Display
contour1Display.Specular = 0.25

# Properties modified on contour1Display
contour1Display.Specular = 0.28

# Properties modified on contour1Display
contour1Display.Specular = 0.3

# Properties modified on contour1Display
contour1Display.Specular = 0.33

# Properties modified on contour1Display
contour1Display.Specular = 0.36

# Properties modified on contour1Display
contour1Display.Specular = 0.38

# Properties modified on contour1Display
contour1Display.Specular = 0.39

# Properties modified on contour1Display
contour1Display.Specular = 0.44

# Properties modified on contour1Display
contour1Display.Specular = 0.48

# Properties modified on contour1Display
contour1Display.Specular = 0.58

# Properties modified on contour1Display
contour1Display.Specular = 0.71

# Properties modified on contour1Display
contour1Display.Specular = 0.84

# Properties modified on contour1Display
contour1Display.Specular = 0.88

# Properties modified on contour1Display
contour1Display.Specular = 1.0

# Properties modified on contour1Display
contour1Display.Specular = 0.5

# Properties modified on contour1Display
contour1Display.Specular = 0.75

# Properties modified on contour1Display
contour1Display.Opacity = 0.5

# Properties modified on contour1Display
contour1Display.Opacity = 0.75

# Properties modified on contour1Display
contour1Display.Opacity = 0.76

# Properties modified on contour1Display
contour1Display.Opacity = 0.77

# Properties modified on contour1Display
contour1Display.Opacity = 0.78

# Properties modified on contour1Display
contour1Display.Opacity = 0.79

# Properties modified on contour1Display
contour1Display.Opacity = 0.82

# Properties modified on contour1Display
contour1Display.Opacity = 0.85

# Properties modified on contour1Display
contour1Display.Opacity = 0.9

# Properties modified on contour1Display
contour1Display.Opacity = 0.95

# Properties modified on contour1Display
contour1Display.Opacity = 0.98

# Properties modified on contour1Display
contour1Display.Opacity = 1.0

# Properties modified on contour1Display
contour1Display.Opacity = 0.99

# Properties modified on contour1Display
contour1Display.Opacity = 0.98

# Properties modified on contour1Display
contour1Display.Opacity = 0.95

# Properties modified on contour1Display
contour1Display.Opacity = 0.93

# Properties modified on contour1Display
contour1Display.Opacity = 0.88

# Properties modified on contour1Display
contour1Display.Opacity = 0.87

# Properties modified on contour1Display
contour1Display.Opacity = 0.82

# Properties modified on contour1Display
contour1Display.Opacity = 0.81

# Properties modified on contour1Display
contour1Display.Opacity = 0.76

# Properties modified on contour1Display
contour1Display.Opacity = 0.75

# Properties modified on contour1Display
contour1Display.Opacity = 0.72

# Properties modified on contour1Display
contour1Display.Opacity = 0.71

# Properties modified on contour1Display
contour1Display.Opacity = 0.69

# Properties modified on contour1Display
contour1Display.Opacity = 0.62

# Properties modified on contour1Display
contour1Display.Opacity = 0.59

# Properties modified on contour1Display
contour1Display.Opacity = 0.54

# Properties modified on contour1Display
contour1Display.Opacity = 0.53

# Properties modified on contour1Display
contour1Display.Opacity = 0.5

# Properties modified on contour1Display
contour1Display.Opacity = 0.48

# Properties modified on contour1Display
contour1Display.Opacity = 0.44

# Properties modified on contour1Display
contour1Display.Opacity = 0.42

# Properties modified on contour1Display
contour1Display.Opacity = 0.38

# Properties modified on contour1Display
contour1Display.Opacity = 0.36

# Properties modified on contour1Display
contour1Display.Opacity = 0.34

# Properties modified on contour1Display
contour1Display.Opacity = 0.28

# Properties modified on contour1Display
contour1Display.Opacity = 0.26

# Properties modified on contour1Display
contour1Display.Opacity = 0.27

# Properties modified on contour1Display
contour1Display.Opacity = 0.36

# Properties modified on contour1Display
contour1Display.Opacity = 0.39

# Properties modified on contour1Display
contour1Display.Opacity = 0.48

# Properties modified on contour1Display
contour1Display.Opacity = 0.52

# Properties modified on contour1Display
contour1Display.Opacity = 0.55

# Properties modified on contour1Display
contour1Display.Opacity = 0.56

# Properties modified on contour1Display
contour1Display.Opacity = 0.62

# Properties modified on contour1Display
contour1Display.Opacity = 0.65

# Properties modified on contour1Display
contour1Display.Opacity = 0.68

# Properties modified on contour1Display
contour1Display.Opacity = 0.71

# Properties modified on contour1Display
contour1Display.Opacity = 0.74

# Properties modified on contour1Display
contour1Display.Opacity = 0.76

# Properties modified on contour1Display
contour1Display.Opacity = 0.77

# Properties modified on contour1Display
contour1Display.Opacity = 0.79

# Properties modified on contour1Display
contour1Display.Opacity = 0.81

# Properties modified on contour1Display
contour1Display.Opacity = 0.86

# Properties modified on contour1Display
contour1Display.Opacity = 0.88

# Properties modified on contour1Display
contour1Display.Opacity = 0.92

# Properties modified on contour1Display
contour1Display.Opacity = 0.93

# Properties modified on contour1Display
contour1Display.Opacity = 0.95

# Properties modified on contour1Display
contour1Display.Opacity = 0.97

# Properties modified on contour1Display
contour1Display.Opacity = 0.98

# Properties modified on contour1Display
contour1Display.Opacity = 0.97

# Properties modified on contour1Display
contour1Display.Opacity = 0.96

# Properties modified on contour1Display
contour1Display.Opacity = 0.95

# Properties modified on contour1Display
contour1Display.Opacity = 0.92

# Properties modified on contour1Display
contour1Display.Opacity = 0.9

# Properties modified on contour1Display
contour1Display.Opacity = 0.88

# Properties modified on contour1Display
contour1Display.Opacity = 0.87

# Properties modified on contour1Display
contour1Display.Opacity = 0.86

# Properties modified on contour1Display
contour1Display.Opacity = 0.85

# Properties modified on contour1Display
contour1Display.Opacity = 0.84

# Properties modified on contour1Display
contour1Display.Opacity = 0.78

# Properties modified on contour1Display
contour1Display.Opacity = 0.77

# Properties modified on contour1Display
contour1Display.Opacity = 0.75

# Properties modified on contour1Display
contour1Display.Opacity = 0.74

# Properties modified on contour1Display
contour1Display.Opacity = 0.73

# Properties modified on contour1Display
contour1Display.Opacity = 0.72

# Properties modified on contour1Display
contour1Display.Opacity = 0.62

# Properties modified on contour1Display
contour1Display.Opacity = 0.54

# Properties modified on contour1Display
contour1Display.Opacity = 0.51

# Properties modified on contour1Display
contour1Display.Opacity = 0.5

# Properties modified on contour1Display
contour1Display.Opacity = 0.48

# Properties modified on contour1Display
contour1Display.Opacity = 0.5

# Properties modified on contour1Display
contour1Display.Opacity = 0.55

# Properties modified on contour1Display
contour1Display.Opacity = 0.57

# Properties modified on contour1Display
contour1Display.Opacity = 0.6

# Properties modified on contour1Display
contour1Display.Opacity = 0.62

# Properties modified on contour1Display
contour1Display.Opacity = 0.65

# Properties modified on contour1Display
contour1Display.Opacity = 0.66

# Properties modified on contour1Display
contour1Display.Opacity = 0.65

# Properties modified on contour1Display
contour1Display.Opacity = 0.62

# Properties modified on contour1Display
contour1Display.Opacity = 0.6

# Properties modified on contour1Display
contour1Display.Opacity = 0.57

# Properties modified on contour1Display
contour1Display.Opacity = 0.55

# Properties modified on contour1Display
contour1Display.Opacity = 0.52

# Properties modified on contour1Display
contour1Display.Opacity = 0.51

# Properties modified on contour1Display
contour1Display.Opacity = 0.48

# Properties modified on contour1Display
contour1Display.Opacity = 0.46

# Properties modified on contour1Display
contour1Display.Opacity = 0.42

# Properties modified on contour1Display
contour1Display.Opacity = 0.4

# Properties modified on contour1Display
contour1Display.Opacity = 0.37

# Properties modified on contour1Display
contour1Display.Opacity = 0.35

# Properties modified on contour1Display
contour1Display.Opacity = 0.34

# Properties modified on contour1Display
contour1Display.Opacity = 0.35

# Properties modified on contour1Display
contour1Display.Opacity = 0.42

# Properties modified on contour1Display
contour1Display.Opacity = 0.45

# Properties modified on contour1Display
contour1Display.Opacity = 0.51

# Properties modified on contour1Display
contour1Display.Opacity = 0.52

# Properties modified on contour1Display
contour1Display.Opacity = 0.53

# Properties modified on contour1Display
contour1Display.Opacity = 0.54

# Properties modified on contour1Display
contour1Display.Opacity = 0.57

# Properties modified on contour1Display
contour1Display.Opacity = 0.59

# Properties modified on contour1Display
contour1Display.Opacity = 0.6

# Properties modified on contour1Display
contour1Display.Opacity = 0.62

# Properties modified on contour1Display
contour1Display.Opacity = 0.63

# Properties modified on contour1Display
contour1Display.Opacity = 0.65

# Properties modified on contour1Display
contour1Display.Opacity = 0.67

# Properties modified on contour1Display
contour1Display.Opacity = 0.71

# Properties modified on contour1Display
contour1Display.Opacity = 0.72

# Properties modified on contour1Display
contour1Display.Opacity = 0.74

# Properties modified on contour1Display
contour1Display.Opacity = 0.77

# Properties modified on contour1Display
contour1Display.Opacity = 0.79

# Properties modified on contour1Display
contour1Display.Opacity = 0.8

# Properties modified on contour1Display
contour1Display.Opacity = 0.81

# Properties modified on contour1Display
contour1Display.Opacity = 0.83

# Properties modified on contour1Display
contour1Display.Opacity = 0.85

# Properties modified on contour1Display
contour1Display.Opacity = 0.87

# Properties modified on contour1Display
contour1Display.Opacity = 0.86

# Properties modified on contour1Display
contour1Display.Opacity = 0.82

# Properties modified on contour1Display
contour1Display.Opacity = 0.76

# Properties modified on contour1Display
contour1Display.Opacity = 0.69

# Properties modified on contour1Display
contour1Display.Opacity = 0.67

# Properties modified on contour1Display
contour1Display.Opacity = 0.6

# Properties modified on contour1Display
contour1Display.Opacity = 0.58

# Properties modified on contour1Display
contour1Display.Opacity = 0.5

# Properties modified on contour1Display
contour1Display.Opacity = 0.43

# Properties modified on contour1Display
contour1Display.Opacity = 0.41

# Properties modified on contour1Display
contour1Display.Opacity = 0.35

# Properties modified on contour1Display
contour1Display.Opacity = 0.34

# Properties modified on contour1Display
contour1Display.Opacity = 0.28

# Properties modified on contour1Display
contour1Display.Opacity = 0.26

# Properties modified on contour1Display
contour1Display.Opacity = 0.23

# Properties modified on contour1Display
contour1Display.Opacity = 0.22

# Properties modified on contour1Display
contour1Display.Opacity = 0.23

# Properties modified on contour1Display
contour1Display.Opacity = 0.28

# Properties modified on contour1Display
contour1Display.Opacity = 0.29

# Properties modified on contour1Display
contour1Display.Opacity = 0.31

# Properties modified on contour1Display
contour1Display.Opacity = 0.32

# Properties modified on contour1Display
contour1Display.Opacity = 0.34

# Properties modified on contour1Display
contour1Display.Opacity = 0.35

# Properties modified on contour1Display
contour1Display.Opacity = 0.36

# Properties modified on contour1Display
contour1Display.Opacity = 0.38

# Properties modified on contour1Display
contour1Display.Opacity = 0.39

# Properties modified on contour1Display
contour1Display.Opacity = 0.41

# Properties modified on contour1Display
contour1Display.Opacity = 0.5

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on contour1Display
contour1Display.Specular = 0.5

# Properties modified on contour1Display
contour1Display.Specular = 0.75

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Sphere'
sphere1 = Sphere()

# Properties modified on sphere1
sphere1.Center = [11.99, -12.02, 38.73]
sphere1.Radius = 2.0
sphere1.ThetaResolution = 256
sphere1.PhiResolution = 256

# show data in view
sphere1Display = Show(sphere1, renderView1)
# trace defaults for the display properties.
sphere1Display.ColorArrayName = [None, '']

# Properties modified on sphere1
sphere1.Radius = 1.5

# change solid color
sphere1Display.DiffuseColor = [0.6666666666666666, 1.0, 0.0]

# Properties modified on sphere1Display
sphere1Display.Specular = 0.75

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data bounds
renderView1.ResetCamera(10.4900283813, 13.4899711609, -13.5199718475, -10.5200281143, 37.2299995422, 40.2299995422)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Text'
text1 = Text()

# Properties modified on text1
text1.Text = ' t = 150 ps'

# show data in view
text1Display = Show(text1, renderView1)

# Properties modified on text1Display
text1Display.Position = [0.05, 0.0]

# Properties modified on text1Display
text1Display.Position = [0.05, 0.9]

# current camera placement for renderView1
renderView1.CameraPosition = [0.06381607055664064, -180.98031727344997, 21.138175679438508]
renderView1.CameraFocalPoint = [0.06381607055664064, 0.4417545066442732, 5.26580108529381]
renderView1.CameraViewUp = [0.0, 0.08715574274765818, 0.9961946980917455]
renderView1.CameraParallelScale = 47.13484971220008

# save screenshot
SaveScreenshot('/users/p1039/barranco/tmpdir/Rb5pP_3o2-relaxation/eta-0.1/Paraview/test.png', magnification=2, quality=100, view=renderView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.06381607055664064, -180.98031727344997, 21.138175679438508]
renderView1.CameraFocalPoint = [0.06381607055664064, 0.4417545066442732, 5.26580108529381]
renderView1.CameraViewUp = [0.0, 0.08715574274765818, 0.9961946980917455]
renderView1.CameraParallelScale = 47.13484971220008

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
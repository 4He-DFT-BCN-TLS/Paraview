#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2150, 1366]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 0
renderView1.Background = [0.0, 0.0, 0.0]

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(0, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Legacy VTK Reader'
denxyzvtk = LegacyVTKReader(FileNames=['denxyz.vtk'])

# show data in view
denxyzvtkDisplay = Show(denxyzvtk, renderView1)
# trace defaults for the display properties.
denxyzvtkDisplay.Representation = 'Outline'
denxyzvtkDisplay.ColorArrayName = ['POINTS', '']
denxyzvtkDisplay.OSPRayScaleArray = 'dens'
denxyzvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
denxyzvtkDisplay.SelectOrientationVectors = 'None'
denxyzvtkDisplay.ScaleFactor = 12.75999984741211
denxyzvtkDisplay.SelectScaleArray = 'dens'
denxyzvtkDisplay.GlyphType = 'Arrow'
denxyzvtkDisplay.GlyphTableIndexArray = 'dens'
denxyzvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
denxyzvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
denxyzvtkDisplay.ScalarOpacityUnitDistance = 0.6928203147425989

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(denxyzvtk, renderView1)

# create a new 'Contour'
contour1 = Contour(Input=denxyzvtk)
contour1.ContourBy = ['POINTS', 'dens']
contour1.Isosurfaces = [0.0109]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 7.872348403930665
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'
contour1Display.GaussianRadius = 3.9361742019653323
contour1Display.SetScaleArray = [None, '']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = [None, '']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
contour1Display.SetRepresentationType('Wireframe')

# Properties modified on contour1Display
contour1Display.Opacity = 0.3

# create a new 'Text'
text1 = Text()

# Properties modified on text1
text1.Text = 't = 123 ps'

# show data in view
text1Display = Show(text1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text1Display
text1Display.FontFamily = 'Times'

# Properties modified on text1Display
text1Display.FontSize = 18


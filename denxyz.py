from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ViewSize = [1842, 1080]
renderView1.Background = [0.0, 0.0, 0.0]
denxyzvtk = LegacyVTKReader(FileNames=['denxyz.vtk'])
denxyzvtkDisplay = Show(denxyzvtk, renderView1)
denxyzvtkDisplay.Representation = 'Outline'
denxyzvtkDisplay.ColorArrayName = ['POINTS', '']
denxyzvtkDisplay.ScalarOpacityUnitDistance = 0.7028044359875067
renderView1.ResetCamera()
Hide(denxyzvtk, renderView1)
contour1 = Contour(Input=denxyzvtk)
contour1.ContourBy = ['POINTS', 'dens']
contour1.Isosurfaces = [0.013057181611657143]
contour1.PointMergeMethod = 'Uniform Binning'
contour1.Isosurfaces = [0.0109]
contour1Display = Show(contour1, renderView1)
contour1Display.ColorArrayName = [None, '']
renderView1.ResetCamera()
contour1Display.DiffuseColor = [0.3333333333333333, 0.6666666666666666, 1.0]
contour1Display.Opacity = 0.45
contour1Display.Specular = 0.5
sphere1 = Sphere()
sphere1.Center = [0.0, 0.0, -43.2]
sphere1.Radius = 2.0
sphere1.ThetaResolution = 128
sphere1.PhiResolution = 128
sphere1Display = Show(sphere1, renderView1)
sphere1Display.ColorArrayName = [None, '']
renderView1.ResetCamera()
renderView1.ResetCamera()
renderView1.ResetCamera()
renderView1.ResetCamera()
renderView1.ResetCamera()
renderView1.ResetCamera()
renderView1.ResetCamera()
SetActiveSource(contour1)
renderView1.ResetCamera(-19.4281272888, 19.4281272888, -23.030090332, 23.030090332, -34.2293167114, 11.8216047287)
SetActiveSource(sphere1)
sphere1Display.Specular = 0.5
sphere1Display.DiffuseColor = [0.6666666666666666, 1.0, 0.0]
text1 = Text()
text1.Text = 't = 0 ps'
text1Display = Show(text1, renderView1)
text1Display.Position = [0.05, 0.0]
text1Display.Position = [0.05, 0.9]
renderView1.CameraPosition = [-1.1901590823981678e-13, -121.08747024551023, -11.203855991363525]
renderView1.CameraFocalPoint = [0.0, 0.0, -11.203855991363525]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 37.92108954161494
SaveScreenshot('denxyz.png', magnification=1, quality=100, view=renderView1)
renderView1.CameraPosition = [-1.1901590823981678e-13, -121.08747024551023, -11.203855991363525]
renderView1.CameraFocalPoint = [0.0, 0.0, -11.203855991363525]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 37.92108954161494

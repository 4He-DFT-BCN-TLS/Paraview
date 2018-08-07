# update the view to ensure updated data information
renderView1.Update()

# current camera placement for renderView1
renderView1.CameraPosition = [-0.007884979248046875, -175.0, 0.0009860992431640625]
renderView1.CameraFocalPoint = [-0.007884979248046875, -0.0003490447998046875, 0.0009860992431640625]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 67.99687199981692

# save screenshot
SaveScreenshot('denxyz.png', renderView1, ImageResolution=[2833, 1800])

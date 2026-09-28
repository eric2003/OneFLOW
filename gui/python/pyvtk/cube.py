import vtk

cube = vtk.vtkCubeSource()
cube.Update()

cube_mapper = vtk.vtkPolyDataMapper()
cube_mapper.SetInputData(cube.GetOutput())

cube_actor = vtk.vtkActor()
cube_actor.SetMapper(cube_mapper)

cube_actor.GetProperty().SetColor(1.0, 0.0, 0.0)

renderer = vtk.vtkRenderer()
renderer.SetBackground(0.0, 0.0, 1.0)
renderer.AddActor(cube_actor)

render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window.Render()
render_window.SetSize(400,400)
render_window.SetWindowName("OneFLOW CFD VTK Cube Test")

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)
interactor.Initialize()
interactor.Start()
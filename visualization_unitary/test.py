import vtk

# Create a renderer
renderer = vtk.vtkRenderer()

# Create a render window
render_window = vtk.vtkRenderWindow()
render_window.SetSize(1080, 720)  # Resized the window
render_window.AddRenderer(renderer)

# Create a render window interactor
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

# Create a list of colors for the actors
colors = [
    (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0),
    (1, 0, 1), (0, 1, 1), (1, 0.5, 0), (0.5, 1, 0),
    (0, 1, 0.5), (0.5, 0, 1), (1, 1, 1), (0, 0, 0),
    (0.5, 0.5, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.2, 0.2, 0.2)
]

# Create actors in a loop
for i, color in enumerate(colors):
    # Create a regular polygon source
    polygon_source = vtk.vtkRegularPolygonSource()
    polygon_source.SetNumberOfSides(6)  # Set number of sides (hexagon)
    polygon_source.SetRadius(1)

    # Create a mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(polygon_source.GetOutputPort())

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)  # Set actor color

    # Position the actors in the scene
    x = i % 4
    y = (i // 4) % 4
    z = i // 16
    actor.SetPosition(x * 3, y * 3, z * 3)

    # Add actor to the renderer
    renderer.AddActor(actor)

# Set background color
renderer.SetBackground(0.1, 0.2, 0.4)

# Render
render_window.Render()

# Start the interaction
render_window_interactor.Start()
